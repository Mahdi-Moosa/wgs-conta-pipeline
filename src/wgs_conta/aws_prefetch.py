from __future__ import annotations
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, Optional

from .acquisition import PrefetchedSRA,prefetch_aws_sra
from .util import eprint

@dataclass
class AwsReady:
    prefetched: PrefetchedSRA
    assembly_future: object = None

class AwsPrefetchManager:
    """Bounded producer for validated AWS SRA objects.

    Download workers are intentionally independent of analysis workers.  The
    manager keeps at most ``prefetch_depth`` accession futures staged ahead,
    while ``download_workers`` controls actual concurrent S3 transfers.
    """
    def __init__(self,s,accessions:Iterable[str],assembly_submitter:Optional[Callable]=None,download_func=None):
        self.s=s
        self.order=list(dict.fromkeys(accessions))
        self.assembly_submitter=assembly_submitter
        self.download_func=download_func or prefetch_aws_sra
        self.executor=ThreadPoolExecutor(max_workers=max(1,int(s.aws_download_workers)),thread_name_prefix='aws-download')
        self.lock=threading.Lock()
        self.futures:Dict[str,object]={}
        self.consumed=set()
        self.cursor=0
        self.depth=max(int(s.aws_download_workers),int(s.aws_prefetch_depth))
        self._fill_locked()

    def _task(self,acc):
        holder={'assembly_future':None}
        def metadata_hook(md):
            if holder['assembly_future'] is None and self.assembly_submitter is not None:
                holder['assembly_future']=self.assembly_submitter(acc,dict(md))
        pf=self.download_func(
            acc,self.s.sra_root,self.s.aws_cli,self.s.aws_bucket,self.s.retries,
            metadata_hook=metadata_hook
        )
        return AwsReady(pf,holder['assembly_future'])

    def _schedule_locked(self,acc):
        if acc in self.futures or acc in self.consumed: return
        self.futures[acc]=self.executor.submit(self._task,acc)
        eprint(f'+ [AWS prefetch queue] {acc}')

    def _fill_locked(self):
        active=sum(1 for a in self.futures if a not in self.consumed)
        while active<self.depth and self.cursor<len(self.order):
            acc=self.order[self.cursor]; self.cursor+=1
            if acc in self.futures or acc in self.consumed: continue
            self._schedule_locked(acc); active+=1

    def get(self,acc):
        with self.lock:
            if acc not in self.futures:
                # A fallback/reprioritization can cause the DB to hand us an
                # accession outside the current look-ahead window. Schedule it
                # immediately rather than waiting for the cursor to catch up.
                self._schedule_locked(acc)
            fut=self.futures[acc]
        try:
            return fut.result()
        finally:
            with self.lock:
                self.consumed.add(acc)
                self.futures.pop(acc,None)
                self._fill_locked()

    def close(self):
        # The bounded queue prevents an unbounded tail of speculative transfers.
        # Active S3 subprocesses are allowed to finish cleanly; queued futures are
        # cancelled once the controller has no more scientific jobs.
        self.executor.shutdown(wait=True,cancel_futures=True)
