from __future__ import annotations

import argparse
import configparser
import fcntl
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path

from . import db
from .aws_prefetch import AwsPrefetchManager
from .config import load_config
from .master import MasterStore,parse_registry
from .util import eprint,ensure_executable
from .worker import _prefetch_exact_assembly,finalize_committed,worker_loop

SCIENTIFIC_METHOD_VERSION='2.3.0'


from .sra_cache import iter_deterministic_sra_files

def _version(root):
    return (Path(__file__).with_name('VERSION')).read_text().strip()


def _available_sources(mode):
    if mode=='aws': return ['aws']
    if mode=='ncbi': return ['ncbi']
    if mode=='ena': return ['ena']
    if mode=='both': return ['ncbi','ena']
    if mode=='all': return ['aws','ncbi','ena']
    raise ValueError(f'Unsupported source mode: {mode}')


def _cached_sra_source_map(root):
    """Return accession -> preferred local normalized-SRA backend.

    AWS cache objects live under ``aws_cache``. Other local SRA files retain the
    historical NCBI/local preference so v2.4.3 does not relabel existing files.
    """
    root=Path(root); out={}
    if not root.exists(): return out
    for p in iter_deterministic_sra_files(root):
        if not p.is_file(): continue
        acc=p.stem
        try: rel=p.relative_to(root)
        except ValueError: rel=p
        src='aws' if rel.parts and rel.parts[0]=='aws_cache' else 'ncbi'
        out.setdefault(acc,src)
    return out


def _selection(s,registry,store):
    available=_available_sources(s.source_mode)
    cache=_cached_sra_source_map(s.sra_root)
    pending=[]; skipped=[]
    for r in registry:
        complete,reason,_=store.completion_detail(r['accession'])
        if complete:
            skipped.append((r,'COMPLETE')); continue
        pref=r['preferred_source'] if r['preferred_source'] in available else ''
        cached=cache.get(r['accession'])
        if cached in available:
            pref=cached
        elif not pref and 'aws' in available:
            # AWS is the primary bulk acquisition path in v2.4.3. NCBI/ENA
            # remain immediate fallbacks after the AWS retry budget exhausts.
            pref='aws'
        pending.append((r,pref,reason))
    return available,pending,skipped


def main(argv=None):
    ap=argparse.ArgumentParser(description='WGS contamination: validated strict calls, acquisition and exact-assembly challenge')
    ap.add_argument('--accessions',type=Path,default=None,help='Master accession registry TSV; defaults to <project_root>/master/accessions.tsv')
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--source',choices=['aws','ncbi','ena','ebi','both','all'])
    ap.add_argument('--state-db',type=Path,default=None)
    ap.add_argument('--retry-errors',action='store_true',help='Deprecated compatibility flag; scientifically incomplete error states retry automatically')
    ap.add_argument('--status',action='store_true')
    ap.add_argument('--dry-run',action='store_true',help='Show master skip/queue decisions only; do not modify state or run tools')
    a=ap.parse_args(argv)
    from .preflight import validate_runtime,guard_resume_identity
    s=load_config(a.config)
    validate_runtime(s, require_resources=not a.dry_run, config_path=a.config)
    import os
    cfg=configparser.ConfigParser();cfg.read(a.config)
    os.environ.setdefault('WGS_CONTA_GLOBAL_ANALYSIS_SLOTS',cfg.get('workers','global_analysis_slots',fallback='1'))
    from . import worker
    worker._GLOBAL_ANALYSIS_SLOTS=worker._configured_global_analysis_slots()
    worker._GLOBAL_ANALYSIS_SEMAPHORE=threading.BoundedSemaphore(worker._GLOBAL_ANALYSIS_SLOTS)
    if a.source: s=replace(s,source_mode='ena' if a.source=='ebi' else a.source)
    if a.state_db: s=replace(s,state_db=a.state_db)

    guard_resume_identity(s,a.config,dry_run=a.dry_run)
    pipeline_version=_version(Path(__file__).resolve().parents[1])
    registry_path=a.accessions or (s.project_root/'master'/'accessions.tsv')
    # Keep deterministic cache lookup bound to this invocation, never a research default.
    os.environ['WGS_CONTA_CONFIG']=str(a.config.resolve())
    os.environ['WGS_CONTA_REGISTRY']=str(registry_path.resolve())
    store=MasterStore(s.project_root/'master',pipeline_version,SCIENTIFIC_METHOD_VERSION)
    registry=parse_registry(registry_path)
    import re
    if any(not re.fullmatch(r"(?:SRR|ERR|DRR)[0-9]+",r["accession"]) for r in registry):
        raise ValueError("registry requires atomic SRA/ENA run accessions")
    available,pending,skipped=_selection(s,registry,store)

    if a.dry_run:
        print(f'Pipeline_version\t{pipeline_version}')
        print(f'Scientific_method_version\t{SCIENTIFIC_METHOD_VERSION}')
        print(f'Active_accessions\t{len(registry)}')
        print(f'Already_complete\t{len(skipped)}')
        print(f'Would_queue\t{len(pending)}')
        for r,pref,reason in pending:
            print(f"QUEUE\t{r['accession']}\tpreferred={pref or 'AUTO'}\tpriority={r['priority']}\treason={reason}")
        for r,_ in skipped:
            print(f"SKIP\t{r['accession']}\tcompatible_method={SCIENTIFIC_METHOD_VERSION}")
        return

    tools=['seqkit','minimap2','samtools','kraken2','datasets','esearch','esummary','xtract','curl']
    if any(x in available for x in {'aws','ncbi'}): tools+=['fasterq-dump','vdb-dump','vdb-validate']
    if 'ncbi' in available: tools+=['prefetch','efetch']
    for t in sorted(set(tools)): ensure_executable(t)
    if 'aws' in available:
        if not s.aws_cli.exists(): raise SystemExit(f'AWS CLI not found: {s.aws_cli}')
        if not s.aws_cli.is_file(): raise SystemExit(f'AWS CLI path is not a file: {s.aws_cli}')
    if 'ena' in available and s.ena_transport in {'aspera','auto'}:
        if s.ena_transport=='aspera' and not s.ena_ascp.exists(): raise SystemExit(f'Aspera ascp not found: {s.ena_ascp}')
        if s.ena_transport=='aspera' and not s.ena_ascp_key.exists(): raise SystemExit(f'Aspera key not found: {s.ena_ascp_key}')
    if not s.human_index.exists(): raise SystemExit(f'Human index not found: {s.human_index}')
    if not s.kraken_db.exists(): raise SystemExit(f'Kraken DB not found: {s.kraken_db}')
    for p in [s.results_root,s.evidence_root,s.logs_root,s.work_root,s.assembly_root,store.root,s.sra_root]: p.mkdir(parents=True,exist_ok=True)

    lock_path=Path(str(s.state_db)+'.controller.lock'); lock_path.parent.mkdir(parents=True,exist_ok=True); lock_fh=lock_path.open('w')
    try: fcntl.flock(lock_fh,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise SystemExit(f'Another controller is already using state DB: {s.state_db}')

    db.init_db(s.state_db)
    if a.status:
        print(db.status_counts(s.state_db,selected_only=True)); return

    reg_by_acc={r['accession']:r for r in registry}
    rows=[(r['accession'],r['label'],pref,r['priority']) for r,pref,_ in pending]
    for r,pref,reason in pending:
        eprint(f"MASTER QUEUE {r['accession']}: {reason}; preferred={pref or 'AUTO'} priority={r['priority']}")

    db.prepare_master_jobs(s.state_db,rows,pipeline_version)
    if a.retry_errors:
        eprint('NOTE: --retry-errors is no longer required; incomplete runtime error states reconcile automatically')
    finalize_committed(s,master_store=store,registry_by_accession=reg_by_acc)

    # Shared assembly pool: independent from read acquisition/analysis workers.
    assembly_executor=ThreadPoolExecutor(
        max_workers=max(1,s.assembly_prefetch_workers),thread_name_prefix='assembly-prefetch'
    )
    def assembly_submitter(acc,md):
        biosample=(md.get('BioSample') or '').strip()
        if not biosample: return None
        return assembly_executor.submit(_prefetch_exact_assembly,s,acc,md)

    aws_manager=None
    if 'aws' in available:
        aws_order=[r['accession'] for r,pref,_ in pending if pref=='aws']
        aws_manager=AwsPrefetchManager(s,aws_order,assembly_submitter=assembly_submitter)

    specs=[]
    if 'aws' in available:
        for i in range(s.aws_analysis_workers): specs.append(('aws',f'aws-{i+1}'))
    if 'ncbi' in available:
        for i in range(s.ncbi_workers): specs.append(('ncbi',f'ncbi-{i+1}'))
    if 'ena' in available:
        for i in range(s.ena_workers): specs.append(('ena',f'ena-{i+1}'))
    available_set={b for b,_ in specs}

    eprint(f'Master registry: {registry_path}')
    eprint(f'Pipeline version: {pipeline_version}; scientific method: {SCIENTIFIC_METHOD_VERSION}')
    eprint(f'Active registry accessions: {len(registry)}; already complete: {len(skipped)}; queued: {len(pending)}')
    eprint(f'Source mode: {s.source_mode}; workers: {specs}')
    if aws_manager is not None:
        eprint(f'AWS pipeline: downloads={s.aws_download_workers}; analysis={s.aws_analysis_workers}; prefetch_depth={s.aws_prefetch_depth}')
    eprint(f'Assembly prefetch workers: {s.assembly_prefetch_workers}; persistent cache: {s.assembly_root}')
    eprint(f'Runtime state DB: {s.state_db}')
    eprint(f'Master results: {store.results}')

    ts=[]
    try:
        for backend,name in specs:
            t=threading.Thread(
                target=worker_loop,
                args=(s,backend,name,available_set,store,reg_by_acc,assembly_executor,aws_manager),
                name=name
            )
            t.start(); ts.append(t)
        for t in ts: t.join()
    finally:
        if aws_manager is not None: aws_manager.close()
        assembly_executor.shutdown(wait=True,cancel_futures=True)

    finalize_committed(s,master_store=store,registry_by_accession=reg_by_acc)
    eprint('Final runtime status:',db.status_counts(s.state_db,selected_only=True))
    complete=sum(1 for r in registry if store.is_complete(r['accession']))
    eprint(f'Master complete: {complete}/{len(registry)} active accessions')

if __name__=='__main__': main()
