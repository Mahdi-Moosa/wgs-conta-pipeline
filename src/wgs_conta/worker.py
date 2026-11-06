from __future__ import annotations
import configparser,csv,gzip,json,os,shutil,threading,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from . import db
from .acquisition import UnsupportedRunError,acquire_aws_prefetched,acquire_ena,acquire_ncbi,load_manifest
from .analysis import bacterial_challenge,extract_strict,human_map,kraken_strict,seqkit_pair_stats,strict_human_from_bam
from .ion_analysis import bacterial_challenge_single,extract_strict_single,human_map_single,kraken_strict_single,seqkit_single_stats,strict_human_single_from_bam
from .assemblies import download_assembly,resolve_exact_assembly
from .util import eprint,now_iso

PAIRED_PLATFORMS={'ILLUMINA','BGISEQ','DNBSEQ'}
ION_PLATFORMS={'ION_TORRENT','IONTORRENT','ION TORRENT'}

# v2.4.8 operational safety:
# acquisition/download and assembly prefetch remain concurrent, but the
# memory-heavy scientific analysis path is globally bounded across every
# backend worker in this master process.
#
# Default=1 is appropriate for the current 32-GB workstation after observed
# Kraken/minimap2 OOM events.  It may be overridden explicitly at launch:
# WGS_CONTA_GLOBAL_ANALYSIS_SLOTS=2 ./run_global.sh
def _configured_global_analysis_slots():
    """
    Precedence:
      1. WGS_CONTA_GLOBAL_ANALYSIS_SLOTS environment variable
      2. [workers] global_analysis_slots in config.ini
      3. hard fallback = 2
    """
    raw=os.environ.get("WGS_CONTA_GLOBAL_ANALYSIS_SLOTS")

    if raw is None:
        cfg=configparser.ConfigParser()
        cfg_path=os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "config.ini"
        )
        cfg.read(cfg_path)
        raw=cfg.get(
            "workers",
            "global_analysis_slots",
            fallback="2"
        )

    try:
        value=int(raw)
    except (TypeError,ValueError):
        raise RuntimeError(
            "global analysis slots must be an integer >=1"
        )

    if value < 1:
        raise RuntimeError(
            "global analysis slots must be >=1"
        )

    return value


_GLOBAL_ANALYSIS_SLOTS=(
    _configured_global_analysis_slots()
)

_GLOBAL_ANALYSIS_SEMAPHORE=threading.BoundedSemaphore(
    _GLOBAL_ANALYSIS_SLOTS
)

def norm(md):
    return {'BioProject':md.get('BioProject',''),'BioSample':md.get('BioSample',''),'Organism':md.get('ScientificName','') or md.get('Organism',''),'Platform':(md.get('Platform','') or '').upper().replace('SRA_PLATFORM_',''),'Layout':(md.get('LibraryLayout','') or '').upper(),'Metadata_source':md.get('_metadata_source','UNKNOWN')}


def _sra_reject_value(acquired):
    return acquired.reject_pairs_removed if acquired.read_filter_policy=='EXCLUDE_REJECT_KEEP_PASS_CRITERIA' else 'NA'



def _prefetch_exact_assembly(s,acc,metadata):
    """Speculatively resolve/download exact assembly from metadata.

    This optimization must never weaken the scientific path: callers only use
    this result for a positive sample if it completed successfully. Any failure
    is logged and the normal blocking resolver is retried in the analysis path.
    """
    md=norm(metadata)
    biosample=md.get('BioSample','')
    if not biosample:
        return None,None
    eprint(f'+ [assembly prefetch] START {acc}/{biosample}')
    res=resolve_exact_assembly(acc,biosample,s.curated_run_assembly,s.assembly_query_retries,s.ena_search_api)
    if not res.lookup_complete:
        raise RuntimeError(f'prefetch exact-assembly lookup incomplete: {res.lookup_status}')
    if not res.accession:
        eprint(f'+ [assembly prefetch] NO_EXACT_ASSEMBLY {acc}/{biosample}')
        return res,None
    bf=download_assembly(res,s.assembly_root,s.curated_assembly_paths,s.ena_ascp,s.ena_ascp_key,s.ena_ascp_port,s.ena_ascp_rate,s.retries)
    eprint(f'+ [assembly prefetch] READY {acc}: {res.accession}')
    return res,bf

def _assembly_for_positive(s,acc,biosample,assembly_future=None):
    """Return (resolution, fasta, error); assembly QC is non-fatal."""
    if assembly_future is not None:
        try:
            res,bf=assembly_future.result()
            if res is not None:
                return res,bf,None
        except Exception as exc:
            eprint(f'WARNING: speculative assembly prefetch failed for {acc}: {exc}; retrying synchronously')
    try:
        res=resolve_exact_assembly(acc,biosample,s.curated_run_assembly,s.assembly_query_retries,s.ena_search_api)
    except Exception as exc:
        return None,None,exc
    if not res.lookup_complete:
        return res,None,RuntimeError(
            f'Exact-assembly lookup incomplete for {acc}/{biosample}: {res.lookup_status}'
        )
    if not res.accession:
        return res,None,None
    try:
        bf=download_assembly(
            res,s.assembly_root,s.curated_assembly_paths,
            s.ena_ascp,s.ena_ascp_key,s.ena_ascp_port,
            s.ena_ascp_rate,s.retries
        )
        return res,bf,None
    except Exception as exc:
        return res,None,exc

def _record_assembly_failure(evidence,acc,exc):
    (evidence/'assembly_check_error.txt').write_text(str(exc)+'\n')
    eprint(f'WARNING: exact bacterial challenge unavailable for {acc}: {exc}')

def gzip_candidate(path):
    path=Path(path)
    if not path.exists(): return
    gz=Path(str(path)+'.gz')
    with path.open('rb') as src,gzip.open(gz,'wb',6) as dst: shutil.copyfileobj(src,dst)
    path.unlink()

def acquire_one(s,backend,acc,work,metadata_hook=None):
    # The fallback keeps internal/third-party acquisition shims written against
    # pre-v2.4.2 signatures working; production acquire_* functions accept the
    # metadata hook directly.
    try:
        return acquire_ncbi(acc,work,s.sra_root,s.threads_per_worker,s.retries,metadata_hook=metadata_hook) if backend=='ncbi' else acquire_ena(acc,work,s.ena_api,s.retries,s.ena_verify_md5,s.ena_transport,s.ena_ascp,s.ena_ascp_key,s.ena_ascp_port,s.ena_ascp_rate,s.ena_https_fallback,metadata_hook=metadata_hook)
    except TypeError as exc:
        if 'metadata_hook' not in str(exc): raise
        return acquire_ncbi(acc,work,s.sra_root,s.threads_per_worker,s.retries) if backend=='ncbi' else acquire_ena(acc,work,s.ena_api,s.retries,s.ena_verify_md5,s.ena_transport,s.ena_ascp,s.ena_ascp_key,s.ena_ascp_port,s.ena_ascp_rate,s.ena_https_fallback)

def _call_acquire_one(s,backend,acc,work,metadata_hook=None):
    try:
        return acquire_one(s,backend,acc,work,metadata_hook=metadata_hook)
    except TypeError as exc:
        if 'metadata_hook' not in str(exc): raise
        return acquire_one(s,backend,acc,work)

def _call_analyze_one(s,backend,acc,label,work,acquired,t0,assembly_future=None):
    try:
        return analyze_one(s,backend,acc,label,work,acquired,t0,assembly_future=assembly_future)
    except TypeError as exc:
        if 'assembly_future' not in str(exc): raise
        return analyze_one(s,backend,acc,label,work,acquired,t0)

def analyze_paired(s,backend,acc,label,work,acquired,t0,assembly_future=None):
    md=norm(acquired.metadata)
    if md['Platform'] not in PAIRED_PLATFORMS: raise RuntimeError(f"Unsupported paired platform {md['Platform']}")
    if md['Layout']!='PAIRED' or acquired.r2 is None: raise RuntimeError(f"Unsupported paired layout {md['Layout']}; expected PAIRED with two read files")

    # Evidence may be regenerated on an analysis retry; downloaded reads are not touched.
    evidence=s.evidence_root/acc; shutil.rmtree(evidence,ignore_errors=True); evidence.mkdir(parents=True,exist_ok=True)
    stats=seqkit_pair_stats(acquired.r1,acquired.r2)
    bam=work/f'{acc}.human.bam'; human_map(s.human_index,acquired.r1,acquired.r2,bam,evidence/'human_minimap2.log',s.threads_per_worker,s.sort_threads)
    calls=strict_human_from_bam(bam,evidence); n=int(calls['unique'])
    kh=0; assembly=''; assembly_source='NA'; lookup_status='NOT_NEEDED_NEGATIVE'; acheck='NOT_NEEDED_NEGATIVE'; bact={'ge1_95':0,'both_95':0,'ge1_98':0,'both_98':0}
    if n:
        sr1,sr2=extract_strict(calls['qnames'],acquired.r1,acquired.r2,evidence)
        kh,_=kraken_strict(s.kraken_db,sr1,sr2,evidence,s.kraken_threads)
        res,bf,aerr=_assembly_for_positive(s,acc,md['BioSample'],assembly_future)
        if res is not None:
            assembly=res.accession
            assembly_source=res.source
            lookup_status=res.lookup_status
        else:
            lookup_status='ASSEMBLY_LOOKUP_FAILED'
        if aerr is not None:
            acheck='FAILED'
            _record_assembly_failure(evidence,acc,aerr)
        elif assembly:
            try:
                bact=bacterial_challenge(
                    bf,sr1,sr2,evidence,
                    s.threads_per_worker,s.sort_threads
                )
                acheck='PASS' if bact['ge1_95']==0 else 'BACTERIAL_MATCH_PRESENT'
            except Exception as exc:
                acheck='FAILED'
                _record_assembly_failure(evidence,acc,exc)
        else:
            acheck='NO_EXACT_ASSEMBLY'
        gzip_candidate(sr1); gzip_candidate(sr2)

    with (evidence/'run_metadata.tsv').open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t')
        w.writerow(['Accession','Label','Analysis_branch','Quantification_unit','Acquisition_source','Acquisition_format','Read_filter_policy','SRA_reject_units_removed','SRA_reject_pairs_removed','SRA_reject_reads_removed','BioProject','BioSample','Organism','Platform','Layout','Assembly','Assembly_source','Assembly_lookup_status'])
        rej=_sra_reject_value(acquired)
        w.writerow([acc,label,'PAIRED_SHORTREAD','PAIRS',acquired.source,acquired.fmt,acquired.read_filter_policy,rej,rej,(2*rej if isinstance(rej,int) else 'NA'),md['BioProject'],md['BioSample'],md['Organism'],md['Platform'],md['Layout'],assembly,assembly_source,lookup_status])
    rate=n/stats['pairs']*1_000_000 if stats['pairs'] else 0; kfrac=kh/n if n else 0; bc_ok=acheck in {'PASS','BACTERIAL_MATCH_PRESENT'}
    return {'Accession':acc,'Label':label,'Analysis_branch':'PAIRED_SHORTREAD','Quantification_unit':'PAIRS','Acquisition_source':acquired.source,'Acquisition_format':acquired.fmt,'Acquisition_bytes':acquired.bytes_downloaded,'Read_filter_policy':acquired.read_filter_policy,'SRA_reject_units_removed':_sra_reject_value(acquired),'SRA_reject_pairs_removed':_sra_reject_value(acquired),'SRA_reject_reads_removed':(2*acquired.reject_pairs_removed if _sra_reject_value(acquired)!='NA' else 'NA'),'Metadata_source':md['Metadata_source'],'BioProject':md['BioProject'],'BioSample':md['BioSample'],'Organism':md['Organism'],'Platform':md['Platform'],'Layout':md['Layout'],'Assembly':assembly or 'NA','Assembly_source':assembly_source,'Assembly_lookup_status':lookup_status,'Total_reads':stats['pairs']*2,'Total_pairs':stats['pairs'],'Total_bases':stats['bases'],'Mean_read_length':f"{stats['mean_len']:.3f}",'Q20_pct':f"{stats['q20']:.3f}",'Q30_pct':f"{stats['q30']:.3f}",'GC_pct':f"{stats['gc']:.3f}",'MAPQ30_reads':calls['mapq30'],'Cov80_ID95':calls['id95'],'Cov80_ID98':calls['id98'],'Cov80_ID99':calls['id99'],'Strict_pairs_ID98':calls['strict_pairs'],'Strict_reads_ID98':'NA','Unique_strict_fragments':n,'Unique_strict_reads':'NA','Unique_strict_units':n,'Strict_per_million_pairs':f'{rate:.6f}','Strict_per_million_reads':'NA','Strict_per_million_units':f'{rate:.6f}','Kraken_strict_human':kh,'Kraken_human_fraction':f'{kfrac:.6f}','N_chromosomes':calls['nchr'],'Chromosome_distribution':calls['distribution'],'Top_chromosome':calls['top_chr'],'Top_count':calls['top_n'],'Top_fraction':f"{calls['top_frac']:.6f}",'Bact_ge1_ID95':bact['ge1_95'] if bc_ok else ('NA' if n else 0),'Bact_both_ID95':bact['both_95'] if bc_ok else ('NA' if n else 0),'Bact_ge1_ID98':bact['ge1_98'] if bc_ok else ('NA' if n else 0),'Bact_both_ID98':bact['both_98'] if bc_ok else ('NA' if n else 0),'Human_positive':'YES' if n else 'NO','Assembly_check':acheck,'Runtime_seconds':f'{time.time()-t0:.1f}','Completed_at':now_iso()}


def analyze_ion(s,backend,acc,label,work,acquired,t0,assembly_future=None):
    md=norm(acquired.metadata)
    platform=md['Platform'].replace('-','_')
    if platform not in {x.replace('-','_') for x in ION_PLATFORMS}: raise RuntimeError(f"Unsupported Ion platform {md['Platform']}")
    if md['Layout']!='SINGLE' or acquired.r2 is not None: raise RuntimeError(f"Ion branch expects SINGLE layout with one read file; got {md['Layout']}")

    evidence=s.evidence_root/acc; shutil.rmtree(evidence,ignore_errors=True); evidence.mkdir(parents=True,exist_ok=True)
    stats=seqkit_single_stats(acquired.r1)
    bam=work/f'{acc}.human.bam'; human_map_single(s.human_index,acquired.r1,bam,evidence/'human_minimap2.log',s.threads_per_worker,s.sort_threads)
    calls=strict_human_single_from_bam(bam,evidence); n=int(calls['unique'])
    kh=0; assembly=''; assembly_source='NA'; lookup_status='NOT_NEEDED_NEGATIVE'; acheck='NOT_NEEDED_NEGATIVE'; bact={'ge1_95':0,'both_95':'NA','ge1_98':0,'both_98':'NA'}
    if n:
        sr=extract_strict_single(calls['qnames'],acquired.r1,evidence)
        kh,_=kraken_strict_single(s.kraken_db,sr,evidence,s.kraken_threads)
        res,bf,aerr=_assembly_for_positive(s,acc,md['BioSample'],assembly_future)
        if res is not None:
            assembly=res.accession
            assembly_source=res.source
            lookup_status=res.lookup_status
        else:
            lookup_status='ASSEMBLY_LOOKUP_FAILED'
        if aerr is not None:
            acheck='FAILED'
            _record_assembly_failure(evidence,acc,aerr)
        elif assembly:
            try:
                bact=bacterial_challenge_single(
                    bf,sr,evidence,
                    s.threads_per_worker,s.sort_threads
                )
                acheck='PASS' if bact['ge1_95']==0 else 'BACTERIAL_MATCH_PRESENT'
            except Exception as exc:
                acheck='FAILED'
                _record_assembly_failure(evidence,acc,exc)
        else:
            acheck='NO_EXACT_ASSEMBLY'
        gzip_candidate(sr)

    rej=_sra_reject_value(acquired)
    with (evidence/'run_metadata.tsv').open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t')
        w.writerow(['Accession','Label','Analysis_branch','Quantification_unit','Acquisition_source','Acquisition_format','Read_filter_policy','SRA_reject_units_removed','SRA_reject_pairs_removed','SRA_reject_reads_removed','BioProject','BioSample','Organism','Platform','Layout','Assembly','Assembly_source','Assembly_lookup_status'])
        w.writerow([acc,label,'ION_SINGLE','READS',acquired.source,acquired.fmt,acquired.read_filter_policy,rej,'NA',rej,md['BioProject'],md['BioSample'],md['Organism'],md['Platform'],md['Layout'],assembly,assembly_source,lookup_status])
    rate=n/stats['reads']*1_000_000 if stats['reads'] else 0; kfrac=kh/n if n else 0; bc_ok=acheck in {'PASS','BACTERIAL_MATCH_PRESENT'}
    return {'Accession':acc,'Label':label,'Analysis_branch':'ION_SINGLE','Quantification_unit':'READS','Acquisition_source':acquired.source,'Acquisition_format':acquired.fmt,'Acquisition_bytes':acquired.bytes_downloaded,'Read_filter_policy':acquired.read_filter_policy,'SRA_reject_units_removed':rej,'SRA_reject_pairs_removed':'NA','SRA_reject_reads_removed':rej,'Metadata_source':md['Metadata_source'],'BioProject':md['BioProject'],'BioSample':md['BioSample'],'Organism':md['Organism'],'Platform':md['Platform'],'Layout':md['Layout'],'Assembly':assembly or 'NA','Assembly_source':assembly_source,'Assembly_lookup_status':lookup_status,'Total_reads':stats['reads'],'Total_pairs':'NA','Total_bases':stats['bases'],'Mean_read_length':f"{stats['mean_len']:.3f}",'Q20_pct':f"{stats['q20']:.3f}",'Q30_pct':f"{stats['q30']:.3f}",'GC_pct':f"{stats['gc']:.3f}",'MAPQ30_reads':calls['mapq30'],'Cov80_ID95':calls['id95'],'Cov80_ID98':calls['id98'],'Cov80_ID99':calls['id99'],'Strict_pairs_ID98':'NA','Strict_reads_ID98':calls['strict_reads'],'Unique_strict_fragments':'NA','Unique_strict_reads':n,'Unique_strict_units':n,'Strict_per_million_pairs':'NA','Strict_per_million_reads':f'{rate:.6f}','Strict_per_million_units':f'{rate:.6f}','Kraken_strict_human':kh,'Kraken_human_fraction':f'{kfrac:.6f}','N_chromosomes':calls['nchr'],'Chromosome_distribution':calls['distribution'],'Top_chromosome':calls['top_chr'],'Top_count':calls['top_n'],'Top_fraction':f"{calls['top_frac']:.6f}",'Bact_ge1_ID95':bact['ge1_95'] if bc_ok else ('NA' if n else 0),'Bact_both_ID95':'NA','Bact_ge1_ID98':bact['ge1_98'] if bc_ok else ('NA' if n else 0),'Bact_both_ID98':'NA','Human_positive':'YES' if n else 'NO','Assembly_check':acheck,'Runtime_seconds':f'{time.time()-t0:.1f}','Completed_at':now_iso()}


def analyze_one(s,backend,acc,label,work,acquired,t0,assembly_future=None):
    md=norm(acquired.metadata); platform=md['Platform'].replace('-','_')
    if md['Layout']=='PAIRED' and md['Platform'] in PAIRED_PLATFORMS:
        return analyze_paired(s,backend,acc,label,work,acquired,t0,assembly_future)
    if md['Layout']=='SINGLE' and platform in {x.replace('-','_') for x in ION_PLATFORMS}:
        return analyze_ion(s,backend,acc,label,work,acquired,t0,assembly_future)
    raise RuntimeError(f"Unsupported platform/layout combination: platform={md['Platform']} layout={md['Layout']}; supported: ILLUMINA/BGISEQ/DNBSEQ+PAIRED, ION_TORRENT+SINGLE")

def cleanup_success(s,backend,acc):
    """Called only after RESULT_COMMITTED. Never called after an analysis failure."""
    work=s.work_root/f'{backend}_{acc}'
    a=load_manifest(work,acc,backend.upper())
    if a and a.source in {'NCBI','AWS'} and a.created_sra and s.delete_new_sra_after_success:
        try:
            a.created_sra.unlink(missing_ok=True); p=a.created_sra.parent
            while p!=s.sra_root and p.exists() and not any(p.iterdir()): p.rmdir(); p=p.parent
        except Exception as exc: eprint(f'WARNING cleanup SRA failed for {acc}: {exc}')
    # Remove large analysis/download work only after the result is committed.
    shutil.rmtree(work)
    # If alternate backends have failed partial work trees, they are safe to
    # remove now because this accession has a committed scientific result.
    for other in {'aws','ncbi','ena'}-{backend}:
        shutil.rmtree(s.work_root/f'{other}_{acc}',ignore_errors=True)

def finalize_committed(s,master_store=None,registry_by_accession=None):
    """Finish cleanup after a controller crash.

    In master mode, a RESULT_COMMITTED runtime row is cleaned only after the
    permanent master TSVs are verified complete. If needed, the compact evidence
    is used to repair the master commit first.
    """
    registry_by_accession=registry_by_accession or {}
    for row in db.result_committed_jobs(s.state_db):
        acc=row['accession']; backend=(row.get('backend') or '').lower()
        if backend not in {'aws','ncbi','ena'}:
            eprint(f'WARNING: cannot cleanup committed {acc}: unknown backend {backend}'); continue
        try:
            if master_store is not None and not master_store.is_complete(acc):
                raw=row.get('result_json') or ''
                if not raw: raise RuntimeError(f'{acc}: RESULT_COMMITTED has no result_json for master repair')
                result=json.loads(raw)
                result['Pipeline_version']=master_store.pipeline_version
                result['Scientific_method_version']=master_store.scientific_method_version
                rec=registry_by_accession.get(acc,{'cohort':'','notes':''})
                master_store.commit(result,rec,s.evidence_root/acc)
                eprint(f'[resume-master] COMMITTED {acc}')
            cleanup_success(s,backend,acc); db.mark_done(s.state_db,acc); eprint(f'[resume-cleanup] DONE {acc}')
        except FileNotFoundError:
            if master_store is not None and not master_store.is_complete(acc):
                eprint(f'WARNING: cannot mark {acc} DONE: work absent and master incomplete')
                continue
            db.mark_done(s.state_db,acc); eprint(f'[resume-cleanup] DONE {acc} (work already absent)')
        except Exception as exc:
            eprint(f'WARNING: cleanup/master finalization remains pending for committed {acc}: {exc}')


def worker_loop(s,backend,name,available_backends,master_store=None,registry_by_accession=None,assembly_executor=None,aws_prefetch_manager=None):
    registry_by_accession=registry_by_accession or {}
    eprint(f'[{name}] started backend={backend}')
    # Master mode supplies one shared assembly executor, so assembly resolution
    # is genuinely independent of read workers. Legacy batch mode falls back to
    # one local background thread for compatibility.
    if assembly_executor is None:
        assembly_executor=ThreadPoolExecutor(max_workers=1,thread_name_prefix=f'{name}-assembly')
    while True:
        job=db.claim_job(s.state_db,backend,name)
        if job is None:
            # Important for cross-source failover: stay alive while another worker may
            # return a backend-eligible job to the queue.
            if db.has_unfinished(s.state_db): time.sleep(2); continue
            eprint(f'[{name}] queue empty; exiting'); return
        acc=job['accession']; label=job['label']; work=s.work_root/f'{backend}_{acc}'; work.mkdir(parents=True,exist_ok=True); t0=time.time()
        eprint(f'[{name}] START {acc}')
        try:
            try:
                holder={'future':None}
                if backend=='aws':
                    if aws_prefetch_manager is None:
                        raise RuntimeError('AWS backend requires AwsPrefetchManager')
                    ready=aws_prefetch_manager.get(acc)
                    holder['future']=ready.assembly_future
                    acquired=acquire_aws_prefetched(acc,work,ready.prefetched,s.threads_per_worker)
                else:
                    def metadata_hook(md):
                        if holder['future'] is None and norm(md).get('BioSample'):
                            holder['future']=assembly_executor.submit(_prefetch_exact_assembly,s,acc,dict(md))
                    acquired=_call_acquire_one(s,backend,acc,work,metadata_hook=metadata_hook)
            except UnsupportedRunError as exc:
                db.mark_unsupported(s.state_db,acc,str(exc))
                eprint(f'[{name}] SKIPPED UNSUPPORTED {acc}: {exc}')
                continue
            except Exception as exc:
                status=db.mark_acquisition_exhausted(s.state_db,acc,backend,str(exc),available_backends)
                eprint(f'[{name}] ACQUISITION EXHAUSTED {acc}: {exc}; state={status}')
                continue
            db.mark_stage(s.state_db,acc,'RUNNING_ANALYSIS')
            try:
                eprint(
                    f'[{name}] WAIT ANALYSIS SLOT {acc} '
                    f'(global_slots={_GLOBAL_ANALYSIS_SLOTS})'
                )
                with _GLOBAL_ANALYSIS_SEMAPHORE:
                    eprint(f'[{name}] ENTER ANALYSIS SLOT {acc}')
                    result=_call_analyze_one(
                        s,backend,acc,label,work,acquired,t0,
                        holder.get('future')
                    )
                    eprint(f'[{name}] LEAVE ANALYSIS SLOT {acc}')
                if master_store is not None:
                    result['Pipeline_version']=master_store.pipeline_version
                    result['Scientific_method_version']=master_store.scientific_method_version
                    rec=registry_by_accession.get(acc,{'cohort':'','notes':''})
                    info=master_store.commit(result,rec,s.evidence_root/acc)
                    eprint(f"[{name}] MASTER COMMIT {acc}: result_added={int(info['result_added'])} chromosome_rows_added={info['chromosome_rows_added']}")
            except Exception as exc:
                db.mark_analysis_error(s.state_db,acc,backend,str(exc))
                eprint(f'[{name}] ANALYSIS/MASTER ERROR {acc}: {exc}; downloaded reads retained in {work}')
                continue
            # Permanent master commit precedes cleanup in master mode.
            db.mark_result_committed(s.state_db,acc,result)
            if master_store is None: db.export_master(s.state_db,s.master_summary)
            try:
                cleanup_success(s,backend,acc); db.mark_done(s.state_db,acc)
            except FileNotFoundError:
                db.mark_done(s.state_db,acc)
            except Exception as exc:
                eprint(f'[{name}] RESULT COMMITTED {acc}, but cleanup pending: {exc}')
                continue
            if master_store is None: db.export_master(s.state_db,s.master_summary)
            eprint(f"[{name}] END {acc} branch={result.get('Analysis_branch','NA')} unique={result.get('Unique_strict_units','NA')} rate={result.get('Strict_per_million_units','NA')}/M-{str(result.get('Quantification_unit','units')).lower()}")
        except Exception as exc:
            # Defensive catch: retain local work and make the failure retryable on same backend.
            db.mark_analysis_error(s.state_db,acc,backend,str(exc)); eprint(f'[{name}] UNEXPECTED ERROR {acc}: {exc}; work retained')
