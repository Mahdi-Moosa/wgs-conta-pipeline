from __future__ import annotations
import csv,json,re,shutil,subprocess,time,urllib.parse
from dataclasses import asdict,dataclass
from pathlib import Path
from typing import Dict,List,Optional
from .util import capture,eprint,md5_file,retry_run,retry_shell_capture,run

@dataclass
class AcquiredReads:
    accession:str; r1:Path; r2:Optional[Path]; metadata:Dict[str,str]; source:str; fmt:str; bytes_downloaded:int; cleanup_paths:List[Path]; created_sra:Optional[Path]=None; reject_pairs_removed:int=0; read_filter_policy:str='NONE'



class UnsupportedRunError(RuntimeError):
    pass

@dataclass
class PrefetchedSRA:
    accession:str
    sra:Path
    metadata:Dict[str,str]
    bytes_downloaded:int
    created:bool


def _norm_platform(x):
    return (x or '').upper().replace('SRA_PLATFORM_','').replace('-','_')

def _supported_combo(md):
    p=_norm_platform(md.get('Platform','')); layout=(md.get('LibraryLayout') or '').upper()
    return (layout=='PAIRED' and p in {'ILLUMINA','BGISEQ','DNBSEQ'}) or (layout=='SINGLE' and p in {'ION_TORRENT','IONTORRENT','ION TORRENT'})

def _known_combo(md):
    return bool(md.get('Platform')) and bool(md.get('LibraryLayout'))

def _manifest_path(work): return Path(work)/'acquisition_manifest.json'
def _save_manifest(work,a:AcquiredReads):
    d={'accession':a.accession,'r1':str(a.r1),'r2':str(a.r2) if a.r2 else '','metadata':a.metadata,'source':a.source,'fmt':a.fmt,'bytes_downloaded':a.bytes_downloaded,'cleanup_paths':[str(x) for x in a.cleanup_paths],'created_sra':str(a.created_sra) if a.created_sra else '','reject_pairs_removed':int(a.reject_pairs_removed),'read_filter_policy':a.read_filter_policy}
    p=_manifest_path(work); p.parent.mkdir(parents=True,exist_ok=True); tmp=p.with_suffix('.tmp'); tmp.write_text(json.dumps(d,indent=2,sort_keys=True)); tmp.replace(p)
def load_manifest(work,accession=None,source=None):
    p=_manifest_path(work)
    if not p.exists(): return None
    try:
        d=json.loads(p.read_text()); r1=Path(d['r1']); r2=Path(d['r2']) if d.get('r2') else None
        if accession and d.get('accession')!=accession: return None
        if source and d.get('source')!=source: return None
        if not (r1.exists() and r1.stat().st_size>0): return None
        if r2 is not None and not (r2.exists() and r2.stat().st_size>0): return None
        cs=Path(d['created_sra']) if d.get('created_sra') else None
        return AcquiredReads(d['accession'],r1,r2,d.get('metadata',{}),d['source'],d.get('fmt',''),int(d.get('bytes_downloaded',0)),[Path(x) for x in d.get('cleanup_paths',[])],cs,int(d.get('reject_pairs_removed',0)),d.get('read_filter_policy','NONE'))
    except Exception as exc:
        eprint(f'WARNING: ignoring invalid acquisition manifest {p}: {exc}'); return None

def parse_runinfo_csv(text,accession):
    rows=list(csv.DictReader(text.splitlines()))
    if not rows: return {}
    row=next((x for x in rows if x.get('Run')==accession),rows[0])
    return {k:(v or '').strip() for k,v in row.items()}

def ncbi_runinfo(accession,retries=4):
    cmd=f"esearch -db sra -query {accession} </dev/null | efetch -format runinfo"
    return parse_runinfo_csv(retry_shell_capture(cmd,tries=retries),accession)

def find_sra(root,accession):
    for p in Path(root).rglob(f'{accession}.sra'): return p
    return None

def _vdb_platform_layout(sra):
    platform='UNKNOWN'; layout='UNKNOWN'
    p=subprocess.run(['vdb-dump',str(sra),'--info'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
    if p.returncode==0:
        for line in p.stdout.splitlines():
            if line.lower().startswith('platf'):
                platform=line.split(':',1)[-1].strip().replace('SRA_PLATFORM_',''); break
    p=subprocess.run(['vdb-dump',str(sra),'-R1','-C','READ_TYPE'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
    txt=p.stdout if p.returncode==0 else ''
    if ',' in txt: layout='PAIRED'
    elif 'BIOLOGICAL' in txt: layout='SINGLE'
    return platform,layout

def ena_metadata(accession,api='https://www.ebi.ac.uk/ena/portal/api/filereport',retries=4):
    fields=['run_accession','study_accession','sample_accession','scientific_name','library_layout','instrument_platform','fastq_ftp','fastq_aspera','fastq_md5','fastq_bytes']
    q=urllib.parse.urlencode({'accession':accession,'result':'read_run','fields':','.join(fields),'format':'tsv'})
    txt=capture(['curl','-fsSL','--retry',str(retries),'--retry-all-errors','--connect-timeout','30',f'{api}?{q}'])
    rows=list(csv.DictReader(txt.splitlines(),delimiter='\t'))
    if not rows: raise RuntimeError(f'ENA returned no read_run record for {accession}')
    r=rows[0]
    return {'Run':r.get('run_accession',''),'BioProject':r.get('study_accession',''),'BioSample':r.get('sample_accession',''),'ScientificName':r.get('scientific_name',''),'LibraryLayout':r.get('library_layout',''),'Platform':r.get('instrument_platform',''),'fastq_ftp':r.get('fastq_ftp',''),'fastq_aspera':r.get('fastq_aspera',''),'fastq_md5':r.get('fastq_md5',''),'fastq_bytes':r.get('fastq_bytes',''),'_metadata_source':'ENA'}

def _sra_reject_spots(sra, reject_path):
    """Return spot IDs where any biological read is explicitly READ_FILTER_REJECT.

    PASS and CRITERIA are retained. For paired-fragment analysis, a spot is removed
    as a unit if either mate is REJECT so the NCBI branch cannot gain orphaned or
    rejected fragments relative to ENA archive-generated paired FASTQs.
    """
    reject_path=Path(reject_path); reject_path.parent.mkdir(parents=True,exist_ok=True)
    p=subprocess.Popen(['vdb-dump',str(sra),'-C','SPOT_ID,READ_TYPE,READ_FILTER','-f','tab'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    rejects=[]; rows=0
    for raw in p.stdout:
        line=raw.rstrip('\n')
        if not line: continue
        f=line.split('\t')
        if len(f)<3: continue
        rows+=1
        read_types=[x.strip() for x in f[1].split(',')]
        read_filters=[x.strip() for x in f[2].split(',')]
        # Consider only biological reads; technical/index reads are not emitted by the
        # paired fasterq workflow and therefore must not cause a biological pair drop.
        biological_reject=False
        for i,flt in enumerate(read_filters):
            typ=read_types[i] if i < len(read_types) else ''
            if 'BIOLOGICAL' in typ and 'SRA_READ_FILTER_REJECT' in flt:
                biological_reject=True; break
        if biological_reject:
            try: rejects.append(int(f[0].strip()))
            except ValueError: raise RuntimeError(f'Could not parse SRA SPOT_ID from READ_FILTER row: {line[:200]}')
    err=p.stderr.read(); rc=p.wait()
    if rc!=0: raise RuntimeError(f'vdb-dump READ_FILTER failed for {sra} (exit={rc}): {err[-2000:]}')
    if rows==0: raise RuntimeError(f'vdb-dump returned no READ_FILTER rows for {sra}')
    with reject_path.open('w') as fh:
        for spot in rejects: fh.write(f'{spot}\n')
    eprint(f'+ SRA READ_FILTER audit: {rows} spots; {len(rejects)} spot(s) with rejected biological mate(s) to exclude; biological PASS+CRITERIA retained')
    return rejects,rows

def _sra_biological_read_indices(sra):
    """Return 1-based SRA read indices marked BIOLOGICAL.

    fasterq-dump --split-files preserves the original SRA read index.
    Therefore a legacy spot such as:

        BIOLOGICAL, TECHNICAL, BIOLOGICAL

    produces biological FASTQs ``_1.fastq`` and ``_3.fastq`` rather
    than ``_1.fastq`` and ``_2.fastq``.

    Likewise:

        TECHNICAL, BIOLOGICAL

    produces ``_2.fastq`` for a biologically single-end run.
    """
    p=subprocess.run(
        [
            'vdb-dump',str(sra),
            '-R','1',
            '-C','READ_TYPE',
            '-f','tab'
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if p.returncode!=0:
        raise RuntimeError(
            f'vdb-dump READ_TYPE failed for {sra} '
            f'(exit={p.returncode}): {p.stderr[-2000:]}'
        )

    line=next(
        (x.strip() for x in p.stdout.splitlines() if x.strip()),
        ''
    )

    if not line:
        raise RuntimeError(
            f'vdb-dump returned no READ_TYPE row for {sra}'
        )

    read_types=[x.strip() for x in line.split(',')]

    indices=[
        i+1
        for i,typ in enumerate(read_types)
        if 'BIOLOGICAL' in typ
    ]

    if not indices:
        raise RuntimeError(
            f'No BIOLOGICAL reads detected in SRA READ_TYPE for '
            f'{sra}: {read_types}'
        )

    return indices


def _select_sra_biological_fastqs(accession,sra,fq,layout):
    """Select fasterq-dump output using actual SRA biological read indices."""
    fq=Path(fq)
    bio_idx=_sra_biological_read_indices(sra)

    if layout=='PAIRED':
        if len(bio_idx)!=2:
            raise RuntimeError(
                f'PAIRED run {accession} has {len(bio_idx)} '
                f'biological read indices {bio_idx}; expected exactly 2'
            )

        r1=fq/f'{accession}_{bio_idx[0]}.fastq'
        r2=fq/f'{accession}_{bio_idx[1]}.fastq'

        missing=[
            x.name
            for x in (r1,r2)
            if not x.exists() or x.stat().st_size==0
        ]

        if missing:
            produced=sorted(
                x.name
                for x in fq.glob(f'{accession}*.fastq')
                if x.is_file()
            )

            raise RuntimeError(
                f'Expected biological paired FASTQs not produced for '
                f'{accession}; biological_indices={bio_idx}; '
                f'missing={missing}; produced={produced}'
            )

        return r1,r2,bio_idx

    if layout=='SINGLE':
        if len(bio_idx)!=1:
            raise RuntimeError(
                f'SINGLE run {accession} has {len(bio_idx)} '
                f'biological read indices {bio_idx}; expected exactly 1'
            )

        # Depending on SRA structure fasterq-dump may emit a truly
        # single biological read as accession.fastq, accession_1.fastq,
        # or (when preceded by a technical read) accession_N.fastq.
        indexed=fq/f'{accession}_{bio_idx[0]}.fastq'
        plain=fq/f'{accession}.fastq'

        candidates=[]
        for x in (plain,indexed):
            if x not in candidates:
                candidates.append(x)

        present=[
            x for x in candidates
            if x.exists() and x.stat().st_size>0
        ]

        if len(present)!=1:
            produced=sorted(
                x.name
                for x in fq.glob(f'{accession}*.fastq')
                if x.is_file()
            )

            raise RuntimeError(
                f'Expected exactly one biological single-end FASTQ for '
                f'{accession}; biological_indices={bio_idx}; '
                f'candidate_matches={[x.name for x in present]}; '
                f'produced={produced}'
            )

        return present[0],None,bio_idx

    raise RuntimeError(
        f'Unsupported/unknown SRA library layout for '
        f'{accession}: {layout}'
    )


def _seqkit_num_seqs(path):
    """Return the number of FASTQ records reported by seqkit stats."""
    txt=capture(['seqkit','stats','-T',str(path)])
    lines=[x for x in txt.splitlines() if x.strip()]
    if len(lines)<2: raise RuntimeError(f'Unexpected seqkit stats output for {path}')
    h=lines[0].split('\t'); row=dict(zip(h,lines[1].split('\t')))
    try: return int(float(str(row['num_seqs']).replace(',','')))
    except Exception as exc: raise RuntimeError(f'Could not parse seqkit num_seqs for {path}: {row}') from exc


def _filter_reject_spots(accession,r1,r2,reject_spots,work):
    """Remove whole rejected spots from fasterq-dump split FASTQs.

    fasterq-dump --split-files uses the same spot identifier for both mates
    (e.g. ``SYNTHETIC_RUN.12345``); mate number is in the FASTQ description,
    not a ``/1`` or ``/2`` suffix on the sequence ID.  Therefore one exact
    sequence-ID pattern per rejected spot removes that spot from both mate files.
    """
    r1=Path(r1); r2=Path(r2)
    raw1=_seqkit_num_seqs(r1); raw2=_seqkit_num_seqs(r2)
    if raw1!=raw2: raise RuntimeError(f'Unequal pre-normalization R1/R2 counts for {accession}: {raw1}, {raw2}')
    if not reject_spots:
        eprint(f'+ SRA reject normalization complete: no rejected paired spots in {accession}; {raw1} pairs retained')
        return raw1
    rf=Path(work)/'read_filter'; rf.mkdir(parents=True,exist_ok=True)
    patterns=rf/'reject_qnames.txt'
    with patterns.open('w') as fh:
        for spot in reject_spots:
            fh.write(f'{accession}.{spot}\n')
    for src in (r1,r2):
        tmp=src.with_suffix(src.suffix+'.no_reject.tmp')
        run(['seqkit','grep','-v','-f',str(patterns),str(src),'-o',str(tmp)])
        if not tmp.exists() or tmp.stat().st_size==0: raise RuntimeError(f'SRA reject normalization produced empty/missing FASTQ: {tmp}')
        tmp.replace(src)
    norm1=_seqkit_num_seqs(r1); norm2=_seqkit_num_seqs(r2)
    if norm1!=norm2: raise RuntimeError(f'Unequal post-normalization R1/R2 counts for {accession}: {norm1}, {norm2}')
    removed=raw1-norm1
    if removed!=len(reject_spots):
        raise RuntimeError(
            f'SRA reject normalization count mismatch for {accession}: '
            f'expected to remove {len(reject_spots)} paired spots, actually removed {removed} '
            f'(raw={raw1}, normalized={norm1}). Refusing to analyze an unverified source representation.'
        )
    eprint(f'+ SRA reject normalization complete: removed {removed} paired spot(s) from {accession}; {norm1} pairs retained')
    return norm1


def _filter_reject_single(accession,reads,reject_spots,work):
    """Remove SRA READ_FILTER_REJECT biological spots from a single-end FASTQ."""
    reads=Path(reads); raw=_seqkit_num_seqs(reads)
    if not reject_spots:
        eprint(f'+ SRA reject normalization complete: no rejected single-end spots in {accession}; {raw} reads retained')
        return raw
    rf=Path(work)/'read_filter'; rf.mkdir(parents=True,exist_ok=True)
    patterns=rf/'reject_qnames.txt'
    with patterns.open('w') as fh:
        for spot in reject_spots: fh.write(f'{accession}.{spot}\n')
    tmp=reads.with_suffix(reads.suffix+'.no_reject.tmp')
    run(['seqkit','grep','-v','-f',str(patterns),str(reads),'-o',str(tmp)])
    if not tmp.exists() or tmp.stat().st_size==0: raise RuntimeError(f'SRA reject normalization produced empty/missing FASTQ: {tmp}')
    tmp.replace(reads)
    norm=_seqkit_num_seqs(reads); removed=raw-norm
    if removed!=len(reject_spots):
        raise RuntimeError(
            f'SRA single-end reject normalization count mismatch for {accession}: '
            f'expected to remove {len(reject_spots)} spots, actually removed {removed} '
            f'(raw={raw}, normalized={norm}). Refusing to analyze an unverified source representation.'
        )
    eprint(f'+ SRA reject normalization complete: removed {removed} single-end spot(s) from {accession}; {norm} reads retained')
    return norm

def acquire_ncbi(accession,work,sra_root,threads,retries,metadata_hook=None):
    cached=load_manifest(work,accession,'NCBI')
    required_policy='EXCLUDE_REJECT_KEEP_PASS_CRITERIA'
    valid_norm_versions={'2.2.3','2.3.0','2.4.0','2.4.1','2.4.2','2.4.3','2.4.4'}
    if cached and cached.read_filter_policy==required_policy and cached.metadata.get('_sra_normalization_version') in valid_norm_versions:
        desc=cached.r1.name if cached.r2 is None else f'{cached.r1.name}, {cached.r2.name}'
        eprint(f'+ [resume] reusing normalized NCBI reads for {accession}: {desc}')
        if metadata_hook:
            try: metadata_hook(cached.metadata)
            except Exception as exc: eprint(f'WARNING: metadata hook failed for cached {accession}: {exc}')
        return cached
    if cached:
        eprint(f'+ [upgrade] discarding incompatible NCBI normalization work for {accession}')
        shutil.rmtree(Path(work)/'reads',ignore_errors=True); shutil.rmtree(Path(work)/'fasterq_tmp',ignore_errors=True); shutil.rmtree(Path(work)/'read_filter',ignore_errors=True)
        _manifest_path(work).unlink(missing_ok=True)

    # Cheap metadata preflight before a potentially large NCBI download. If the
    # platform/layout is authoritatively known and unsupported, skip acquisition.
    md={}; msrc='LOCAL_ONLY'
    # AWS bulk acquisition can run four transfers concurrently. Prefer ENA's
    # lightweight run metadata endpoint here so those four preflights do not
    # multiply NCBI E-utility traffic from one public IP.
    try:
        md=ena_metadata(accession,retries=retries); msrc='ENA'
    except Exception as exc:
        eprint(f'WARNING: ENA metadata preflight failed for {accession}: {exc}')
        try: md=ncbi_runinfo(accession,retries); msrc='NCBI_FALLBACK'
        except Exception: md={}; msrc='LOCAL_ONLY'
    if _known_combo(md) and not _supported_combo(md):
        raise UnsupportedRunError(f"Unsupported platform/layout before NCBI acquisition for {accession}: platform={md.get('Platform','')} layout={md.get('LibraryLayout','')}")
    if metadata_hook:
        try: metadata_hook(md)
        except Exception as exc: eprint(f'WARNING: metadata hook failed for {accession}: {exc}')

    sra_root=Path(sra_root); sra_root.mkdir(parents=True,exist_ok=True)
    sra=find_sra(sra_root,accession); created=None
    if sra is None:
        retry_run(['prefetch',accession,'-O',str(sra_root)],tries=retries)
        sra=find_sra(sra_root,accession)
        if sra is None: raise RuntimeError(f'prefetch completed but {accession}.sra was not found')
        created=sra
    vp,vl=_vdb_platform_layout(sra)
    if not md.get('Platform'): md['Platform']=vp
    if not md.get('LibraryLayout'): md['LibraryLayout']=vl
    md['_metadata_source']=msrc
    layout=(md.get('LibraryLayout') or vl or '').upper()
    if _known_combo(md) and not _supported_combo(md):
        raise UnsupportedRunError(f"Unsupported platform/layout from NCBI/local metadata for {accession}: platform={md.get('Platform','')} layout={layout}")

    rf=Path(work)/'read_filter'; rejects,total_spots=_sra_reject_spots(sra,rf/'reject_spot_ids.txt')
    fq=Path(work)/'reads'; tmp=Path(work)/'fasterq_tmp'; fq.mkdir(parents=True,exist_ok=True); tmp.mkdir(parents=True,exist_ok=True)
    run(['fasterq-dump',str(sra),'--split-files','--threads',str(threads),'--temp',str(tmp),'--outdir',str(fq)])

    r1,r2,bio_idx=_select_sra_biological_fastqs(
        accession,sra,fq,layout
    )

    eprint(
        f'+ SRA biological read indices for {accession}: '
        f'{",".join(map(str,bio_idx))}'
    )

    if layout=='PAIRED':
        normalized_units=_filter_reject_spots(
            accession,r1,r2,rejects,work
        )
        fmt='SRA_REJECT_NORMALIZED'; unit='PAIRS'
    elif layout=='SINGLE':
        normalized_units=_filter_reject_single(
            accession,r1,rejects,work
        )
        fmt='SRA_REJECT_NORMALIZED_SINGLE'; unit='READS'
    else:
        raise RuntimeError(
            f'Unsupported/unknown NCBI library layout for '
            f'{accession}: {layout}'
        )

    md['_sra_read_filter_policy']='EXCLUDE_REJECT_KEEP_PASS_CRITERIA'
    md['_sra_biological_read_indices']=','.join(map(str,bio_idx))
    md['_sra_reject_units_removed']=str(len(rejects))
    md['_sra_total_spots_before_filter']=str(total_spots)
    md['_sra_normalized_units']=str(normalized_units)
    md['_sra_normalization_unit']=unit
    md['_sra_normalization_version']='2.4.4'
    a=AcquiredReads(accession,r1,r2,md,'NCBI',fmt,sra.stat().st_size,[fq,tmp,rf],created,len(rejects),'EXCLUDE_REJECT_KEEP_PASS_CRITERIA')
    _save_manifest(work,a); return a

def _validate_sra(path):
    """Run the SRA toolkit's internal integrity validation."""
    p=subprocess.run(['vdb-validate',str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    if p.returncode!=0:
        msg=(p.stderr or p.stdout or '')[-4000:]
        raise RuntimeError(f'vdb-validate failed for {path}: {msg}')
    return True


def _aws_metadata_preflight(accession,retries,metadata_hook=None):
    """Resolve run metadata before AWS transfer so assembly prefetch can overlap it."""
    md={}; msrc='LOCAL_ONLY'
    try:
        md=ncbi_runinfo(accession,retries); msrc='NCBI'
    except Exception as exc:
        eprint(f'WARNING: NCBI RunInfo preflight failed for {accession}: {exc}')
        try: md=ena_metadata(accession,retries=retries); msrc='ENA_FALLBACK'
        except Exception: md={}; msrc='LOCAL_ONLY'
    if _known_combo(md) and not _supported_combo(md):
        raise UnsupportedRunError(
            f"Unsupported platform/layout before AWS acquisition for {accession}: "
            f"platform={md.get('Platform','')} layout={md.get('LibraryLayout','')}"
        )
    md['_metadata_source']=msrc
    if metadata_hook:
        try: metadata_hook(md)
        except Exception as exc: eprint(f'WARNING: metadata hook failed for {accession}: {exc}')
    return md


def prefetch_aws_sra(accession,sra_root,aws_cli,aws_bucket,retries,metadata_hook=None):
    """Download one normalized public SRA object from the AWS Open Data mirror.

    The object is staged atomically under ``sra_root/aws_cache`` and validated
    with ``vdb-validate`` before becoming visible to downstream normalization.
    Known-bad partial/invalid objects are removed before a fresh retry.
    """
    md=_aws_metadata_preflight(accession,retries,metadata_hook=metadata_hook)
    root=Path(sra_root)/'aws_cache'/accession
    root.mkdir(parents=True,exist_ok=True)
    final=root/f'{accession}.sra'
    if final.exists() and final.stat().st_size>0:
        try:
            _validate_sra(final)
            eprint(f'+ [AWS cache] VALID {accession}: {final}')
            return PrefetchedSRA(accession,final,md,final.stat().st_size,False)
        except Exception as exc:
            eprint(f'+ [AWS cache] INVALID {accession}: {exc}; removing before fresh download')
            final.unlink(missing_ok=True)

    cli=str(Path(aws_cli).expanduser())
    obj=f"{str(aws_bucket).rstrip('/')}/{accession}/{accession}"
    last=''
    for attempt in range(1,retries+1):
        tmp=final.with_name(f'{accession}.download.sra')
        tmp.unlink(missing_ok=True)
        eprint(f'+ [AWS SRA attempt {attempt}/{retries}] {obj}')
        p=subprocess.run([
            cli,'s3','cp',obj,str(tmp),'--no-sign-request','--only-show-errors'
        ],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if p.returncode==0 and tmp.exists() and tmp.stat().st_size>0:
            try:
                _validate_sra(tmp)
                tmp.replace(final)
                eprint(f'+ AWS SRA VALID {accession}: {final.stat().st_size/1024/1024:.1f} MiB')
                return PrefetchedSRA(accession,final,md,final.stat().st_size,True)
            except Exception as exc:
                last=str(exc)
                eprint(f'+ AWS SRA VALIDATION FAIL {accession}: {exc}; removing corrupt object')
                tmp.unlink(missing_ok=True)
        else:
            last=(p.stderr or p.stdout or f'exit={p.returncode}')[-4000:]
            tmp.unlink(missing_ok=True)
        if attempt<retries: time.sleep(min(2*attempt,8))
    raise RuntimeError(f'AWS SRA acquisition failed after {retries} attempts for {accession}: {last}')


def acquire_aws_prefetched(accession,work,prefetched,threads):
    """Normalize a validated AWS SRA object using the frozen v2.2.3 READ_FILTER policy."""
    cached=load_manifest(work,accession,'AWS')
    required_policy='EXCLUDE_REJECT_KEEP_PASS_CRITERIA'
    valid_norm_versions={'2.2.3','2.3.0','2.4.0','2.4.1','2.4.2','2.4.3','2.4.4'}
    if cached and cached.read_filter_policy==required_policy and cached.metadata.get('_sra_normalization_version') in valid_norm_versions:
        desc=cached.r1.name if cached.r2 is None else f'{cached.r1.name}, {cached.r2.name}'
        eprint(f'+ [resume] reusing normalized AWS reads for {accession}: {desc}')
        return cached
    if cached:
        eprint(f'+ [upgrade] discarding incompatible AWS normalization work for {accession}')
        shutil.rmtree(Path(work)/'reads',ignore_errors=True)
        shutil.rmtree(Path(work)/'fasterq_tmp',ignore_errors=True)
        shutil.rmtree(Path(work)/'read_filter',ignore_errors=True)
        _manifest_path(work).unlink(missing_ok=True)

    sra=Path(prefetched.sra); md=dict(prefetched.metadata)
    if not sra.exists() or sra.stat().st_size==0:
        raise RuntimeError(f'Prefetched AWS SRA is missing/empty for {accession}: {sra}')
    vp,vl=_vdb_platform_layout(sra)
    if not md.get('Platform'): md['Platform']=vp
    if not md.get('LibraryLayout'): md['LibraryLayout']=vl
    md['_metadata_source']=md.get('_metadata_source') or 'NCBI'
    layout=(md.get('LibraryLayout') or vl or '').upper()
    if _known_combo(md) and not _supported_combo(md):
        raise UnsupportedRunError(
            f"Unsupported platform/layout from AWS SRA metadata for {accession}: "
            f"platform={md.get('Platform','')} layout={layout}"
        )

    rf=Path(work)/'read_filter'; rejects,total_spots=_sra_reject_spots(sra,rf/'reject_spot_ids.txt')
    fq=Path(work)/'reads'; tmp=Path(work)/'fasterq_tmp'
    fq.mkdir(parents=True,exist_ok=True); tmp.mkdir(parents=True,exist_ok=True)
    run(['fasterq-dump',str(sra),'--split-files','--threads',str(threads),'--temp',str(tmp),'--outdir',str(fq)])

    r1,r2,bio_idx=_select_sra_biological_fastqs(
        accession,sra,fq,layout
    )

    eprint(
        f'+ SRA biological read indices for {accession}: '
        f'{",".join(map(str,bio_idx))}'
    )

    if layout=='PAIRED':
        normalized_units=_filter_reject_spots(
            accession,r1,r2,rejects,work
        )
        fmt='AWS_SRA_REJECT_NORMALIZED'; unit='PAIRS'
    elif layout=='SINGLE':
        normalized_units=_filter_reject_single(
            accession,r1,rejects,work
        )
        fmt='AWS_SRA_REJECT_NORMALIZED_SINGLE'; unit='READS'
    else:
        raise RuntimeError(
            f'Unsupported/unknown AWS SRA library layout for '
            f'{accession}: {layout}'
        )

    md['_sra_read_filter_policy']=required_policy
    md['_sra_biological_read_indices']=','.join(map(str,bio_idx))
    md['_sra_reject_units_removed']=str(len(rejects))
    md['_sra_total_spots_before_filter']=str(total_spots)
    md['_sra_normalized_units']=str(normalized_units)
    md['_sra_normalization_unit']=unit
    md['_sra_normalization_version']='2.4.4'
    a=AcquiredReads(
        accession,r1,r2,md,'AWS',fmt,prefetched.bytes_downloaded,[fq,tmp,rf],
        (sra if prefetched.created else None),len(rejects),required_policy
    )
    _save_manifest(work,a)
    return a


def _ena_https(x):
    x=x.strip()
    if x.startswith('ftp://'): return 'https://'+x[6:]
    if x.startswith(('http://','https://')): return x
    return 'https://'+x

def _ena_aspera_source(x):
    x=x.strip()
    if not x: return x
    if '@' in x.split(':',1)[0]: return x
    return 'era-fasp@'+x

def _pair_rank(name):
    if re.search(r'_1(?:\.f(?:ast)?q)?\.gz$',name): return 1
    if re.search(r'_2(?:\.f(?:ast)?q)?\.gz$',name): return 2
    return 9

def _download_https(url,dest,retries):
    t0=time.time(); retry_run(['curl','-fL','--retry',str(retries),'--retry-all-errors','--connect-timeout','30','-o',str(dest),url],tries=1)
    sec=max(time.time()-t0,0.001); mb=dest.stat().st_size/1024/1024 if dest.exists() else 0
    eprint(f'+ HTTPS complete {dest.name}: {mb:.1f} MiB in {sec:.1f}s ({mb*8/sec:.1f} Mib/s)')

def _download_aspera(source,dest_dir,ascp,key,port,rate,retries,expected_name=''):
    ascp=Path(ascp).expanduser(); key=Path(key).expanduser(); dest_dir=Path(dest_dir)
    if not ascp.exists(): raise RuntimeError(f'Aspera ascp not found: {ascp}')
    if not key.exists(): raise RuntimeError(f'Aspera bypass key not found: {key}')
    last=None
    for attempt in range(1,retries+1):
        eprint(f'+ [Aspera attempt {attempt}/{retries}] {source}')
        t0=time.time()
        p=subprocess.run([str(ascp),'-T','-P',str(port),'-l',str(rate),'-i',str(key),source,str(dest_dir) + '/'])
        if p.returncode==0:
            name=expected_name or Path(source.split(':',1)[-1]).name; dest=dest_dir/name; sec=max(time.time()-t0,0.001)
            if dest.exists():
                mb=dest.stat().st_size/1024/1024; eprint(f'+ Aspera complete {dest.name}: {mb:.1f} MiB in {sec:.1f}s ({mb*8/sec:.1f} Mib/s)')
            return
        last=p.returncode
        if attempt<retries: time.sleep(min(5*attempt,20))
    raise RuntimeError(f'Aspera transfer failed after {retries} attempts (exit={last}): {source}')

def _download_ena_file_verified(accession,name,https_url,aspera_source,dest,fq,retries,verify_md5,expected_md5,chosen,ascp,ascp_key,ascp_port,ascp_rate,https_fallback):
    """Download one ENA FASTQ and verify integrity before it can be used.

    Transport failures retain the existing Aspera retry/resume behaviour.  A
    *successful* transfer whose MD5 is wrong is different: the local file is
    known to be corrupt, so it is unlinked immediately and the same mate is
    acquired again from scratch.  A valid mate already downloaded earlier in
    the pair is left untouched.
    """
    md5_attempts=retries if verify_md5 and expected_md5 else 1
    # A previous interrupted/failed run may have left an unmanifested FASTQ.
    # If ENA supplied an expected MD5, validate it before any new transfer: a
    # complete valid file can be reused, while a stale/corrupt one is removed
    # before the first acquisition attempt.
    if dest.exists() and verify_md5 and expected_md5:
        got_existing=md5_file(dest)
        if got_existing.lower()==expected_md5.lower():
            eprint(f'+ [resume] existing ENA FASTQ MD5 PASS {dest.name}; reusing file')
            return dest.stat().st_size
        eprint(f'+ [resume] existing ENA FASTQ MD5 FAIL {dest.name}: got {got_existing}, expected {expected_md5}; removing before fresh acquisition')
        dest.unlink(missing_ok=True)
    for integrity_attempt in range(1,md5_attempts+1):
        if integrity_attempt>1:
            # A checksum-failed file must never be resumed/reused.
            dest.unlink(missing_ok=True)
        if chosen=='aspera':
            try:
                _download_aspera(aspera_source,fq,ascp,ascp_key,ascp_port,ascp_rate,retries,name)
            except Exception:
                if not https_fallback: raise
                eprint(f'WARNING: Aspera failed for {accession}/{name}; falling back to HTTPS')
                dest.unlink(missing_ok=True); _download_https(https_url,dest,retries)
        else:
            _download_https(https_url,dest,retries)
        if not dest.exists(): raise RuntimeError(f'ENA download completed but file is missing: {dest}')
        if verify_md5 and expected_md5:
            got=md5_file(dest)
            if got.lower()!=expected_md5.lower():
                eprint(f'+ MD5 FAIL {dest.name}: got {got}, expected {expected_md5}')
                # Delete immediately: this path is now known-bad and must not
                # survive a failed accession or contaminate a later resume.
                dest.unlink(missing_ok=True)
                if integrity_attempt<md5_attempts:
                    eprint(f'+ removed corrupt {dest.name}; fresh acquisition attempt {integrity_attempt+1}/{md5_attempts}')
                    time.sleep(min(2*integrity_attempt,8))
                    continue
                raise RuntimeError(
                    f'ENA MD5 mismatch for {dest.name} after {md5_attempts} fresh acquisition attempt(s): '
                    f'got {got}, expected {expected_md5}; corrupt file removed'
                )
            eprint(f'+ MD5 PASS {dest.name}')
        return dest.stat().st_size
    raise RuntimeError(f'Unreachable ENA integrity retry state for {accession}/{name}')


def acquire_ena(accession,work,api,retries,verify_md5=True,transport='aspera',ascp='~/.aspera/sdk/ascp',ascp_key='~/.aspera/sdk/aspera_bypass_rsa.pem',ascp_port=33001,ascp_rate='300m',https_fallback=False,metadata_hook=None):
    cached=load_manifest(work,accession,'ENA')
    if cached:
        if cached.read_filter_policy=='NONE':
            cached.read_filter_policy='ENA_ARCHIVE_GENERATED_FASTQ'; _save_manifest(work,cached)
        desc=cached.r1.name if cached.r2 is None else f'{cached.r1.name}, {cached.r2.name}'
        eprint(f'+ [resume] reusing ENA reads for {accession}: {desc}')
        if metadata_hook:
            try: metadata_hook(cached.metadata)
            except Exception as exc: eprint(f'WARNING: metadata hook failed for cached {accession}: {exc}')
        return cached
    md=ena_metadata(accession,api,retries)
    if metadata_hook:
        try: metadata_hook(md)
        except Exception as exc: eprint(f'WARNING: metadata hook failed for {accession}: {exc}')
    layout=(md.get('LibraryLayout') or '').upper()
    if not _supported_combo(md):
        raise UnsupportedRunError(f"Unsupported platform/layout from ENA metadata for {accession}: platform={md.get('Platform','')} layout={layout}")
    ftps=[x for x in md.get('fastq_ftp','').split(';') if x]; asperas=[x for x in md.get('fastq_aspera','').split(';') if x]; md5s=[x for x in md.get('fastq_md5','').split(';') if x]
    expected=2 if layout=='PAIRED' else 1
    if len(ftps)!=expected: raise RuntimeError(f'ENA expected {expected} generated FASTQ file(s) for {accession}; got {len(ftps)}')
    if asperas and len(asperas)!=expected: raise RuntimeError(f'ENA returned unexpected Aspera file count for {accession}: {len(asperas)}')
    fq=Path(work)/'reads'; fq.mkdir(parents=True,exist_ok=True); outputs=[]; downloaded=0
    mode=(transport or 'aspera').lower()
    if mode not in {'aspera','https','auto'}: raise RuntimeError(f'Unsupported ENA transport: {mode}')
    aspera_available=Path(ascp).expanduser().exists() and Path(ascp_key).expanduser().exists() and len(asperas)==expected
    if mode=='aspera' and not aspera_available: raise RuntimeError('ENA Aspera requested but ascp/key/fastq_aspera is unavailable')
    chosen='aspera' if (mode=='aspera' or (mode=='auto' and aspera_available)) else 'https'
    # Sequential files are intentional: concurrent FASP sessions were empirically unstable on this route.
    for i,u0 in enumerate(ftps):
        https_url=_ena_https(u0); name=Path(urllib.parse.urlparse(https_url).path).name or (f'{accession}_{i+1}.fastq.gz' if layout=='PAIRED' else f'{accession}.fastq.gz'); dest=fq/name
        source=_ena_aspera_source(asperas[i]) if chosen=='aspera' else ''
        expected_md5=md5s[i] if i<len(md5s) else ''
        downloaded+=_download_ena_file_verified(
            accession,name,https_url,source,dest,fq,retries,verify_md5,expected_md5,
            chosen,ascp,ascp_key,ascp_port,ascp_rate,https_fallback
        )
        outputs.append(dest)
    if layout=='PAIRED':
        outputs=sorted(outputs,key=lambda p:_pair_rank(p.name))
        if len(outputs)!=2 or _pair_rank(outputs[0].name)!=1 or _pair_rank(outputs[1].name)!=2: raise RuntimeError(f'Could not unambiguously order ENA paired FASTQs for {accession}: {[p.name for p in outputs]}')
        r1,r2=outputs[0],outputs[1]; fmt=f'FASTQ_GZ_{chosen.upper()}'
    else:
        r1,r2=outputs[0],None; fmt=f'FASTQ_GZ_{chosen.upper()}_SINGLE'
    md['_ena_transport']=chosen.upper(); a=AcquiredReads(accession,r1,r2,md,'ENA',fmt,downloaded,[fq],None,0,'ENA_ARCHIVE_GENERATED_FASTQ')
    _save_manifest(work,a); return a

