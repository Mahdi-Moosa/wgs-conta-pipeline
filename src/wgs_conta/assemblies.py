from __future__ import annotations
import csv,fcntl,gzip,hashlib,json,re,shutil,subprocess,time,threading,urllib.parse,zipfile
from dataclasses import dataclass
from pathlib import Path
from .acquisition import _download_aspera,_ena_aspera_source,_ena_https
from .util import eprint,md5_file,run

# NCBI E-utilities allow 3 requests/second per IP without an API key.
# Multiple assembly-prefetch threads share this conservative process-wide
# pacer so concurrency does not become abusive request burstiness.
_NCBI_EUTILS_LOCK=threading.Lock()
_NCBI_EUTILS_LAST=0.0
_NCBI_EUTILS_MIN_INTERVAL=0.36  # ~2.78 request starts/second

def _pace_ncbi_eutils():
    global _NCBI_EUTILS_LAST
    with _NCBI_EUTILS_LOCK:
        now=time.monotonic()
        wait=_NCBI_EUTILS_MIN_INTERVAL-(now-_NCBI_EUTILS_LAST)
        if wait>0: time.sleep(wait)
        _NCBI_EUTILS_LAST=time.monotonic()

@dataclass
class AssemblyResolution:
    accession:str=''
    source:str='NONE'          # CURATED | NCBI | ENA | NONE
    lookup_status:str='NOT_FOUND'
    ena_ftp:str=''
    ena_aspera:str=''
    ena_md5:str=''
    lookup_complete:bool=True

def lookup_curated_run_assembly(path,accession):
    path=Path(path)
    if not path.exists(): return ''
    with path.open() as f:
        for row in csv.reader(f,delimiter='\t'):
            if row and row[0]==accession:
                gcfs=[x for x in row if re.fullmatch(r'GCF_\d+\.\d+',x or '')]; gcas=[x for x in row if re.fullmatch(r'GCA_\d+\.\d+',x or '')]
                return gcfs[0] if gcfs else (gcas[0] if gcas else '')
    return ''

def lookup_curated_fasta(path,assembly):
    path=Path(path)
    if not path.exists(): return None
    with path.open() as f:
        for row in csv.reader(f,delimiter='\t'):
            if row and row[0]==assembly and len(row)>=2:
                p=Path(row[1])
                if p.exists() and p.stat().st_size>0: return p
    return None

def _splitsemi(x): return [v for v in (x or '').split(';') if v]
def _is_fasta_name(x): return bool(re.search(r'\.(?:fa|fna|fasta)(?:\.gz)?$',x,re.I))

def _ena_search_request(api,params,retries,label):
    """Query ENA with retries only for technical request failures.

    A successful response containing zero rows is authoritative and is returned
    immediately; callers decide whether that means "no record" or whether a
    different ENA result type should be queried.
    """
    url=api.rstrip('/')+'?'+urllib.parse.urlencode(params)
    last=''
    for attempt in range(1,retries+1):
        eprint(f'+ [ENA {label} attempt {attempt}/{retries}]')
        p=subprocess.run(['curl','-fsSL','--connect-timeout','30',url],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if p.returncode==0: return p.stdout
        last=p.stderr[-2000:]
        if attempt<retries: time.sleep(min(4*attempt,20))
    raise RuntimeError(f'ENA query failed after {retries} attempts for {label}: {last}')


def _ena_analysis_lookup(biosample,api,retries=4):
    """Resolve an exact-BioSample assembly from ENA.

    First use ENA's assembly result (normally returns the INSDC GCA accession).
    If that successful query has no exact hit, fall back immediately to a
    SEQUENCE_ASSEMBLY analysis with a submitted FASTA. Both routes are
    constrained to the exact BioSample.
    """
    # Preferred ENA path: assembly database -> GCA accession.
    txt=_ena_search_request(api,{
        'result':'assembly','query':f'sample_accession="{biosample}"',
        'fields':'accession,sample_accession','format':'tsv','limit':'100'
    },retries,f'assembly-result {biosample}')
    rows=list(csv.DictReader(txt.splitlines(),delimiter='\t'))
    accs=[]
    for r in rows:
        if (r.get('sample_accession') or '')!=biosample: continue
        x=(r.get('accession') or '').strip()
        if re.fullmatch(r'GCA_\d+\.\d+',x): accs.append(x)
    if accs:
        return AssemblyResolution(accession=accs[0],source='ENA',lookup_status='ENA_EXACT_BIOSAMPLE_GCA')

    # Secondary ENA path: assembly analysis record -> submitted FASTA file.
    fields='accession,analysis_type,sample_accession,submitted_ftp,submitted_aspera,submitted_md5'
    query=f'sample_accession="{biosample}" AND analysis_type="SEQUENCE_ASSEMBLY"'
    txt=_ena_search_request(api,{'result':'analysis','query':query,'fields':fields,'format':'tsv','limit':'100'},retries,f'analysis-result {biosample}')
    rows=list(csv.DictReader(txt.splitlines(),delimiter='\t'))
    for r in rows:
        if (r.get('sample_accession') or '')!=biosample: continue
        if 'ASSEMBLY' not in (r.get('analysis_type') or '').upper(): continue
        ftps=_splitsemi(r.get('submitted_ftp','')); asperas=_splitsemi(r.get('submitted_aspera','')); md5s=_splitsemi(r.get('submitted_md5',''))
        for i,f in enumerate(ftps):
            name=Path(urllib.parse.urlparse(_ena_https(f)).path).name
            if _is_fasta_name(name):
                return AssemblyResolution(
                    accession=r.get('accession','') or 'ENA_SEQUENCE_ASSEMBLY', source='ENA', lookup_status='ENA_EXACT_BIOSAMPLE_ANALYSIS',
                    ena_ftp=f, ena_aspera=(asperas[i] if i<len(asperas) else ''), ena_md5=(md5s[i] if i<len(md5s) else '')
                )
    return AssemblyResolution(source='NONE',lookup_status='ENA_NO_EXACT_ASSEMBLY')


def _ncbi_eutils_json(endpoint,params,retries,label):
    """Direct NCBI EUtils GET returning decoded JSON.

    Retries are reserved for technical failures: curl/HTTP failure, empty or
    malformed JSON, or an explicit EUtils error payload. A valid biological
    result with count=0 is *not* retried.
    """
    base=f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}'
    params=dict(params); params.setdefault('tool','wgs_conta')
    last=''
    for attempt in range(1,retries+1):
        _pace_ncbi_eutils()
        eprint(f'+ [NCBI EUtils {label} attempt {attempt}/{retries}]')
        cmd=['curl','--http1.1','-fsSL','--connect-timeout','30','--get']
        for k,v in params.items():
            cmd += ['--data-urlencode',f'{k}={v}']
        cmd.append(base)
        p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if p.returncode==0:
            try:
                data=json.loads(p.stdout)
                if not isinstance(data,dict): raise ValueError('JSON root is not an object')
                if data.get('error'): raise ValueError(f"EUtils error: {data.get('error')}")
                return data
            except Exception as exc:
                last=f'invalid EUtils JSON: {exc}; body={p.stdout[:1000]!r}'
        else:
            last=p.stderr[-3000:]
        if attempt<retries: time.sleep(min(4*attempt,20))
    raise RuntimeError(f'NCBI EUtils {label} failed after {retries} attempts: {last}')


def _ncbi_assembly_lookup(biosample,retries=4):
    """Resolve exact-BioSample NCBI Assembly accession via direct EUtils JSON.

    The successful no-hit state is represented explicitly as NCBI_NO_RECORD and
    triggers immediate ENA fallback. It is not a retry condition.
    """
    search=_ncbi_eutils_json('esearch.fcgi',{
        'db':'assembly','term':f'{biosample}[BioSample]','retmode':'json','retmax':'100'
    },retries,f'esearch {biosample}')
    er=(search.get('esearchresult') or {})
    try: count=int(er.get('count','0') or 0)
    except Exception: raise RuntimeError(f'NCBI EUtils malformed count for {biosample}: {er.get("count")!r}')
    if count==0:
        return AssemblyResolution(source='NONE',lookup_status='NCBI_NO_RECORD')
    ids=[str(x) for x in (er.get('idlist') or []) if str(x)]
    if not ids:
        raise RuntimeError(f'NCBI EUtils returned count={count} but no assembly UIDs for {biosample}')

    summary=_ncbi_eutils_json('esummary.fcgi',{
        'db':'assembly','id':','.join(ids),'retmode':'json'
    },retries,f'esummary {biosample}')
    result=summary.get('result') or {}
    uids=[str(x) for x in (result.get('uids') or ids)]
    accs=[]
    for uid in uids:
        rec=result.get(uid) or {}
        # ESummary Assembly JSON uses assemblyaccession; tolerate a few
        # spelling variants defensively without weakening the exact-sample query.
        x=(rec.get('assemblyaccession') or rec.get('assembly_accession') or rec.get('assemblyacc') or '').strip()
        if re.fullmatch(r'GC[AF]_\d+\.\d+',x): accs.append(x)
    gcfs=[x for x in accs if x.startswith('GCF_')]
    gcas=[x for x in accs if x.startswith('GCA_')]
    x=gcfs[0] if gcfs else (gcas[0] if gcas else '')
    if x:
        return AssemblyResolution(accession=x,source='NCBI',lookup_status='NCBI_EXACT_BIOSAMPLE')
    # Search genuinely found Assembly records but no usable accession was in the
    # summary payload. This is a successful archive response, so fall back to ENA
    # rather than pointlessly repeating the same biological query.
    return AssemblyResolution(source='NONE',lookup_status='NCBI_RECORD_WITHOUT_ACCESSION')


def resolve_exact_assembly(accession,biosample,curated_run_assembly,retries=4,ena_search_api='https://www.ebi.ac.uk/ena/portal/api/search'):
    x=lookup_curated_run_assembly(curated_run_assembly,accession)
    if x: return AssemblyResolution(x,'CURATED','CURATED_EXACT_RUN')
    if not biosample: return AssemblyResolution(source='NONE',lookup_status='NO_BIOSAMPLE')

    try:
        nr=_ncbi_assembly_lookup(biosample,retries)
        if nr.source=='NCBI': return nr
        ncbi_status=nr.lookup_status
        eprint(f'+ [assembly fallback] {accession}/{biosample}: {ncbi_status}; trying ENA')
    except Exception as exc:
        ncbi_status='NCBI_LOOKUP_FAILED'
        eprint(f'WARNING: NCBI assembly lookup exhausted for {accession}/{biosample}: {exc}; trying ENA')

    try:
        r=_ena_analysis_lookup(biosample,ena_search_api,retries)
        if r.source=='ENA':
            r.lookup_status=f'{ncbi_status}__{r.lookup_status}'
            return r
        # If NCBI itself failed as a service, ENA saying 'no exact assembly'
        # cannot prove that no exact NCBI assembly exists. Keep the lookup
        # incomplete so the analysis is retained for a later retry.
        complete = (ncbi_status != 'NCBI_LOOKUP_FAILED')
        return AssemblyResolution(source='NONE',lookup_status=f'{ncbi_status}__{r.lookup_status}',lookup_complete=complete)
    except Exception as exc:
        eprint(f'WARNING: ENA assembly lookup also failed for {accession}/{biosample}: {exc}')
        return AssemblyResolution(source='NONE',lookup_status=f'{ncbi_status}__ENA_LOOKUP_FAILED',lookup_complete=False)

def _download_ncbi_assembly(assembly,root,curated_paths):
    cur=lookup_curated_fasta(curated_paths,assembly)
    if cur: return cur
    root=Path(root)/assembly; root.mkdir(parents=True,exist_ok=True); locks=root.parent/'.locks'; locks.mkdir(exist_ok=True)
    with (locks/f'{assembly}.lock').open('w') as lf:
        fcntl.flock(lf,fcntl.LOCK_EX)
        hits=list(root.rglob('*_genomic.fna'))
        if hits: return hits[0]
        z=root/'dataset.zip'
        try:
            run(['datasets','download','genome','accession',assembly,'--include','genome','--filename',str(z),'--no-progressbar'])
            with zipfile.ZipFile(z) as zh:
                bad=zh.testzip()
                if bad: raise RuntimeError(f'invalid zip member {bad}')
                zh.extractall(root/'package')
        except Exception as exc:
            eprint(f'WARNING: direct Datasets package failed for {assembly}: {exc}; using dehydrated mode')
            shutil.rmtree(root/'package',ignore_errors=True); dz=root/'dehydrated.zip'
            run(['datasets','download','genome','accession',assembly,'--include','genome','--dehydrated','--filename',str(dz),'--no-progressbar'])
            dd=root/'dehydrated'; shutil.rmtree(dd,ignore_errors=True)
            with zipfile.ZipFile(dz) as zh: zh.extractall(dd)
            run(['datasets','rehydrate','--directory',str(dd),'--no-progressbar'])
        hits=list(root.rglob('*_genomic.fna'))
        if not hits: raise RuntimeError(f'No genomic FASTA found for {assembly}')
        return hits[0]

def _download_ena_assembly(r:AssemblyResolution,root,ascp,key,port,rate,retries):
    """Download an exact ENA assembly with integrity-aware retry.

    Transport failures retain normal Aspera retry/resume behavior.  A completed
    transfer whose MD5 disagrees with ENA is known-bad: delete only that source
    file and reacquire it from scratch.  Existing validated assembly.fna files
    from earlier releases remain reusable because v2.4.1 only created them after
    a successful source-file checksum.
    """
    d=Path(root)/f'ENA_{r.accession}'; d.mkdir(parents=True,exist_ok=True)
    locks=Path(root)/'.locks'; locks.mkdir(parents=True,exist_ok=True)
    with (locks/f'ENA_{r.accession}.lock').open('w') as lf:
        fcntl.flock(lf,fcntl.LOCK_EX)

        final=d/'assembly.fna'
        if final.exists() and final.stat().st_size>0:
            eprint(f'+ [assembly cache] reusing {final}')
            return final

        url=_ena_https(r.ena_ftp)
        srcname=Path(urllib.parse.urlparse(url).path).name
        raw=d/srcname
        expected=(r.ena_md5 or '').strip().lower()

        # A failed v2.4.0/v2.4.1 attempt may have left an unvalidated source
        # file.  Verify before reuse; delete immediately if known-corrupt.
        if raw.exists() and raw.stat().st_size>0 and expected:
            got=md5_file(raw).lower()
            if got==expected:
                eprint(f'+ [assembly resume] existing ENA source MD5 PASS {raw.name}; reusing file')
            else:
                eprint(f'+ [assembly resume] existing ENA source MD5 FAIL {raw.name}: got {got}, expected {expected}; removing before fresh acquisition')
                raw.unlink(missing_ok=True)

        integrity_attempts=retries if expected else 1
        for integrity_attempt in range(1,integrity_attempts+1):
            if not raw.exists() or raw.stat().st_size==0:
                if r.ena_aspera and Path(ascp).expanduser().exists() and Path(key).expanduser().exists():
                    try:
                        _download_aspera(_ena_aspera_source(r.ena_aspera),d,ascp,key,port,rate,retries,srcname)
                    except Exception as exc:
                        # Assembly files are small relative to reads; a direct
                        # HTTPS fallback is preferable to turning a transient
                        # FASP failure into an analysis error.
                        eprint(f'WARNING: ENA assembly Aspera failed for {r.accession}: {exc}; falling back to HTTPS')
                        raw.unlink(missing_ok=True)
                        run(['curl','-fL','--retry',str(retries),'--retry-all-errors','--connect-timeout','30','-o',str(raw),url])
                else:
                    run(['curl','-fL','--retry',str(retries),'--retry-all-errors','--connect-timeout','30','-o',str(raw),url])

            if not raw.exists() or raw.stat().st_size==0:
                raise RuntimeError(f'ENA assembly transfer completed but source file is missing/empty: {raw}')

            if expected:
                got=md5_file(raw).lower()
                if got!=expected:
                    eprint(f'+ ASSEMBLY MD5 FAIL {raw.name}: got {got}, expected {expected}')
                    raw.unlink(missing_ok=True)
                    if integrity_attempt<integrity_attempts:
                        eprint(f'+ removed corrupt assembly source {raw.name}; fresh acquisition attempt {integrity_attempt+1}/{integrity_attempts}')
                        time.sleep(min(2*integrity_attempt,8))
                        continue
                    raise RuntimeError(
                        f'ENA assembly MD5 mismatch after {integrity_attempts} fresh acquisition attempt(s): '
                        f'{got} != {expected}; corrupt source removed'
                    )
                eprint(f'+ ASSEMBLY MD5 PASS {raw.name}')

            tmp=d/'assembly.fna.tmp'
            tmp.unlink(missing_ok=True)
            if raw.suffix=='.gz':
                with gzip.open(raw,'rb') as fi,tmp.open('wb') as fo:
                    shutil.copyfileobj(fi,fo)
            else:
                shutil.copy2(raw,tmp)
            if not tmp.exists() or tmp.stat().st_size==0:
                tmp.unlink(missing_ok=True)
                raise RuntimeError(f'ENA assembly decompression/copy produced empty FASTA: {raw}')
            tmp.replace(final)
            return final

        raise RuntimeError(f'Unreachable ENA assembly integrity retry state for {r.accession}')

def download_assembly(resolution:AssemblyResolution,root,curated_paths,ena_ascp,ena_key,ena_port,ena_rate,retries):
    if not resolution.accession: raise RuntimeError('No assembly resolution available')
    if resolution.source in {'CURATED','NCBI'}:
        return _download_ncbi_assembly(resolution.accession,root,curated_paths)
    if resolution.source=='ENA':
        if re.fullmatch(r'GCA_\d+\.\d+',resolution.accession) and not resolution.ena_ftp:
            return _download_ncbi_assembly(resolution.accession,root,curated_paths)
        return _download_ena_assembly(resolution,root,ena_ascp,ena_key,ena_port,ena_rate,retries)
    raise RuntimeError(f'Unsupported assembly source: {resolution.source}')
