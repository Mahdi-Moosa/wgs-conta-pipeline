"""Bounded runtime checks around the unchanged scientific callers."""
import configparser,hashlib,re
from pathlib import Path

def file_sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def validate_runtime(s,require_resources=True,config_path=None):
    if min(s.ncbi_workers,s.ena_workers,s.aws_analysis_workers,s.aws_download_workers)<1:
        raise ValueError('backend worker counts must be positive')
    if min(s.threads_per_worker,s.sort_threads,s.kraken_threads,s.retries)<1:
        raise ValueError('threads and retry counts must be positive')
    installed=Path(__file__).resolve().parent
    outputs=[s.project_root,s.sra_root,s.metadata_root,s.results_root,s.evidence_root,
             s.logs_root,s.work_root,s.assembly_root,s.state_db.parent]
    for p in outputs:
        p=Path(p).resolve()
        if p==Path(p.anchor) or p==Path.home() or p==installed or installed.is_relative_to(p) or p.is_relative_to(installed):
            raise ValueError('unsafe output location')
        if s.human_index.is_relative_to(p) or s.kraken_db.is_relative_to(p):
            raise ValueError('reference resources must be outside output/cleanup roots')
    if require_resources:
        if not s.human_index.is_file():raise FileNotFoundError('human minimap2 index does not exist: '+str(s.human_index))
        if not s.kraken_db.is_dir():raise FileNotFoundError('Kraken2 database does not exist: '+str(s.kraken_db))
        if config_path is None:raise ValueError('configuration required for reference identity validation')
        c=configparser.ConfigParser();c.read(config_path)
        fasta=Path(c.get('reference_identity','human_fasta',fallback=''))
        if not fasta.is_absolute():fasta=Path(config_path).resolve().parent/fasta
        if not fasta.is_file():raise FileNotFoundError('CHM13 masked-Y/rCRS FASTA is required')
        expected_fasta='a5128da5a5729d1ef92bf2040494092d784a5ff07e4d08f7b9f385b5498d84d3'
        if file_sha256(fasta)!=expected_fasta:raise ValueError('FASTA does not match the validated CHM13 masked-Y/rCRS analysis set')
        expected=c.get('reference_identity','human_index_sha256',fallback='')
        if not re.fullmatch('[0-9a-f]{64}',expected):raise ValueError('set reference_identity.human_index_sha256 to your verified index hash')
        if file_sha256(s.human_index)!=expected:raise ValueError('human index checksum mismatch')
    return True


def guard_resume_identity(s,config_path,dry_run=False):
    """Reject implicit reuse across references, classifier identity or science code changes."""
    import json,subprocess,shutil
    from . import analysis,ion_analysis
    c=configparser.ConfigParser();c.read(config_path)
    identity_path=Path(c.get('reference_identity','kraken_identity_file',fallback=''))
    if not identity_path.is_absolute():identity_path=Path(config_path).resolve().parent/identity_path
    if not identity_path.is_file() and not dry_run:raise FileNotFoundError('provide a Kraken database version/checksum manifest')
    # The database manifest must identify/checksum the actual database; its verification
    # remains an external preparation requirement, not merely this fingerprint.
    versions={}
    if not dry_run:
        for name in ('minimap2','samtools','kraken2'):
            result=subprocess.run([name,'--version'],capture_output=True,text=True,check=True)
            versions[name]=(result.stdout or result.stderr).splitlines()[0]
    settings=dict(human_index_sha256=c.get('reference_identity','human_index_sha256',fallback=''),
                  kraken_database=str(s.kraken_db),kraken_identity_sha256=file_sha256(identity_path) if identity_path.is_file() else None,
                  science_modules={Path(m.__file__).name:file_sha256(m.__file__) for m in (analysis,ion_analysis)},versions=versions)
    marker=s.project_root/'master'/'release_scientific_identity.json'
    if marker.exists():
        old=json.loads(marker.read_text())
        if dry_run:
            old.pop('versions',None);settings.pop('versions',None)
        if old!=settings:raise ValueError('scientific configuration changed; use a new isolated workspace')
    elif (s.project_root/'master'/'master_results.tsv').exists():
        raise ValueError('legacy master lacks release identity; use a new isolated workspace')
    elif not dry_run:
        marker.parent.mkdir(parents=True,exist_ok=True)
        with marker.open('x') as f:json.dump(settings,f,indent=2)
