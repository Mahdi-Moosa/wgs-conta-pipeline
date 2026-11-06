from __future__ import annotations
import configparser
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class Settings:
    project_root: Path; sra_root: Path; metadata_root: Path; results_root: Path; evidence_root: Path; logs_root: Path; work_root: Path; assembly_root: Path
    state_db: Path; master_summary: Path; human_index: Path; kraken_db: Path; curated_run_assembly: Path; curated_assembly_paths: Path
    source_mode: str; ncbi_workers: int; ena_workers: int; aws_analysis_workers: int; aws_download_workers: int; aws_prefetch_depth: int; assembly_prefetch_workers: int
    threads_per_worker: int; sort_threads: int; kraken_threads: int; retries: int
    delete_new_sra_after_success: bool; ena_api: str; ena_search_api: str; ena_verify_md5: bool; ena_transport: str; ena_ascp: Path; ena_ascp_key: Path; ena_ascp_port: int; ena_ascp_rate: str; ena_https_fallback: bool; assembly_query_retries: int
    aws_cli: Path; aws_bucket: str

def _p(c,s,k):
    p=Path(c.get(s,k)).expanduser()
    return (c._release_base/p).resolve() if not p.is_absolute() else p.resolve()
def load_config(path: Path):
    c=configparser.ConfigParser()
    if not c.read(path): raise FileNotFoundError(path)
    c._release_base=Path(path).resolve().parent
    mode=c.get('acquisition','source_mode',fallback='both').strip().lower()
    if mode=='ebi': mode='ena'
    if mode not in {'ncbi','ena','both','aws','all'}:
        raise ValueError('source_mode must be ncbi, ena/ebi, both, aws, or all')
    transport=c.get('acquisition','ena_transport',fallback='aspera').strip().lower()
    if transport not in {'aspera','https','auto'}: raise ValueError('ena_transport must be aspera, https, or auto')
    return Settings(
        _p(c,'paths','project_root'),_p(c,'paths','sra_root'),_p(c,'paths','metadata_root'),_p(c,'paths','results_root'),_p(c,'paths','evidence_root'),_p(c,'paths','logs_root'),_p(c,'paths','work_root'),_p(c,'paths','assembly_root'),
        _p(c,'paths','state_db'),_p(c,'paths','master_summary'),_p(c,'references','human_index'),_p(c,'references','kraken_db'),_p(c,'references','curated_run_assembly'),_p(c,'references','curated_assembly_paths'),
        mode,
        c.getint('workers','ncbi_workers',fallback=1),
        c.getint('workers','ena_workers',fallback=1),
        c.getint('workers','aws_analysis_workers',fallback=2),
        c.getint('workers','aws_download_workers',fallback=4),
        c.getint('workers','aws_prefetch_depth',fallback=8),
        c.getint('workers','assembly_prefetch_workers',fallback=4),
        c.getint('workers','threads_per_worker',fallback=5),
        c.getint('workers','sort_threads',fallback=2),
        c.getint('workers','kraken_threads',fallback=5),
        c.getint('acquisition','retries',fallback=4),
        c.getboolean('cleanup','delete_new_sra_after_success',fallback=True),
        c.get('acquisition','ena_api',fallback='https://www.ebi.ac.uk/ena/portal/api/filereport'),
        c.get('assembly','ena_search_api',fallback='https://www.ebi.ac.uk/ena/portal/api/search'),
        c.getboolean('acquisition','ena_verify_md5',fallback=True),
        transport,
        _p(c,'acquisition','ena_ascp'),
        _p(c,'acquisition','ena_ascp_key'),
        c.getint('acquisition','ena_ascp_port',fallback=33001),
        c.get('acquisition','ena_ascp_rate',fallback='300m').strip(),
        c.getboolean('acquisition','ena_https_fallback',fallback=False),
        c.getint('assembly','query_retries',fallback=4),
        Path(__import__('shutil').which(c.get('acquisition','aws_cli',fallback='aws')) or c.get('acquisition','aws_cli',fallback='aws')).expanduser(),
        c.get('acquisition','aws_bucket',fallback='s3://sra-pub-run-odp/sra').rstrip('/'),
    )
