from __future__ import annotations
import csv, json, sqlite3, threading
from pathlib import Path

SUMMARY_COLUMNS=[
'Accession','Label','Analysis_branch','Quantification_unit','Acquisition_source','Acquisition_format','Acquisition_bytes','Read_filter_policy','SRA_reject_units_removed','SRA_reject_pairs_removed','SRA_reject_reads_removed','Metadata_source','BioProject','BioSample','Organism','Platform','Layout',
'Assembly','Assembly_source','Assembly_lookup_status',
'Total_reads','Total_pairs','Total_bases','Mean_read_length','Q20_pct','Q30_pct','GC_pct','MAPQ30_reads','Cov80_ID95','Cov80_ID98','Cov80_ID99','Strict_pairs_ID98','Strict_reads_ID98','Unique_strict_fragments','Unique_strict_reads','Unique_strict_units','Strict_per_million_pairs','Strict_per_million_reads','Strict_per_million_units',
'Kraken_strict_human','Kraken_human_fraction','N_chromosomes','Chromosome_distribution','Top_chromosome','Top_count','Top_fraction','Bact_ge1_ID95','Bact_both_ID95','Bact_ge1_ID98','Bact_both_ID98','Human_positive','Assembly_check','Runtime_seconds','Completed_at']

TERMINAL={'DONE','ERROR_ACQUISITION','ERROR_ANALYSIS','SKIPPED_UNSUPPORTED'}

def connect(path: Path):
    path.parent.mkdir(parents=True,exist_ok=True)
    c=sqlite3.connect(path,timeout=60,check_same_thread=False); c.row_factory=sqlite3.Row
    c.execute('PRAGMA journal_mode=WAL'); c.execute('PRAGMA synchronous=NORMAL'); return c

def _columns(c): return {r['name'] for r in c.execute('PRAGMA table_info(jobs)').fetchall()}
def _addcol(c,name,decl):
    if name not in _columns(c): c.execute(f'ALTER TABLE jobs ADD COLUMN {name} {decl}')

def init_db(path: Path):
    c=connect(path)
    c.execute("""CREATE TABLE IF NOT EXISTS jobs(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        accession TEXT NOT NULL UNIQUE,
        label TEXT NOT NULL DEFAULT '',
        status TEXT NOT NULL DEFAULT 'PENDING',
        backend TEXT, worker TEXT, attempts INTEGER NOT NULL DEFAULT 0,
        started_at TEXT, completed_at TEXT, error TEXT, result_json TEXT,
        aws_exhausted INTEGER NOT NULL DEFAULT 0,
        ncbi_exhausted INTEGER NOT NULL DEFAULT 0,
        ena_exhausted INTEGER NOT NULL DEFAULT 0,
        aws_error TEXT, ncbi_error TEXT, ena_error TEXT,
        resume_backend TEXT,
        preferred_backend TEXT,
        priority INTEGER NOT NULL DEFAULT 100,
        pipeline_version TEXT,
        selected INTEGER NOT NULL DEFAULT 1
    )""")
    # In-place migration from older state DBs.
    _addcol(c,'aws_exhausted','INTEGER NOT NULL DEFAULT 0')
    _addcol(c,'ncbi_exhausted','INTEGER NOT NULL DEFAULT 0')
    _addcol(c,'ena_exhausted','INTEGER NOT NULL DEFAULT 0')
    _addcol(c,'aws_error','TEXT')
    _addcol(c,'ncbi_error','TEXT')
    _addcol(c,'ena_error','TEXT')
    _addcol(c,'resume_backend','TEXT')
    _addcol(c,'preferred_backend','TEXT')
    _addcol(c,'priority','INTEGER NOT NULL DEFAULT 100')
    _addcol(c,'pipeline_version','TEXT')
    _addcol(c,'selected','INTEGER NOT NULL DEFAULT 1')
    c.execute('CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status,priority,id)')
    c.commit(); c.close()

def import_jobs(path: Path, rows, force=False):
    """Legacy batch-mode importer retained for run_pipeline.py."""
    c=connect(path)
    with c:
        for acc,label in rows:
            c.execute("INSERT OR IGNORE INTO jobs(accession,label,status,priority,selected) VALUES(?,?,'PENDING',100,1)",(acc,label))
            c.execute('UPDATE jobs SET selected=1 WHERE accession=?',(acc,))
            if label: c.execute('UPDATE jobs SET label=? WHERE accession=?',(label,acc))
            if force:
                c.execute("""UPDATE jobs SET status='PENDING',backend=NULL,worker=NULL,error=NULL,result_json=NULL,
                          started_at=NULL,completed_at=NULL,aws_exhausted=0,ncbi_exhausted=0,ena_exhausted=0,aws_error=NULL,ncbi_error=NULL,ena_error=NULL,
                          resume_backend=NULL,preferred_backend=NULL,priority=100,selected=1
                          WHERE accession=?""",(acc,))
    c.close()

def prepare_master_jobs(path:Path,rows,pipeline_version):
    """Synchronize scientifically incomplete registry accessions into runtime.

    The permanent master is authoritative. Any selected accession that the
    master says is incomplete is made retryable automatically, regardless of
    stale runtime states from a prior controller invocation. This makes the
    SQLite DB disposable execution state rather than something users must
    manually rename/reset between runs.

    rows: (accession,label,preferred_backend,priority)
    """
    c=connect(path)
    with c:
        c.execute('UPDATE jobs SET selected=0')
        for acc,label,pref,priority in rows:
            old=c.execute('SELECT pipeline_version,status,backend,resume_backend FROM jobs WHERE accession=?',(acc,)).fetchone()
            if old is None:
                c.execute("""INSERT INTO jobs(accession,label,status,preferred_backend,priority,pipeline_version,selected)
                             VALUES(?,?,'PENDING',?,?,?,1)""",(acc,label,pref or None,int(priority),pipeline_version))
                continue

            version_changed=(old['pipeline_version'] or '')!=pipeline_version
            status=old['status'] or 'PENDING'
            if status=='ERROR_ANALYSIS':
                # Acquisition succeeded previously. Preserve the source pin and
                # retained normalized reads even across an orchestration-version
                # upgrade; acquisition.py itself decides whether a cached
                # normalization manifest remains scientifically compatible.
                resume=(old['resume_backend'] or old['backend'] or None)
                c.execute("""UPDATE jobs SET label=?,status='PENDING',worker=NULL,error=NULL,completed_at=NULL,
                             resume_backend=?,preferred_backend=?,priority=?,pipeline_version=?,selected=1 WHERE accession=?""",
                          (label,resume,pref or None,int(priority),pipeline_version,acc))
            elif status in {'RUNNING','RUNNING_ACQUIRE','RUNNING_ANALYSIS'}:
                # Controller lock guarantees there is no live peer using this DB.
                # Resume on the same source so retained work can be reused.
                resume=(old['backend'] or old['resume_backend'] or None)
                c.execute("""UPDATE jobs SET label=?,status='PENDING',worker=NULL,error=NULL,completed_at=NULL,
                             resume_backend=?,preferred_backend=?,priority=?,pipeline_version=?,selected=1 WHERE accession=?""",
                          (label,resume,pref or None,int(priority),pipeline_version,acc))
            elif status=='ERROR_ACQUISITION':
                # Both archive attempts may have been exhausted in a previous
                # invocation; a fresh controller gets a fresh retry budget.
                c.execute("""UPDATE jobs SET label=?,status='PENDING',backend=NULL,worker=NULL,error=NULL,completed_at=NULL,
                             aws_exhausted=0,ncbi_exhausted=0,ena_exhausted=0,aws_error=NULL,ncbi_error=NULL,ena_error=NULL,resume_backend=NULL,
                             preferred_backend=?,priority=?,pipeline_version=?,selected=1 WHERE accession=?""",
                          (label,pref or None,int(priority),pipeline_version,acc))
            elif version_changed or status in {'DONE','RESULT_COMMITTED'}:
                # New orchestration version or stale completion while master is
                # incomplete: rebuild runtime state from scratch.
                c.execute("""UPDATE jobs SET label=?,status='PENDING',backend=NULL,worker=NULL,attempts=0,
                             started_at=NULL,completed_at=NULL,error=NULL,result_json=NULL,aws_exhausted=0,ncbi_exhausted=0,ena_exhausted=0,
                             aws_error=NULL,ncbi_error=NULL,ena_error=NULL,resume_backend=NULL,preferred_backend=?,priority=?,pipeline_version=?,selected=1
                             WHERE accession=?""",(label,pref or None,int(priority),pipeline_version,acc))
            else:
                # PENDING and SKIPPED_UNSUPPORTED retain their semantic state.
                c.execute('UPDATE jobs SET label=?,preferred_backend=?,priority=?,pipeline_version=?,selected=1 WHERE accession=?',
                          (label,pref or None,int(priority),pipeline_version,acc))
    c.close()

def reset_running(path: Path):
    c=connect(path)
    c.execute("""UPDATE jobs SET status='PENDING',resume_backend=backend,worker=NULL
                 WHERE status IN ('RUNNING','RUNNING_ACQUIRE','RUNNING_ANALYSIS')""")
    c.commit(); c.close()

def reset_errors(path: Path):
    c=connect(path)
    c.execute("UPDATE jobs SET status='PENDING',error=NULL,worker=NULL WHERE status='ERROR_ANALYSIS'")
    c.execute("""UPDATE jobs SET status='PENDING',backend=NULL,worker=NULL,error=NULL,completed_at=NULL,
                 aws_exhausted=0,ncbi_exhausted=0,ena_exhausted=0,aws_error=NULL,ncbi_error=NULL,ena_error=NULL,resume_backend=NULL
                 WHERE status='ERROR_ACQUISITION'""")
    c.commit(); c.close()

def _backend_cols(backend):
    if backend not in {'aws','ncbi','ena'}: raise ValueError(f'Unknown backend: {backend}')
    return f'{backend}_exhausted',f'{backend}_error'


def claim_job(path: Path, backend, worker):
    col,_=_backend_cols(backend)
    c=connect(path)
    try:
        c.execute('BEGIN IMMEDIATE')
        r=c.execute(f"""SELECT id,accession,label FROM jobs
                         WHERE selected=1 AND status='PENDING' AND {col}=0
                           AND (resume_backend IS NULL OR resume_backend=? )
                           AND (preferred_backend IS NULL OR preferred_backend=? )
                         ORDER BY priority,id LIMIT 1""",(backend,backend)).fetchone()
        if r is None: c.commit(); return None
        c.execute("""UPDATE jobs SET status='RUNNING_ACQUIRE',backend=?,worker=?,attempts=attempts+1,
                     started_at=datetime('now'),error=NULL WHERE id=?""",(backend,worker,r['id']))
        c.commit(); return dict(r)
    finally: c.close()

def mark_stage(path: Path, acc, status):
    c=connect(path); c.execute('UPDATE jobs SET status=? WHERE accession=?',(status,acc)); c.commit(); c.close()

def mark_acquisition_exhausted(path: Path, acc, backend, error, available_backends):
    """Exhaust one backend and requeue for any remaining acquisition source."""
    available={b for b in available_backends if b in {'aws','ncbi','ena'}}
    own,errcol=_backend_cols(backend)
    c=connect(path)
    try:
        c.execute('BEGIN IMMEDIATE')
        # Clear source preference once the preferred backend has exhausted its
        # retry budget, allowing any remaining configured source to claim it.
        c.execute(f"UPDATE jobs SET {own}=1,{errcol}=?,error=?,resume_backend=NULL,preferred_backend=NULL,worker=NULL WHERE accession=?",
                  (error[:10000],error[:10000],acc))
        r=c.execute('SELECT aws_exhausted,ncbi_exhausted,ena_exhausted FROM jobs WHERE accession=?',(acc,)).fetchone()
        can_other=False
        if r is not None:
            for other in available-{backend}:
                if int(r[f'{other}_exhausted'] or 0)==0:
                    can_other=True; break
        status='PENDING' if can_other else 'ERROR_ACQUISITION'
        if status=='ERROR_ACQUISITION':
            c.execute("UPDATE jobs SET status=?,completed_at=datetime('now') WHERE accession=?",(status,acc))
        else:
            c.execute('UPDATE jobs SET status=?,completed_at=NULL WHERE accession=?',(status,acc))
        c.commit(); return status
    finally: c.close()

def mark_unsupported(path: Path, acc, error):
    c=connect(path)
    c.execute("""UPDATE jobs SET status='SKIPPED_UNSUPPORTED',error=?,worker=NULL,resume_backend=NULL,completed_at=datetime('now')
                 WHERE accession=?""",(error[:10000],acc)); c.commit(); c.close()

def mark_analysis_error(path: Path, acc, backend, error):
    c=connect(path)
    c.execute("""UPDATE jobs SET status='ERROR_ANALYSIS',error=?,resume_backend=?,worker=NULL,completed_at=datetime('now')
                 WHERE accession=?""",(error[:10000],backend,acc)); c.commit(); c.close()

def mark_result_committed(path: Path, acc, result):
    c=connect(path)
    c.execute("""UPDATE jobs SET status='RESULT_COMMITTED',result_json=?,completed_at=datetime('now'),error=NULL
                 WHERE accession=?""",(json.dumps(result,sort_keys=True),acc)); c.commit(); c.close()

def mark_done(path: Path, acc):
    c=connect(path); c.execute("UPDATE jobs SET status='DONE',error=NULL,resume_backend=NULL,worker=NULL WHERE accession=?",(acc,)); c.commit(); c.close()

def result_committed_jobs(path: Path):
    c=connect(path); rows=c.execute("SELECT accession,backend,result_json FROM jobs WHERE status='RESULT_COMMITTED' ORDER BY id").fetchall(); c.close(); return [dict(r) for r in rows]

def has_unfinished(path: Path):
    c=connect(path)
    n=c.execute("SELECT COUNT(*) n FROM jobs WHERE selected=1 AND status IN ('PENDING','RUNNING_ACQUIRE','RUNNING_ANALYSIS')").fetchone()['n']
    c.close(); return bool(n)

def status_counts(path: Path,selected_only=False):
    c=connect(path)
    q='SELECT status,COUNT(*) n FROM jobs' + (' WHERE selected=1' if selected_only else '') + ' GROUP BY status'
    rows=c.execute(q).fetchall(); c.close(); return {r['status']:r['n'] for r in rows}

_export_lock=threading.Lock()
def export_master(path: Path, out: Path):
    """Legacy per-batch TSV export. Master mode commits directly via MasterStore."""
    with _export_lock:
        c=connect(path)
        rows=c.execute("""SELECT result_json FROM jobs
                          WHERE status IN ('RESULT_COMMITTED','DONE') AND result_json IS NOT NULL ORDER BY id""").fetchall(); c.close()
        out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+'.tmp')
        with tmp.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=SUMMARY_COLUMNS,delimiter='\t'); w.writeheader()
            for r in rows:
                d=json.loads(r['result_json']); w.writerow({k:d.get(k,'NA') for k in SUMMARY_COLUMNS})
        tmp.replace(out)
