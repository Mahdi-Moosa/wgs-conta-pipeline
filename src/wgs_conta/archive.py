#!/usr/bin/env python3
"""Resumable, count-checked ENA read_run metadata capture. Standard library only."""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import shutil
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

FIELDS = ('run_accession study_accession secondary_study_accession sample_accession '
          'secondary_sample_accession scientific_name tax_id country strain isolate '
          'isolation_source center_name first_public instrument_platform instrument_model '
          'library_strategy library_source library_layout nominal_length read_count '
          'base_count sample_title study_title fastq_ftp fastq_md5 fastq_bytes '
          'submitted_ftp submitted_md5 submitted_bytes').split()
LAYOUT = {'ILLUMINA': 'PAIRED', 'BGISEQ': 'PAIRED',
          'DNBSEQ': 'PAIRED', 'ION_TORRENT': 'SINGLE'}
BASE = 'tax_tree(2) AND library_strategy="WGS" AND library_source="GENOMIC"'
PLAN_COLUMNS = ['id', 'platform', 'start', 'end', 'query']
MANIFEST_COLUMNS = PLAN_COLUMNS + ['fields', 'status', 'expected_before', 'expected_after',
    'rows', 'unique_runs', 'raw_bytes', 'gzip_bytes', 'raw_sha256', 'gzip_sha256',
    'utc_start', 'utc_end', 'retry_count', 'file']


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_text(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    with open(tmp, 'w', encoding='utf-8', newline='') as f:
        f.write(content)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def tsv_string(rows, columns):
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=columns, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue()


def date(s):
    return dt.date.fromisoformat(s)


def cutoff_exclusive(cfg):
    return (date(cfg['cutoff_date']) + dt.timedelta(days=1)).isoformat()


def intervals(start, end, unit='year'):
    """Contiguous [start, end) intervals, clipped to the requested bounds."""
    cur = date(start)
    last = date(end)
    if cur >= last:
        raise ValueError('start must precede exclusive end')
    while cur < last:
        if unit == 'year':
            nxt = dt.date(cur.year + 1, 1, 1)
        elif unit == 'month':
            nxt = dt.date(cur.year + (cur.month == 12), cur.month % 12 + 1, 1)
        elif unit == 'day':
            nxt = cur + dt.timedelta(days=1)
        else:
            raise ValueError(unit)
        stop = min(last, nxt)
        yield cur.isoformat(), stop.isoformat()
        cur = stop


def platform_query(platform):
    return f'{BASE} AND library_layout="{LAYOUT[platform]}" AND instrument_platform="{platform}"'


def range_query(platform, start, end):
    return f'{platform_query(platform)} AND first_public>={start} AND first_public<{end}'


def union_query():
    return (BASE + ' AND ((library_layout="PAIRED" AND '
            '(instrument_platform="ILLUMINA" OR instrument_platform="BGISEQ" OR '
            'instrument_platform="DNBSEQ")) OR (library_layout="SINGLE" AND '
            'instrument_platform="ION_TORRENT"))')


def make_partition(platform, start, end):
    return dict(id=f'{platform}_{start}_{end}', platform=platform, start=start,
                end=end, query=range_query(platform, start, end))


def initial_plan(cfg):
    return [make_partition(p, a, b) for p in cfg['platforms']
            for a, b in intervals(cfg['start_date'], cutoff_exclusive(cfg))]


def validate_plan(plan, cfg):
    for platform in cfg['platforms']:
        rows = sorted((r for r in plan if r['platform'] == platform), key=lambda x: x['start'])
        if not rows or rows[0]['start'] != cfg['start_date'] or rows[-1]['end'] != cutoff_exclusive(cfg):
            raise ValueError(f'plan boundary gap: {platform}')
        prev = None
        for r in rows:
            if r['id'] != make_partition(platform, r['start'], r['end'])['id'] or r['query'] != range_query(platform, r['start'], r['end']):
                raise ValueError(f'plan query mismatch: {r}')
            if date(r['start']) >= date(r['end']) or (prev and r['start'] != prev['end']):
                raise ValueError(f'plan overlap/gap: {platform}')
            prev = r
    if set(r['platform'] for r in plan) != set(cfg['platforms']):
        raise ValueError('plan platform mismatch')


class ENA:
    def __init__(self, cfg):
        self.root = cfg['api_base'].rstrip('/')
        self.connect_timeout = cfg['connect_timeout_seconds']
        self.transfer_timeout = cfg['transfer_timeout_seconds']
        self.agent = 'wgs-conta-ena-metadata-snapshot/1.0 (research metadata only)'
        self.http_attempts = cfg['max_attempts']
        self.retry_base = cfg['retry_base_seconds']

    def get(self, endpoint, params, target=None):
        url = self.root + '/' + endpoint + '?' + urllib.parse.urlencode(params)
        for attempt in range(self.http_attempts):
            req = urllib.request.Request(url, headers={'User-Agent': self.agent, 'Accept-Encoding': 'identity'})
            start = time.monotonic()
            try:
                with urllib.request.urlopen(req, timeout=self.connect_timeout) as response:
                    if response.status != 200:
                        raise RuntimeError(f'HTTP {response.status}: {url}')
                    if response.headers.get('Content-Encoding', 'identity').lower() != 'identity':
                        raise RuntimeError('unexpected HTTP content encoding')
                    if target is None:
                        return response.read(20_000_000)
                    with open(target, 'wb') as f:
                        while True:
                            if time.monotonic() - start > self.transfer_timeout:
                                raise TimeoutError('transfer deadline exceeded')
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            f.write(block)
                        f.flush()
                        os.fsync(f.fileno())
                    return None
            except (OSError, TimeoutError, urllib.error.URLError) as e:
                if attempt + 1 == self.http_attempts:
                    raise
                record = {'status': 'NETWORK_RETRY', 'endpoint': endpoint,
                          'attempt': attempt + 1, 'error': str(e), 'utc': utc()}
                if getattr(self, 'log_root', None):
                    log_attempt(self.log_root, record)
                else:
                    print(json.dumps(record), flush=True)
                time.sleep(min(self.retry_base * 2 ** attempt, 120))

    def count(self, query):
        data = self.get('count', {'result': 'read_run', 'query': query}).decode('utf-8-sig').strip().splitlines()
        if len(data) != 2 or data[0] != 'count' or not data[1].isdigit():
            raise ValueError(f'invalid ENA count response: {data[:3]}')
        return int(data[1])

    def fields(self):
        data = self.get('returnFields', {'result': 'read_run'}).decode('utf-8-sig')
        rows = csv.DictReader(io.StringIO(data), delimiter='\t')
        names = {r.get('columnId', '') for r in rows}
        missing = set(FIELDS) - names
        if missing:
            raise ValueError(f'ENA returnFields missing: {sorted(missing)}; headers={rows.fieldnames}')
        return hashlib.sha256(data.encode()).hexdigest()

    def download(self, query, path):
        self.get('search', {'result': 'read_run', 'query': query, 'fields': ','.join(FIELDS),
                            'format': 'tsv', 'limit': '0'}, path)


def validate_tsv(path, partition, cutoff):
    n = 0
    seen = set()
    with open(path, 'rb') as raw:
        prefix = raw.read(512).lstrip().lower()
        if prefix.startswith((b'<!doctype html', b'<html', b'<?xml', b'{"error"')):
            raise ValueError('ENA returned an HTML/XML/JSON error body instead of TSV')
        raw.seek(0, os.SEEK_END)
        size = raw.tell()
        if size == 0:
            raise ValueError('empty ENA response')
        raw.seek(-1, os.SEEK_END)
        if raw.read(1) != b'\n':
            raise ValueError('ENA response has no final newline; possible truncation')
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        # ENA Portal TSV uses literal tabs; double quotes in metadata are data.
        # QUOTE_NONE keeps them unchanged. Embedded tabs/newlines are rejected
        # by the exact field-count and row-count checks below.
        reader = csv.reader(f, delimiter='\t', quotechar=None,
                            quoting=csv.QUOTE_NONE, escapechar=None, strict=True)
        if next(reader, None) != FIELDS:
            raise ValueError('header/schema mismatch')
        for row in reader:
            n += 1
            if len(row) != len(FIELDS):
                raise ValueError(f'column count mismatch at row {n}')
            x = dict(zip(FIELDS, row))
            acc = x['run_accession']
            if not acc or acc in seen:
                raise ValueError(f'blank/duplicate run_accession at row {n}: {acc}')
            seen.add(acc)
            if (x['instrument_platform'] != partition['platform'] or
                x['library_layout'] != LAYOUT[partition['platform']] or
                x['library_strategy'] != 'WGS' or x['library_source'] != 'GENOMIC'):
                raise ValueError(f'scope mismatch at row {n}: {acc}')
            try:
                public = date(x['first_public'][:10])
            except ValueError as e:
                raise ValueError(f'missing/unparseable first_public at row {n}: {acc}') from e
            if not date(partition['start']) <= public < min(date(partition['end']), date(cutoff) + dt.timedelta(days=1)):
                raise ValueError(f'date outside partition at row {n}: {acc}')
    return n, len(seen)


def verify_entry(entry, root, cutoff, expected=None):
    if entry.get('fields') != ','.join(FIELDS) or (expected is not None and
            (entry.get('query') != expected['query'] or entry.get('start') != expected['start'] or
             entry.get('end') != expected['end'] or entry.get('platform') != expected['platform'])):
        return False
    path = root / entry['file']
    if not path.is_file() or sha256(path) != entry['gzip_sha256'] or path.stat().st_size != int(entry['gzip_bytes']):
        return False
    temp = root / 'tmp' / (entry['id'] + '.resume.partial')
    try:
        with gzip.open(path, 'rb') as source, open(temp, 'wb') as target:
            shutil.copyfileobj(source, target)
        if sha256(temp) != entry['raw_sha256'] or temp.stat().st_size != int(entry['raw_bytes']):
            return False
        rows, unique = validate_tsv(temp, entry, cutoff)
        return rows == int(entry['rows']) == unique == int(entry['unique_runs']) == int(entry['expected_before']) == int(entry['expected_after'])
    except (OSError, ValueError, csv.Error, UnicodeError, EOFError):
        return False
    finally:
        temp.unlink(missing_ok=True)


def gzip_file(source, target):
    with open(source, 'rb') as src, open(target, 'wb') as dst:
        with gzip.GzipFile(filename='', mode='wb', compresslevel=6, mtime=0, fileobj=dst) as gz:
            shutil.copyfileobj(src, gz)
        dst.flush()
        os.fsync(dst.fileno())


def save_manifest(root, entries, plan):
    rows = [entries[p['id']] for p in plan if p['id'] in entries]
    atomic_text(root / 'manifest.tsv', tsv_string(rows, MANIFEST_COLUMNS))


def load_tsv(path):
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def log_attempt(root, record):
    path = root / 'logs' / 'attempts.jsonl'
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, sort_keys=True) + '\n')
        f.flush()
        os.fsync(f.fileno())
    print(json.dumps(record, sort_keys=True), flush=True)


def capture_one(ena, cfg, root, p):
    temp = root / 'tmp' / (p['id'] + '.tsv.partial')
    zipped = root / 'tmp' / (p['id'] + '.tsv.gz.partial')
    final = root / 'raw' / p['platform'] / (p['id'] + '.tsv.gz')
    final.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(cfg['max_attempts']):
        temp.unlink(missing_ok=True)
        zipped.unlink(missing_ok=True)
        began = utc()
        record = dict(id=p['id'], attempt=attempt + 1, utc_start=began, query=p['query'])
        try:
            before = ena.count(p['query'])
            record['expected_before'] = before
            if before > cfg['max_partition_rows']:
                record['status'] = 'TOO_LARGE'
                log_attempt(root, record)
                return None
            ena.download(p['query'], temp)
            rows, unique = validate_tsv(temp, p, cfg['cutoff_date'])
            after = ena.count(p['query'])
            record.update(expected_after=after, rows=rows, unique_runs=unique)
            if before != rows or rows != unique or unique != after:
                raise ValueError(f'count mismatch: {before}, {rows}, {unique}, {after}')
            raw_hash, raw_bytes = sha256(temp), temp.stat().st_size
            gzip_file(temp, zipped)
            gz_hash, gz_bytes = sha256(zipped), zipped.stat().st_size
            os.replace(zipped, final)
            entry = dict(p, fields=','.join(FIELDS), status='COMPLETE', expected_before=before, expected_after=after,
                rows=rows, unique_runs=unique, raw_bytes=raw_bytes, gzip_bytes=gz_bytes,
                raw_sha256=raw_hash, gzip_sha256=gz_hash, utc_start=began, utc_end=utc(),
                retry_count=attempt, file=str(final.relative_to(root)))
            record.update(status='COMPLETE', utc_end=entry['utc_end'])
            log_attempt(root, record)
            return entry
        except (OSError, ValueError, csv.Error, UnicodeError, TimeoutError, urllib.error.URLError) as e:
            last_error = e
            record.update(status='FAILED', error_type=f'{type(e).__module__}.{type(e).__name__}',
                          error=str(e), utc_end=utc())
            if temp.is_file():
                record['response_bytes'] = temp.stat().st_size
                record['response_sha256'] = sha256(temp)
                quarantine = root / 'logs' / 'failed_responses' / f'{p["id"]}.attempt{attempt + 1}.{time.time_ns()}.tsv'
                quarantine.parent.mkdir(parents=True, exist_ok=True)
                try:
                    os.replace(temp, quarantine)
                    record['quarantine_file'] = str(quarantine.relative_to(root))
                except OSError as quarantine_error:
                    record['quarantine_error'] = str(quarantine_error)
            log_attempt(root, record)
            if attempt + 1 < cfg['max_attempts']:
                time.sleep(min(cfg['retry_base_seconds'] * 2 ** attempt, 120))
        finally:
            temp.unlink(missing_ok=True)
            zipped.unlink(missing_ok=True)
    log_attempt(root, dict(status='FAILED_EXHAUSTED', id=p['id'],
                           attempts=cfg['max_attempts'], utc=utc(),
                           error_type=f'{type(last_error).__module__}.{type(last_error).__name__}',
                           error=str(last_error)))
    raise RuntimeError(f'partition {p["id"]} failed after {cfg["max_attempts"]} attempts; '
                       f'completed partitions are preserved; resume with the same config after checking logs/attempts.jsonl') from last_error


def split_partition(p):
    days = (date(p['end']) - date(p['start'])).days
    if days == 1:
        raise RuntimeError(f'one-day partition cannot be subdivided: {p["id"]}')
    unit = 'month' if days > 31 else 'day'
    pieces = [make_partition(p['platform'], a, b) for a, b in intervals(p['start'], p['end'], unit)]
    if len(pieces) == 1:
        unit = 'day'
        pieces = [make_partition(p['platform'], a, b) for a, b in intervals(p['start'], p['end'], unit)]
    return pieces


def preflight(ena, cfg, root):
    field_hash = ena.fields()
    audit = {'utc': utc(), 'return_fields_sha256': field_hash, 'platforms': {}}
    for p in cfg['platforms']:
        q = platform_query(p)
        before = ena.count(q + f' AND first_public<{cfg["start_date"]}')
        within = ena.count(range_query(p, cfg['start_date'], cutoff_exclusive(cfg)))
        after = ena.count(q + f' AND first_public>={cutoff_exclusive(cfg)}')
        broad = ena.count(q)
        unknown = broad - before - within - after
        item = dict(broad=broad, before_start=before, within=within,
                    after_cutoff=after, unexplained=unknown)
        if unknown:
            # Test whether the Portal query language can address the undated bucket.
            # It cannot establish a first-public cutoff for undated records, so stop.
            missing_query = q + f' AND NOT (first_public<{cfg["start_date"]} OR first_public>={cfg["start_date"]})'
            item['missing_date_query_test'] = missing_query
            try:
                item['missing_date_query_count'] = ena.count(missing_query)
            except Exception as e:
                item['missing_date_query_error'] = str(e)
        audit['platforms'][p] = item
    atomic_text(root / 'logs' / 'preflight.json', json.dumps(audit, indent=2) + '\n')
    if not cfg.get('smoke') and any(x['before_start'] or x['unexplained'] for x in audit['platforms'].values()):
        raise RuntimeError('preflight found records before plan start or missing/unparseable dates; see logs/preflight.json')
    return audit


def global_validate(ena, cfg, root, plan, entries):
    dbpath = root / 'tmp' / 'global_accessions.sqlite.partial'
    dbpath.unlink(missing_ok=True)
    db = sqlite3.connect(dbpath)
    db.execute('CREATE TABLE runs (accession TEXT PRIMARY KEY)')
    total = 0
    by_year, by_platform, by_both = {}, {}, {}
    combined = root / 'combined' / 'ena_supported_wgs_metadata.tsv.gz'
    combined.parent.mkdir(parents=True, exist_ok=True)
    tmp = combined.with_name(combined.name + '.partial')
    try:
        with open(tmp, 'wb') as raw:
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0, compresslevel=6) as out:
                out.write(('\t'.join(FIELDS) + '\n').encode())
                for p in plan:
                    with gzip.open(root / entries[p['id']]['file'], 'rt', encoding='utf-8-sig', newline='') as src:
                        reader = csv.DictReader(src, delimiter='\t', quotechar=None,
                                                quoting=csv.QUOTE_NONE, escapechar=None, strict=True)
                        if reader.fieldnames != FIELDS:
                            raise ValueError(f'combined schema mismatch: {p["id"]}')
                        for row in reader:
                            if None in row or any(v is None for v in row.values()):
                                raise ValueError(f'bad combined row: {p["id"]}')
                            a = row['run_accession']
                            try:
                                db.execute('INSERT INTO runs VALUES (?)', (a,))
                            except sqlite3.IntegrityError as e:
                                raise ValueError(f'global duplicate run: {a}') from e
                            year = row['first_public'][:4]
                            platform = row['instrument_platform']
                            by_year[year] = by_year.get(year, 0) + 1
                            by_platform[platform] = by_platform.get(platform, 0) + 1
                            key = (platform, year)
                            by_both[key] = by_both.get(key, 0) + 1
                            out.write(('\t'.join(row[f] for f in FIELDS) + '\n').encode())
                            total += 1
            raw.flush()
            os.fsync(raw.fileno())
        unique = db.execute('SELECT count(*) FROM runs').fetchone()[0]
        counts = {p: ena.count(range_query(p, cfg['start_date'], cutoff_exclusive(cfg))) for p in cfg['platforms']}
        if cfg.get('smoke'):
            union = sum(counts.values())
        else:
            union = ena.count(union_query() + f' AND first_public>={cfg["start_date"]} AND first_public<{cutoff_exclusive(cfg)}')
        audit = dict(total_rows=total, unique_runs=unique, platform_counts=counts,
                     union_server_count=union, local_platform_counts=by_platform,
                     utc=utc())
        atomic_text(root / 'logs' / 'global_validation.json', json.dumps(audit, indent=2) + '\n')
        if not (total == unique == union == sum(counts.values()) == sum(int(x['rows']) for x in entries.values())):
            raise ValueError('global count mismatch; see logs/global_validation.json')
        if any(by_platform.get(p, 0) != counts[p] for p in cfg['platforms']):
            raise ValueError('global platform mismatch; see logs/global_validation.json')
        os.replace(tmp, combined)
        census = [{'platform': p, 'year': y, 'rows': n} for (p, y), n in sorted(by_both.items())]
        atomic_text(root / 'census_platform_year.tsv', tsv_string(census, ['platform', 'year', 'rows']))
        atomic_text(root / 'census_year.tsv', tsv_string([{'year': y, 'rows': n} for y, n in sorted(by_year.items())], ['year', 'rows']))
        atomic_text(root / 'census_platform.tsv', tsv_string([{'platform': p, 'rows': n} for p, n in sorted(by_platform.items())], ['platform', 'rows']))
        return total, unique, combined, audit
    finally:
        db.close()
        dbpath.unlink(missing_ok=True)
        tmp.unlink(missing_ok=True)


def run(cfg, root, allow_download=True):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock = root / ('SMOKE_COMPLETE.lock' if cfg.get('smoke') else 'SNAPSHOT_COMPLETE.lock')
    if lock.exists():
        raise RuntimeError(f'completed snapshot is immutable: {lock}')
    (root / 'tmp').mkdir(exist_ok=True)
    for p in cfg['platforms']:
        (root / 'raw' / p).mkdir(parents=True, exist_ok=True)
    config_path = root / 'snapshot_config.json'
    config_text = json.dumps(cfg, indent=2, sort_keys=True) + '\n'
    if config_path.exists() and config_path.read_text() != config_text:
        raise RuntimeError('snapshot_config.json differs; refusing scientific query change')
    atomic_text(config_path, config_text)
    plan_path = root / 'partition_plan.tsv'
    plan = load_tsv(plan_path) if plan_path.exists() else initial_plan(cfg)
    validate_plan(plan, cfg)
    if not plan_path.exists():
        atomic_text(plan_path, tsv_string(plan, PLAN_COLUMNS))
    start_path = root / 'logs' / 'retrieval_started_utc.txt'
    if not start_path.exists():
        atomic_text(start_path, utc() + '\n')
    ena = ENA(cfg)
    ena.log_root = root
    log_attempt(root, dict(status='START_OR_RESUME', utc=utc(), snapshot_id=cfg['snapshot_id'], partitions=len(plan)))
    pre = preflight(ena, cfg, root)
    log_attempt(root, dict(status='PREFLIGHT_COMPLETE', utc=utc(), summary=pre['platforms']))
    manifest_path = root / 'manifest.tsv'
    entries = {r['id']: r for r in load_tsv(manifest_path)} if manifest_path.exists() else {}
    for p in plan:
        if p['id'] in entries and not verify_entry(entries[p['id']], root, cfg['cutoff_date'], p):
            log_attempt(root, dict(id=p['id'], status='LOCAL_INVALID', utc=utc()))
            del entries[p['id']]
    save_manifest(root, entries, plan)
    for round_no in range(cfg['final_audit_rounds']):
        i = 0
        while i < len(plan):
            p = plan[i]
            if p['id'] not in entries:
                if not allow_download:
                    raise RuntimeError(f'repaired partition missing or invalid: {p["id"]}; no snapshot lock created')
                result = capture_one(ena, cfg, root, p)
                if result is None:
                    pieces = split_partition(p)
                    plan[i:i + 1] = pieces
                    validate_plan(plan, cfg)
                    atomic_text(plan_path, tsv_string(plan, PLAN_COLUMNS))
                    continue
                entries[p['id']] = result
                save_manifest(root, entries, plan)
            i += 1
        drift = []
        for p in plan:
            current = ena.count(p['query'])
            log_attempt(root, dict(status='FINAL_AUDIT_COUNT', utc=utc(), id=p['id'], server=current))
            if current != int(entries[p['id']]['rows']):
                drift.append(dict(id=p['id'], local=entries[p['id']]['rows'], server=current))
        atomic_text(root / 'logs' / 'final_partition_audit.json', json.dumps({'utc': utc(), 'round': round_no + 1, 'drift': drift}, indent=2) + '\n')
        if not drift:
            break
        if not allow_download:
            raise RuntimeError(f'repaired partition drift: {drift}; no snapshot lock created')
        for item in drift:
            del entries[item['id']]
        save_manifest(root, entries, plan)
        log_attempt(root, dict(status='ARCHIVE_DRIFT', partitions=drift, utc=utc()))
    else:
        raise RuntimeError('partition drift persists; no snapshot lock created')
    log_attempt(root, dict(status='GLOBAL_VALIDATE_START', utc=utc(), partitions=len(plan)))
    total, unique, combined, audit = global_validate(ena, cfg, root, plan, entries)
    log_attempt(root, dict(status='GLOBAL_VALIDATE_COMPLETE', utc=utc(), rows=total, unique=unique))
    if not cfg.get('smoke'):
        # Recheck missing-date arithmetic and cutoff scope immediately before freeze.
        final_pre = preflight(ena, cfg, root)
        for p in cfg['platforms']:
            if final_pre['platforms'][p]['within'] != audit['platform_counts'][p]:
                raise RuntimeError(f'archive drift after global validation: {p}')
    manifest_hash = sha256(manifest_path)
    combined_hash = sha256(combined)
    files = [root / e['file'] for e in entries.values()] + [manifest_path, plan_path, config_path, combined,
        root / 'census_year.tsv', root / 'census_platform.tsv', root / 'census_platform_year.tsv']
    atomic_text(root / 'SHA256SUMS', ''.join(f'{sha256(f)}  {f.relative_to(root)}\n' for f in sorted(files)))
    provenance = dict(snapshot_id=cfg['snapshot_id'], retrieval_start_utc=start_path.read_text().strip(),
        retrieval_end_utc=utc(), cutoff_date=cfg['cutoff_date'], fields=FIELDS,
        fields_return_response_sha256=pre['return_fields_sha256'],
        plan_sha256=sha256(plan_path), manifest_sha256=manifest_hash,
        combined_sha256=combined_hash, global_validation=audit,
        missing_date_assessment=pre['platforms'])
    atomic_text(root / 'provenance.json', json.dumps(provenance, indent=2, sort_keys=True) + '\n')
    atomic_text(root / 'README.md', f'# {cfg["snapshot_id"]}\n\nENA read_run metadata capture. first_public is a timestamp; all intervals are [start, end_exclusive). The cutoff includes the whole day {cfg["cutoff_date"]}, through records before {cutoff_exclusive(cfg)}. See provenance.json, partition_plan.tsv, manifest.tsv and SHA256SUMS. Raw gzip partitions are the archival data; combined/ is a derivative.\n')
    atomic_text(lock, json.dumps(dict(snapshot_id=cfg['snapshot_id'], cutoff_date=cfg['cutoff_date'],
        utc_completion=utc(), total_rows=total, total_unique_runs=unique,
        manifest_sha256=manifest_hash, combined_sha256=combined_hash), indent=2) + '\n')
    print(f'COMPLETE {lock} rows={total} unique={unique}', flush=True)


def main():
    os.umask(0o077)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    required={'platforms','start_date','cutoff_date','api_base','max_attempts','connect_timeout_seconds','transfer_timeout_seconds','retry_base_seconds','max_partition_rows','final_audit_rounds','snapshot_id'}
    if required-set(cfg): raise ValueError('missing archive config keys: '+str(sorted(required-set(cfg))))
    if cfg['platforms'] != list(LAYOUT) and not cfg.get('smoke'):
        raise ValueError('full capture must use all four canonical platform/layout combinations')
    if cfg['start_date'] > cfg['cutoff_date']:
        raise ValueError('start after cutoff')
    run(cfg, args.output)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr, flush=True)
        sys.exit(1)
