from __future__ import annotations

import csv
import fcntl
import os
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

CHROMS=[f'chr{i}' for i in range(1,23)]+['chrX','chrY','mtDNA']
TRUE_WORDS={'1','true','yes','y','on','active'}
FALSE_WORDS={'0','false','no','n','off','inactive'}


def _norm(v):
    return '' if v is None else str(v).strip()


def _read_tsv(path:Path):
    path=Path(path)
    if not path.exists() or path.stat().st_size==0:
        return [],[]
    with path.open(newline='',encoding='utf-8') as fh:
        r=csv.DictReader(fh,delimiter='\t')
        if not r.fieldnames: return [],[]
        return list(r.fieldnames),list(r)


def _write_tsv_atomic(path:Path,fields,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+f'.tmp.{os.getpid()}')
    with tmp.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader()
        for row in rows: w.writerow({f:row.get(f,'') for f in fields})
    tmp.replace(path)


def _canonical_chromosome(name):
    x=_norm(name)
    if x in CHROMS: return x
    if x in {str(i) for i in range(1,23)}: return f'chr{x}'
    if x=='X': return 'chrX'
    if x=='Y': return 'chrY'
    if x in {'chrM','chrMT','MT','M','mitochondrion_genome','NC_012920.1','J01415.2'}: return 'mtDNA'
    raise RuntimeError(f'Unknown human reference sequence name {x!r}')


def _as_int(v):
    x=_norm(v)
    if x in {'','NA','nan','None'}: return None
    return int(float(x))


def _scientific_method(row):
    return _norm(row.get('Scientific_method_version')) or _norm(row.get('Pipeline_version'))


def parse_registry(path:Path):
    """Read a flexible master accession registry.

    Required: Accession.
    Optional: Label, Cohort, Active, Priority, Preferred_source, Notes.
    Unknown columns are retained in ``extra`` for future compatibility.
    """
    fields,rows=_read_tsv(path)
    if not fields: raise RuntimeError(f'Master accession registry is empty: {path}')
    if 'Accession' not in fields: raise RuntimeError(f'{path}: required column Accession is missing')
    out=[]; seen=set()
    for i,r in enumerate(rows,2):
        acc=_norm(r.get('Accession'))
        if not acc or acc.startswith('#'): continue
        if acc in seen: raise RuntimeError(f'{path}: duplicate accession {acc} at line {i}')
        seen.add(acc)
        active=_norm(r.get('Active')).lower()
        if active in FALSE_WORDS: continue
        if active and active not in TRUE_WORDS: raise RuntimeError(f'{path}: invalid Active value for {acc}: {r.get("Active")!r}')
        try: priority=int(_norm(r.get('Priority')) or '100')
        except ValueError as exc: raise RuntimeError(f'{path}: invalid Priority for {acc}: {r.get("Priority")!r}') from exc
        pref=_norm(r.get('Preferred_source')).lower()
        if pref=='ebi': pref='ena'
        if pref and pref not in {'aws','ncbi','ena','auto'}: raise RuntimeError(f'{path}: invalid Preferred_source for {acc}: {pref!r}')
        if pref=='auto': pref=''
        out.append({
            'accession':acc,
            'label':_norm(r.get('Label')) or _norm(r.get('Cohort')),
            'cohort':_norm(r.get('Cohort')),
            'priority':priority,
            'preferred_source':pref,
            'notes':_norm(r.get('Notes')),
            'extra':dict(r),
        })
    out.sort(key=lambda x:(x['priority'],x['accession']))
    return out


class MasterStore:
    def __init__(self,root:Path,pipeline_version:str,scientific_method_version:str):
        self.root=Path(root)
        self.results=self.root/'master_results.tsv'
        self.chromosomes=self.root/'chromosome_counts.tsv'
        self.cohorts=self.root/'cohorts.tsv'
        self.commit_log=self.root/'commit_log.tsv'
        self.lock_path=self.root/'.master.lock'
        self.pipeline_version=str(pipeline_version)
        self.scientific_method_version=str(scientific_method_version)
        self.root.mkdir(parents=True,exist_ok=True)

    @contextmanager
    def locked(self):
        self.lock_path.parent.mkdir(parents=True,exist_ok=True)
        with self.lock_path.open('a+') as fh:
            fcntl.flock(fh,fcntl.LOCK_EX)
            try: yield
            finally: fcntl.flock(fh,fcntl.LOCK_UN)

    def _matching_results(self,accession):
        _,rows=_read_tsv(self.results)
        return [r for r in rows if _norm(r.get('Accession'))==accession and _scientific_method(r)==self.scientific_method_version]

    def completion_detail(self,accession):
        """Return (complete, reason, result_row).

        A scientifically complete record requires one compatible result row and
        all 25 chromosome rows matching that result's Pipeline_version. This
        makes partial multi-file commits repairable on restart.
        """
        matches=self._matching_results(accession)
        if not matches: return False,'NO_COMPATIBLE_RESULT',None
        # Prefer newest row in file; duplicate compatible rows are tolerated but not created by v2.4.
        result=matches[-1]
        pver=_norm(result.get('Pipeline_version'))
        _,cr=_read_tsv(self.chromosomes)
        rows=[r for r in cr if _norm(r.get('Accession'))==accession and _norm(r.get('Pipeline_version'))==pver]
        by={_canonical_chromosome(r.get('Chromosome')):r for r in rows if _norm(r.get('Chromosome'))}
        if set(by)!=set(CHROMS): return False,f'CHROMOSOME_ROWS_{len(by)}_OF_25',result
        branch=_norm(result.get('Analysis_branch'))
        expected=_as_int(result.get('Unique_strict_fragments')) if branch=='PAIRED_SHORTREAD' else _as_int(result.get('Unique_strict_reads'))
        if expected is None: expected=_as_int(result.get('Unique_strict_units'))
        observed=sum(_as_int(by[c].get('Strict_units')) or 0 for c in CHROMS)
        if expected is None or observed!=expected: return False,f'CHROMOSOME_SUM_{observed}_EXPECTED_{expected}',result
        return True,'COMPLETE',result

    def is_complete(self,accession):
        return self.completion_detail(accession)[0]

    def _chromosome_rows_from_evidence(self,result,evidence_dir:Path):
        acc=_norm(result.get('Accession')); branch=_norm(result.get('Analysis_branch'))
        expected=_as_int(result.get('Unique_strict_fragments')) if branch=='PAIRED_SHORTREAD' else _as_int(result.get('Unique_strict_reads'))
        if expected is None: expected=_as_int(result.get('Unique_strict_units'))
        if expected is None: raise RuntimeError(f'{acc}: cannot determine unique strict count for chromosome commit')
        path=Path(evidence_dir)/'chromosome_counts.tsv'
        counts={c:0 for c in CHROMS}
        if path.exists():
            fields,rows=_read_tsv(path)
            if 'Chromosome' not in fields: raise RuntimeError(f'{path}: Chromosome column missing')
            count_field='Strict_fragments' if 'Strict_fragments' in fields else ('Strict_reads' if 'Strict_reads' in fields else None)
            if count_field is None: raise RuntimeError(f'{path}: strict-count column missing')
            for r in rows:
                c=_canonical_chromosome(r.get('Chromosome')); n=_as_int(r.get(count_field))
                if n is None: raise RuntimeError(f'{path}: missing count for {c}')
                counts[c]+=n
        elif expected:
            raise RuntimeError(f'{acc}: positive result but chromosome evidence is missing: {path}')
        if sum(counts.values())!=expected:
            raise RuntimeError(f'{acc}: chromosome evidence sum={sum(counts.values())}, expected unique strict={expected}')
        out=[]
        for c in CHROMS:
            n=counts[c]; frac=(n/expected if expected else 0.0)
            out.append({
                'Accession':acc,
                'Pipeline_version':self.pipeline_version,
                'Scientific_method_version':self.scientific_method_version,
                'Analysis_branch':branch,
                'Chromosome':c,
                'Strict_units':str(n),
                'Fraction_of_strict':f'{frac:.6f}',
                'Evidence_file':str(path),
            })
        return out

    def _merge_result(self,result,existing_fields,existing_rows):
        acc=result['Accession']
        # Idempotence for this exact pipeline version.
        same=[r for r in existing_rows if _norm(r.get('Accession'))==acc and _norm(r.get('Pipeline_version'))==self.pipeline_version]
        if same:
            old=same[-1]
            critical=['Analysis_branch','Unique_strict_units','Total_reads','Total_pairs','Total_bases','Human_positive','Assembly_check']
            diffs=[(k,_norm(old.get(k)),_norm(result.get(k))) for k in critical if _norm(old.get(k))!=_norm(result.get(k))]
            if diffs: raise RuntimeError(f'{acc}: conflicting existing v{self.pipeline_version} master result: {diffs}')
            return existing_fields,existing_rows,False
        fields=list(existing_fields)
        for k in result:
            if k not in fields: fields.append(k)
        return fields,existing_rows+[result],True

    def _merge_chromosomes(self,new_rows,existing_fields,existing_rows):
        fields=list(existing_fields)
        for r in new_rows:
            for k in r:
                if k not in fields: fields.append(k)
        idx={(_norm(r.get('Accession')),_norm(r.get('Pipeline_version')),_norm(r.get('Chromosome'))):r for r in existing_rows}
        added=0
        for r in new_rows:
            key=(_norm(r.get('Accession')),_norm(r.get('Pipeline_version')),_norm(r.get('Chromosome')))
            if key in idx:
                old=idx[key]
                if _norm(old.get('Strict_units'))!=_norm(r.get('Strict_units')):
                    raise RuntimeError(f'Chromosome conflict for {key}: old={old.get("Strict_units")} new={r.get("Strict_units")}')
                continue
            existing_rows.append(r); idx[key]=r; added+=1
        return fields,existing_rows,added

    def _merge_cohorts(self,accession,cohort,notes,fields,rows):
        if not fields: fields=['Accession','Cohort','Notes']
        cohorts=[x.strip() for x in cohort.replace(',', ';').split(';') if x.strip()]
        existing={(_norm(r.get('Accession')),_norm(r.get('Cohort'))) for r in rows}
        for c in cohorts:
            if (accession,c) not in existing:
                rows.append({'Accession':accession,'Cohort':c,'Notes':notes}); existing.add((accession,c))
        return fields,rows

    def commit(self,result,registry_record,evidence_dir:Path):
        """Validate and idempotently commit result + 25 chromosome rows + cohort.

        Multi-file atomicity is approximated with one cross-process file lock and
        idempotent merge semantics. ``is_complete`` requires both result and all
        chromosome rows, so a crash between os.replace operations is repaired on
        the next run rather than being treated as complete.
        """
        result=dict(result)
        result['Pipeline_version']=self.pipeline_version
        result['Scientific_method_version']=self.scientific_method_version
        result['Master_committed_at']=datetime.now().isoformat(timespec='seconds')
        chrom_new=self._chromosome_rows_from_evidence(result,evidence_dir)
        acc=_norm(result.get('Accession'))
        with self.locked():
            rf,rr=_read_tsv(self.results); cf,cr=_read_tsv(self.chromosomes); kf,kr=_read_tsv(self.cohorts)
            rf,rr,result_added=self._merge_result(result,rf,rr)
            cf,cr,chrom_added=self._merge_chromosomes(chrom_new,cf,cr)
            kf,kr=self._merge_cohorts(acc,registry_record.get('cohort',''),registry_record.get('notes',''),kf,kr)
            # Chromosomes first, result second: a crash after chromosomes but before
            # result remains safely incomplete and can be repaired idempotently.
            _write_tsv_atomic(self.chromosomes,cf,cr)
            _write_tsv_atomic(self.results,rf,rr)
            _write_tsv_atomic(self.cohorts,kf,kr)
            complete,reason,_=self.completion_detail(acc)
            if not complete: raise RuntimeError(f'{acc}: post-commit master verification failed: {reason}')
            lf,lr=_read_tsv(self.commit_log)
            if not lf: lf=['Committed_at','Accession','Pipeline_version','Scientific_method_version','Result_row_added','Chromosome_rows_added','Status']
            lr.append({'Committed_at':datetime.now().isoformat(timespec='seconds'),'Accession':acc,'Pipeline_version':self.pipeline_version,'Scientific_method_version':self.scientific_method_version,'Result_row_added':int(result_added),'Chromosome_rows_added':chrom_added,'Status':'COMPLETE'})
            _write_tsv_atomic(self.commit_log,lf,lr)
        return {'result_added':result_added,'chromosome_rows_added':chrom_added}
