"""Nonidentifying aggregate exports; no automatic release of per-run profiles."""
from collections import Counter

def aggregate(rows):
    if not rows:raise ValueError('empty result table')
    required={'Accession','Human_positive','Assembly_check','Quantification_unit','Unique_strict_units'}
    if not required.issubset(rows[0]):raise ValueError('missing result schema fields')
    if len({r['Accession'] for r in rows})!=len(rows):raise ValueError('select exactly one result version per run before aggregation')
    allowed={'PASS','BACTERIAL_MATCH_PRESENT','NO_EXACT_ASSEMBLY','FAILED','NOT_NEEDED_NEGATIVE'}
    counts=Counter();units=Counter();strict=completed=without=0
    for r in rows:
        flag=r['Human_positive'];status=r['Assembly_check'];unit=r['Quantification_unit']
        if flag not in {'YES','NO'} or status not in allowed or unit not in {'PAIRS','READS'}:raise ValueError('unknown endpoint status or quantification unit')
        n=int(r['Unique_strict_units'])
        if n<0 or (n>0)!=(flag=='YES'):raise ValueError('strict positivity/count mismatch')
        positive=flag=='YES'
        if (status=='NOT_NEEDED_NEGATIVE')!= (not positive):raise ValueError('strict/assembly state mismatch')
        counts[status]+=1;units[unit]+=n;strict+=positive
        completed+=positive and status in {'PASS','BACTERIAL_MATCH_PRESENT'}
        without+=positive and status=='PASS'
    return dict(runs=len(rows),strict_positive_runs=strict,assembly_status_counts=dict(sorted(counts.items())),
                positive_runs_with_completed_challenge=completed,positive_runs_without_detected_bacterial_match=without,
                strict_units_by_counting_unit=dict(sorted(units.items())),
                note='Missing/failed exact-assembly challenges remain separate; not classified as confirmed negatives. Do not pool PAIRS and READS.')
