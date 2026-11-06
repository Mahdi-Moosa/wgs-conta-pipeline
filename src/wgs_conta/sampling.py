"""Generic outcome-blind replicate cohort frame, allocation, ranking and selection.

Scientific eligibility is supplied by the frozen `primary_eligible` field.
Geography is constrained by the requested country rule and a caller-supplied
ISO3 -> frozen UN M49 mapping. No archive-specific stratum counts live here.
"""
from __future__ import annotations

import csv
import hashlib
import hmac
import json
import os
import secrets
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

CODE_VERSION='WGS_CONTA_REPLICATE_SELECTOR_V1'
DOMAIN_LABELS={
    'project':'WGS_CONTA_REPLICATE_V1|PROJECT',
    'cell':'WGS_CONTA_REPLICATE_V1|CELL',
    'biosample':'WGS_CONTA_REPLICATE_V1|BIOSAMPLE',
    'run':'WGS_CONTA_REPLICATE_V1|RUN',
    'pseudonym':'WGS_CONTA_REPLICATE_V1|PSEUDONYM',
}


@dataclass(frozen=True)
class RunRecord:
    accession:str
    country_iso3:str
    platform:str
    first_public:str=''


@dataclass
class Cell:
    project:str
    scientific_name:str
    instrument_model:str
    subregion:str
    biosamples:dict[str,list[RunRecord]]=field(default_factory=dict)

    @property
    def stable_id(self):
        return _encode((self.project,self.scientific_name,self.instrument_model,self.subregion))


@dataclass
class Project:
    accession:str
    subregion:str
    broad_region:str
    countries:tuple[str,...]
    cells:dict[tuple[str,str],Cell]=field(default_factory=dict)


@dataclass
class SamplingFrame:
    projects:dict[str,Project]
    summary:dict
    input_counts:dict

    def capacities(self):
        result={}
        for p in self.projects.values():
            row=result.setdefault(p.subregion,dict(m49_subregion=p.subregion,
                m49_broad_region=p.broad_region,qualifying_projects=0,
                qualifying_cells=0,eligible_biosamples_in_cells=0,eligible_runs_in_cells=0))
            row['qualifying_projects']+=1
            for cell in p.cells.values():
                row['qualifying_cells']+=1
                row['eligible_biosamples_in_cells']+=len(cell.biosamples)
                row['eligible_runs_in_cells']+=sum(len(runs) for runs in cell.biosamples.values())
        return dict(sorted(result.items()))


@dataclass(frozen=True)
class Ranked:
    entity:object
    digest:str
    rank:int


def _encode(parts):
    return json.dumps(parts,ensure_ascii=False,separators=(',',':'))


def _clean(value):
    if value is None:return ''
    return str(value).strip()


def _bool(value):
    if isinstance(value,bool):return value
    return str(value).strip().lower() in ('true','1','yes')


def load_m49_mapping(path):
    """Read only the frozen ISO3/region/subregion columns from the supplied TSV."""
    mapping={}
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'iso3','region_un','subregion_un'}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f'M49 map requires columns {sorted(required)}')
        for row in reader:
            iso=_clean(row.get('iso3')).upper()
            sub=_clean(row.get('subregion_un'))
            region=_clean(row.get('region_un'))
            if not iso or not sub:continue
            if len(iso)!=3 or not iso.isalpha():raise ValueError(f'invalid M49 ISO3 key: {iso!r}')
            value=(sub,region)
            if iso in mapping and mapping[iso]!=value:
                raise ValueError(f'conflicting M49 rows for {iso}: {mapping[iso]} vs {value}')
            mapping[iso]=value
    if not mapping:raise ValueError('M49 map contains no usable ISO3/subregion rows')
    return mapping


def build_sampling_frame(records:Iterable[Mapping],m49_mapping:Mapping[str,tuple[str,str]],
                         biosamples_per_project=3,country_rule='MATCH',input_counts=None):
    """Build the project frame from canonical metadata without reclassifying eligibility.

    Every eligible MATCH BioSample must have a mapped country and resolve to
    one M49 subregion. Projects with more than one subregion are excluded;
    distinct countries within one subregion remain allowed. Cells use the
    frozen audit key: BioProject × scientific_name × instrument_model ×
    project M49 subregion, with the requested number of distinct BioSamples.
    """
    n=int(biosamples_per_project)
    if n<1:raise ValueError('biosamples_per_project must be positive')
    if country_rule!='MATCH':raise ValueError('only country_rule=MATCH is supported')
    project_geo={}
    candidate_cells={}
    seen_runs=set()
    stats=dict(input_records=0,primary_eligible_records=0,country_rule_records=0,
        eligible_match_records=0,rows_missing_project_or_biosample=0,
        projects_unmapped=0,projects_multi_subregion=0,projects_no_subregion=0,
        projects_without_qualifying_cell=0,duplicate_eligible_run_accessions=0)
    for raw in records:
        stats['input_records']+=1
        if not _bool(raw.get('primary_eligible')):continue
        stats['primary_eligible_records']+=1
        if _clean(raw.get('country_compare_status')).upper()!=country_rule:continue
        stats['country_rule_records']+=1
        stats['eligible_match_records']+=1
        project=_clean(raw.get('canonical_bioproject'))
        sample=_clean(raw.get('canonical_biosample'))
        accession=_clean(raw.get('run_accession'))
        if not project or not sample or not accession:
            stats['rows_missing_project_or_biosample']+=1
            continue
        if accession in seen_runs:
            stats['duplicate_eligible_run_accessions']+=1
            raise ValueError(f'duplicate eligible run accession: {accession}')
        seen_runs.add(accession)
        country=_clean(raw.get('sample_country_iso3') or raw.get('ena_country_iso3')).upper()
        mapped=m49_mapping.get(country)
        subregion,broad=(mapped if mapped else ('',''))
        geo=project_geo.setdefault(project,dict(samples={},countries=set()))
        sample_geo=geo['samples'].setdefault(sample,dict(subregions=set(),unmapped=False,countries=set()))
        if subregion:sample_geo['subregions'].add(subregion)
        else:sample_geo['unmapped']=True
        if country:sample_geo['countries'].add(country)
        geo['countries'].add(country) if country else None
        species=_clean(raw.get('scientific_name'))
        instrument=_clean(raw.get('instrument_model'))
        if not species or not instrument:continue
        project_cells=candidate_cells.setdefault(project,{})
        sample_runs=project_cells.setdefault((species,instrument),{}).setdefault(sample,[])
        sample_runs.append(RunRecord(accession,country,_clean(raw.get('instrument_platform')),
                                     _clean(raw.get('first_public'))))

    projects={}
    exclusion_status={}
    for project,geo in project_geo.items():
        sample_regions=[]
        unmapped=any(s['unmapped'] or len(s['subregions'])!=1 for s in geo['samples'].values())
        if unmapped:
            stats['projects_unmapped']+=1;exclusion_status[project]='HAS_UNMAPPED_OR_UNRESOLVED_M49_SAMPLE';continue
        for sample_geo in geo['samples'].values():sample_regions.extend(sample_geo['subregions'])
        subregions=set(sample_regions)
        if not subregions:
            stats['projects_no_subregion']+=1;exclusion_status[project]='NO_M49_SUBREGION';continue
        if len(subregions)>1:
            stats['projects_multi_subregion']+=1;exclusion_status[project]='MULTI_M49_SUBREGION';continue
        subregion=next(iter(subregions))
        # Recover the frozen broad-region label from any mapped sample country.
        broad=''
        for iso in sorted(geo['countries']):
            value=m49_mapping.get(iso)
            if value and value[0]==subregion:
                broad=value[1];break
        cells={}
        for (species,instrument),sample_runs in candidate_cells.get(project,{}).items():
            biosamples={sample:sorted(runs,key=lambda r:(r.accession,r.country_iso3,r.platform,r.first_public))
                        for sample,runs in sample_runs.items()}
            if len(biosamples)<n:continue
            key=(species,instrument)
            cells[key]=Cell(project,species,instrument,subregion,biosamples)
        if not cells:
            stats['projects_without_qualifying_cell']+=1;exclusion_status[project]='NO_QUALIFYING_CELL';continue
        countries=tuple(sorted(x for x in geo['countries'] if x))
        projects[project]=Project(project,subregion,broad,countries,cells)
    capacities={}
    for p in projects.values():capacities[p.subregion]=capacities.get(p.subregion,0)+1
    all_cells=sum(len(p.cells) for p in projects.values())
    frame_biosamples=sum(len(c.biosamples) for p in projects.values() for c in p.cells.values())
    frame_runs=sum(sum(len(rs) for rs in c.biosamples.values()) for p in projects.values() for c in p.cells.values())
    summary={**stats,'qualifying_projects':len(projects),'qualifying_cells':all_cells,
        'm49_subregions':len(capacities),'eligible_biosamples_in_qualifying_cells':frame_biosamples,
        'eligible_runs_in_qualifying_cells':frame_runs,'projects_by_geography_status':{
            'IN_SAMPLING_FRAME':len(projects),**_count_values(exclusion_status.values())},
        'projects_excluded_total':len(exclusion_status)}
    if input_counts:summary.update(input_counts)
    return SamplingFrame(projects,summary,input_counts or {})


def _count_values(values):
    out={}
    for x in values:out[x]=out.get(x,0)+1
    return out


def allocate_hamilton(stratum_capacity,projects):
    """Proportional allocation with deterministic largest-remainder ties.

    Tie order matches the frozen M49 audit: larger remainder, larger stratum,
    then lexical M49 subregion. No minimum-one rule is imposed.
    """
    n=int(projects)
    if n<0:raise ValueError('projects cannot be negative')
    raw={str(k):int(v['qualifying_projects'] if isinstance(v,Mapping) else v)
         for k,v in stratum_capacity.items()}
    if any(v<0 for v in raw.values()):raise ValueError('stratum capacities cannot be negative')
    capacities={k:v for k,v in raw.items() if v>0}
    total=sum(capacities.values())
    if n>total:raise ValueError(f'requested {n} projects exceeds frame capacity {total}')
    if n and not total:raise ValueError('cannot allocate from an empty frame')
    rows=[]
    for subregion,nh in sorted(capacities.items()):
        ideal=n*nh/total if total else 0.0
        floor_n=int(ideal//1)
        rows.append(dict(m49_subregion=subregion,eligible_projects=nh,frame_fraction=(nh/total if total else 0),
            ideal_allocation=ideal,floor_allocation=floor_n,fractional_remainder=ideal-floor_n,
            final_allocation=floor_n,realized_selected_project_count=0))
    left=n-sum(x['floor_allocation'] for x in rows)
    ranked=sorted(rows,key=lambda r:(-r['fractional_remainder'],-r['eligible_projects'],r['m49_subregion']))
    for row in ranked[:left]:row['final_allocation']+=1
    if sum(x['final_allocation'] for x in rows)!=n:raise AssertionError('Hamilton allocation failed to sum to target')
    return rows


def _digest(seed,domain,cohort_label,*parts):
    if len(seed)!=32:raise ValueError('selection and pseudonym keys must be exactly 32 bytes')
    message=_encode([CODE_VERSION,DOMAIN_LABELS[domain],str(cohort_label),*parts]).encode('utf-8')
    return hmac.new(seed,message,hashlib.sha256).hexdigest()


def _sort_ranked(entities,seed,domain,label,keyfn):
    items=[(entity,_digest(seed,domain,label,*keyfn(entity)),keyfn(entity)) for entity in entities]
    items.sort(key=lambda x:(x[1],_encode(x[2])))
    return [Ranked(entity,digest,i+1) for i,(entity,digest,_) in enumerate(items)]


def rank_projects(frame,seed,cohort_label):
    ranked=[]
    for subregion in sorted(frame.capacities()):
        projects=[p for p in frame.projects.values() if p.subregion==subregion]
        ranked.extend((subregion,rank) for rank in _sort_ranked(
            projects,seed,'project',cohort_label,lambda p:(p.subregion,p.accession)))
    return ranked


def rank_cells(project,seed,cohort_label):
    return _sort_ranked(project.cells.values(),seed,'cell',cohort_label,
        lambda c:(c.project,c.subregion,c.scientific_name,c.instrument_model))


def rank_biosamples(project,cell,seed,cohort_label):
    return _sort_ranked(sorted(cell.biosamples),seed,'biosample',cohort_label,
        lambda sample:(project.accession,cell.subregion,cell.scientific_name,cell.instrument_model,sample))


def rank_runs(project,cell,biosample,runs,seed,cohort_label):
    return _sort_ranked(runs,seed,'run',cohort_label,
        lambda run:(project.accession,cell.subregion,cell.scientific_name,cell.instrument_model,biosample,run.accession))


def select_cohort(frame,allocation,biosamples_per_project=3,runs_per_biosample=1,
                  seed=b'',cohort_label='B',keep_full_rankings=False):
    bpp,rpb=int(biosamples_per_project),int(runs_per_biosample)
    if bpp<1 or rpb<1:raise ValueError('biosamples_per_project and runs_per_biosample must be positive')
    if len(seed)!=32:raise ValueError('selection seed must contain exactly 32 random bytes')
    capacities={r['m49_subregion']:r['final_allocation'] for r in allocation}
    project_ranks=rank_projects(frame,seed,cohort_label)
    project_rank_rows=[]
    selected_projects=[]
    for subregion,ranked_group in _group_ranked_projects(project_ranks).items():
        quota=capacities.get(subregion,0)
        for item in ranked_group:
            project_rank_rows.append(dict(m49_subregion=subregion,project=item.entity.accession,
                project_rank=item.rank,rank_hmac_sha256=item.digest,selected=item.rank<=quota))
            if item.rank<=quota:selected_projects.append((item.entity,item.rank,item.digest))
    if len(selected_projects)!=sum(capacities.values()):raise AssertionError('project draw did not meet Hamilton allocation')
    membership=[];cell_rows=[];sample_rows=[];run_rows=[]
    selected_cell_ids={}; selected_sample_ids=set(); selected_run_ids=set()
    for project,project_rank,project_digest in selected_projects:
        cells_ranked=rank_cells(project,seed,cohort_label)
        selected_cell=None;chosen_samples=None
        cell_rank_by_id={}
        for cr in cells_ranked:
            cell=cr.entity;cell_rank_by_id[cell.stable_id]=cr.rank
            samples_ranked=rank_biosamples(project,cell,seed,cohort_label)
            viable=[r for r in samples_ranked if len(cell.biosamples[r.entity])>=rpb]
            if selected_cell is None and len(viable)>=bpp:
                selected_cell=cell;chosen_samples=viable[:bpp]
        if selected_cell is None:
            raise ValueError(f'project {project.accession} has no cell with {bpp} BioSamples and {rpb} runs each')
        selected_cell_ids[project.accession]=selected_cell.stable_id
        selected_sample_ids.update((project.accession,selected_cell.stable_id,x.entity) for x in chosen_samples)
        # Project ranks are complete. Cell/sample/run rankings are complete for selected projects.
        for cr in cells_ranked:
            cell=cr.entity;chosen=cell.stable_id==selected_cell.stable_id
            cell_rows.append(dict(project=project.accession,scientific_name=cell.scientific_name,
                instrument_model=cell.instrument_model,m49_subregion=cell.subregion,
                qualifying_biosamples=len(cell.biosamples),cell_rank=cr.rank,
                rank_hmac_sha256=cr.digest,selected=chosen))
            for sr in rank_biosamples(project,cell,seed,cohort_label):
                bio=sr.entity;bio_selected=chosen and any(x.entity==bio for x in chosen_samples)
                sample_rows.append(dict(project=project.accession,scientific_name=cell.scientific_name,
                    instrument_model=cell.instrument_model,m49_subregion=cell.subregion,biosample=bio,
                    qualifying_runs=len(cell.biosamples[bio]),biosample_rank=sr.rank,
                    rank_hmac_sha256=sr.digest,selected=bio_selected))
                for rr in rank_runs(project,cell,bio,cell.biosamples[bio],seed,cohort_label):
                    is_selected=bio_selected and rr.rank<=rpb
                    run_rows.append(dict(project=project.accession,scientific_name=cell.scientific_name,
                        instrument_model=cell.instrument_model,m49_subregion=cell.subregion,biosample=bio,
                        run_accession=rr.entity.accession,country_iso3=rr.entity.country_iso3,
                        platform=rr.entity.platform,first_public=rr.entity.first_public,
                        run_rank=rr.rank,rank_hmac_sha256=rr.digest,selected=is_selected))
                    if is_selected:
                        selected_run_ids.add(rr.entity.accession)
                        membership.append(dict(project=project.accession,biosample=bio,
                            run_accession=rr.entity.accession,scientific_name=cell.scientific_name,
                            instrument_model=cell.instrument_model,platform=rr.entity.platform,
                            country_iso3=rr.entity.country_iso3,m49_subregion=cell.subregion,
                            project_rank=project_rank,cell_rank=cr.rank,biosample_rank=sr.rank,run_rank=rr.rank,
                            project_rank_hmac=project_digest,cell_rank_hmac=cr.digest,
                            biosample_rank_hmac=sr.digest,run_rank_hmac=rr.digest))
    if keep_full_rankings:
        selected_ids={p.accession for p,_,_ in selected_projects}
        # Optional reserve artifact: full within-project order for every
        # unselected project, with selection flags false throughout.
        for project in sorted((p for p in frame.projects.values() if p.accession not in selected_ids),
                              key=lambda p:(p.subregion,p.accession)):
            for cr in rank_cells(project,seed,cohort_label):
                cell=cr.entity
                cell_rows.append(dict(project=project.accession,scientific_name=cell.scientific_name,
                    instrument_model=cell.instrument_model,m49_subregion=cell.subregion,
                    qualifying_biosamples=len(cell.biosamples),cell_rank=cr.rank,
                    rank_hmac_sha256=cr.digest,selected=False))
                for sr in rank_biosamples(project,cell,seed,cohort_label):
                    sample=sr.entity
                    sample_rows.append(dict(project=project.accession,scientific_name=cell.scientific_name,
                        instrument_model=cell.instrument_model,m49_subregion=cell.subregion,biosample=sample,
                        qualifying_runs=len(cell.biosamples[sample]),biosample_rank=sr.rank,
                        rank_hmac_sha256=sr.digest,selected=False))
                    for rr in rank_runs(project,cell,sample,cell.biosamples[sample],seed,cohort_label):
                        run_rows.append(dict(project=project.accession,scientific_name=cell.scientific_name,
                            instrument_model=cell.instrument_model,m49_subregion=cell.subregion,biosample=sample,
                            run_accession=rr.entity.accession,country_iso3=rr.entity.country_iso3,
                            platform=rr.entity.platform,first_public=rr.entity.first_public,
                            run_rank=rr.rank,rank_hmac_sha256=rr.digest,selected=False))
    if len(selected_run_ids)!=len(membership):raise AssertionError('selected run accession duplicated')
    return dict(membership=membership,project_ranking=project_rank_rows,
        cell_ranking=cell_rows,biosample_ranking=sample_rows,run_ranking=run_rows,
        selected_projects=[p.accession for p,_,_ in selected_projects],
        selected_cells=selected_cell_ids)


def _group_ranked_projects(rows):
    groups={}
    for sub,rank in rows:groups.setdefault(sub,[]).append(rank)
    return groups


def pseudonymize_membership(membership,key,cohort_label):
    """Assign independent HMAC-ordered public pseudonyms without accession sorting."""
    if len(key)!=32:raise ValueError('pseudonym key must contain exactly 32 random bytes')
    projects=sorted({r['project'] for r in membership})
    runs=sorted({r['run_accession'] for r in membership})
    project_order=sorted(projects,key=lambda x:(_digest(key,'pseudonym',cohort_label,'PROJECT',x),x))
    run_order=sorted(runs,key=lambda x:(_digest(key,'pseudonym',cohort_label,'RUN',x),x))
    width_project=max(3,len(str(len(project_order))))
    width_run=max(4,len(str(len(run_order))))
    pmap={p:f'BP{i:0{width_project}d}' for i,p in enumerate(project_order,1)}
    rmap={r:f'{cohort_label}{i:0{width_run}d}' for i,r in enumerate(run_order,1)}
    output=[]
    for row in membership:
        output.append({**row,'pseudonym':rmap[row['run_accession']],
                       'project_pseudonym':pmap[row['project']]})
    output.sort(key=lambda r:r['pseudonym'])
    return output


def compare_memberships(selected,prior_path):
    prior={k:set() for k in ('project','biosample','run_accession')}
    with Path(prior_path).open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        aliases={'project':('project','bioproject','BioProject'),
                 'biosample':('biosample','BioSample'),
                 'run_accession':('run_accession','accession','run','Run')}
        names=set(reader.fieldnames or ())
        cols={k:next((x for x in options if x in names),None) for k,options in aliases.items()}
        for row in reader:
            for key,col in cols.items():
                if col and _clean(row.get(col)):prior[key].add(_clean(row[col]))
    current={k:{_clean(r.get(k)) for r in selected if _clean(r.get(k))} for k in prior}
    return {f'shared_{k}s':len(prior[k]&current[k]) for k in prior}


def summarize_selection(selection):
    """Metadata-only realized platform and instrument breadth for selected projects."""
    by_project={}
    for row in selection['membership']:
        item=by_project.setdefault(row['project'],{'platforms':set(),'models':set()})
        item['platforms'].add(row['platform'])
        item['models'].add(row['instrument_model'])
    platform_counts={}
    platform_projects={}
    for project,item in by_project.items():
        for platform in item['platforms']:
            platform_counts[platform]=platform_counts.get(platform,0)+1
            platform_projects.setdefault(platform,set()).add(project)
    platform_counts['BGI_MGI']=len(platform_projects.get('BGISEQ',set())|platform_projects.get('DNBSEQ',set()))
    return {'selected_projects':len(by_project),'projects_by_platform':dict(sorted(platform_counts.items())),
        'instrument_model_breadth':len({m for item in by_project.values() for m in item['models']}),
        'instrument_models':sorted({m for item in by_project.values() for m in item['models']})}


def generate_seed_file(path):
    """Explicit helper for a future authorized invocation; not used by plan-only."""
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    flags=os.O_WRONLY|os.O_CREAT|os.O_EXCL
    fd=os.open(path,flags,0o600)
    try:os.write(fd,secrets.token_bytes(32))
    finally:os.close(fd)
    return path


def read_key_file(path):
    key=Path(path).read_bytes()
    if len(key)!=32:raise ValueError(f'key file must contain exactly 32 random bytes: {path}')
    return key


def seed_sha256(key):
    if len(key)!=32:raise ValueError('key must contain exactly 32 bytes')
    return hashlib.sha256(key).hexdigest()


def verify_sha256_manifest(directory):
    directory=Path(directory)
    checked=[]
    for line in (directory/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split(None,1);name=name.strip()
        p=directory/name
        actual=hashlib.sha256(p.read_bytes()).hexdigest()
        if actual!=digest:raise ValueError(f'hash mismatch for {name}')
        checked.append(name)
    return checked
