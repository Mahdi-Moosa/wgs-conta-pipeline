"""Public commands; research-specific paths and selections are never implicit."""
import argparse,csv,json,os,sys
from pathlib import Path

def read_tsv(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path,rows,fields):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise FileExistsError('Refusing to overwrite '+str(path))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def main(argv=None):
    os.umask(0o077)
    args=list(sys.argv[1:] if argv is None else argv)
    if args and args[0]=='run':
        from .master_cli import main as run
        return run(args[1:])
    p=argparse.ArgumentParser(description='Human-sequence detection in bacterial isolate WGS')
    subs=p.add_subparsers(dest='command',required=True)
    v=subs.add_parser('validate',help='validate configuration and reference identity')
    v.add_argument('--config',required=True,type=Path);v.add_argument('--paths-only',action='store_true')
    c=subs.add_parser('classify',help='metadata-only eligibility; writes local sensitive metadata')
    c.add_argument('--input',required=True);c.add_argument('--output',required=True)
    c.add_argument('--m49-map');c.add_argument('--iso-json');c.add_argument('--taxon-deny-list')
    q=subs.add_parser('select',help='outcome-blind selection on a reviewed canonical sampling frame')
    q.add_argument('--input',required=True);q.add_argument('--m49-map',required=True)
    q.add_argument('--seed-file',required=True);q.add_argument('--projects',type=int,required=True)
    q.add_argument('--biosamples-per-project',type=int,default=3);q.add_argument('--runs-per-biosample',type=int,default=1)
    q.add_argument('--cohort-label',required=True);q.add_argument('--output',required=True)
    n=subs.add_parser('new-seed',help='create a fresh private 32-byte selection key')
    n.add_argument('--output',required=True)
    a=subs.add_parser('summarize',help='aggregate strict and completed assembly endpoints')
    a.add_argument('--input',required=True);a.add_argument('--output',required=True)
    subs.add_parser('run',help='archive acquisition, strict mapping and assembly challenge; run --help for arguments')
    a=p.parse_args(args)
    if a.command=='validate':
        from .config import load_config
        from .preflight import validate_runtime
        validate_runtime(load_config(a.config),not a.paths_only,a.config)
        print('configuration valid'+(' (resource identity not tested)' if a.paths_only else ''))
    elif a.command=='new-seed':
        from .sampling import generate_seed_file
        generate_seed_file(a.output)
    elif a.command=='classify':
        from .metadata_rules import classify,CountryReference,TaxonList,DERIVED_FIELDS
        if bool(a.m49_map)!=bool(a.iso_json):p.error('--m49-map and --iso-json must be supplied together')
        class Unresolved:
            def normalize(self,raw):return dict(country_raw=raw,sample_country_normalized='',sample_country_iso3='',sample_region_candidate='',sample_geo_status='UNRESOLVED',sample_geo_source='NOT_EVALUATED')
        countries=CountryReference(Path(a.m49_map),Path(a.iso_json)) if a.m49_map else Unresolved()
        rows=read_tsv(a.input)
        if not rows:raise ValueError('empty metadata input')
        required={'run_accession','study_accession','sample_accession','instrument_platform','library_layout','library_strategy','library_source','scientific_name','instrument_model'}
        if not required.issubset(rows[0]):raise ValueError('missing raw ENA metadata columns: '+str(sorted(required-set(rows[0]))))
        deny=TaxonList(Path(a.taxon_deny_list)) if a.taxon_deny_list else None
        result=[classify(r,countries,taxon_deny_list=deny) for r in rows]
        write_tsv(a.output,result,list(dict.fromkeys([*rows[0],*DERIVED_FIELDS])))
    elif a.command=='select':
        from .sampling import build_sampling_frame,load_m49_mapping,allocate_hamilton,select_cohort,read_key_file
        rows=read_tsv(a.input)
        required={'run_accession','canonical_bioproject','canonical_biosample','scientific_name','instrument_model','instrument_platform','sample_country_iso3','primary_eligible','country_compare_status'}
        if not rows or not required.issubset(rows[0]):raise ValueError('reviewed frame missing canonical selection fields')
        frame=build_sampling_frame(rows,load_m49_mapping(a.m49_map),a.biosamples_per_project)
        allocation=allocate_hamilton(frame.capacities(),a.projects)
        result=select_cohort(frame,allocation,a.biosamples_per_project,a.runs_per_biosample,read_key_file(a.seed_file),a.cohort_label)
        registry=[dict(Accession=r['run_accession'],Cohort=a.cohort_label,Active='yes',Priority=100,Preferred_source='auto') for r in result['membership']]
        write_tsv(a.output,registry,['Accession','Cohort','Active','Priority','Preferred_source'])
        print(json.dumps(dict(selected_runs=len(registry),selected_projects=len({r['project'] for r in result['membership']}))))
    elif a.command=='summarize':
        from .summary import aggregate
        out=Path(a.output)
        if out.exists():raise FileExistsError(out)
        out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(aggregate(read_tsv(a.input)),indent=2)+'\n')
    else:p.error('unknown command')

if __name__=='__main__':main()
