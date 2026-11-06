"""Artificial reads and metadata only. No reference-derived human sequences."""
import csv,io,json,subprocess,tempfile,unittest,shutil
from pathlib import Path
from unittest.mock import patch
from wgs_conta import analysis as A,ion_analysis as I,acquisition as Q,db
from wgs_conta.master import MasterStore,parse_registry
from wgs_conta.summary import aggregate
from wgs_conta.sampling import build_sampling_frame,allocate_hamilton,select_cohort
from wgs_conta import archive

class Stream:
    def __init__(self,text):self.stdout=io.StringIO(text);self.returncode=0
    def communicate(self):return '', ''
    def wait(self):return 0

# These records are artificial; qualities/bases are uniform, not biological data.
def sam(q,flag,pos,mq=30,cigar='100M',nm=0,chrom='chr1'):
    return '\t'.join(map(str,[q,flag,chrom,pos,mq,cigar,'=',pos+100,0,'A'*100,'I'*100,'NM:i:'+str(nm)]))+'\n'
def pair(q,pos,**kw):return sam(q,99,pos,**kw)+sam(q,147,pos+100,**kw)

class ScientificTests(unittest.TestCase):
    def test_cigar_and_identity_definition(self):
        self.assertEqual(A.cigar_stats('20S75M2I3M1D'),(80,81,79))
        self.assertEqual(A.parse_nm(['']*11+['NM:i:2']),2)
    def test_paired_boundaries_and_coordinate_dedup(self):
        text=pair('pass',100)+pair('duplicate',100)+pair('at_identity',400,nm=2)+pair('at_coverage',700,cigar='20S80M',nm=1)+pair('lowmq',1000,mq=29)+pair('lowcov',1300,cigar='21S79M')+pair('lowidentity',1600,nm=3)+sam('unpaired',65,1900)+sam('unpaired',129,2000)+sam('different',99,2200)+sam('different',147,2300,chrom='chr2')
        with tempfile.TemporaryDirectory() as d,patch.object(A.subprocess,'Popen',return_value=Stream(text)) as pop:
            r=A.strict_human_from_bam('synthetic.bam',d)
            self.assertEqual(r['unique'],3);self.assertEqual(r['strict_pairs'],4)
            self.assertEqual(pop.call_args.args[0][:4],['samtools','view','-F','0x904'])
            with open(Path(d)/'strict_unique_fragments.tsv') as handle:
                self.assertEqual(len(list(csv.DictReader(handle,delimiter='\t'))),3)
    def test_single_read_boundaries_and_strand_dedup(self):
        text=sam('one',0,100)+sam('duplicate_reverse',16,100)+sam('identity',0,400,nm=2)+sam('fail',0,700,nm=3)
        with tempfile.TemporaryDirectory() as d,patch.object(I.subprocess,'Popen',return_value=Stream(text)):
            r=I.strict_human_single_from_bam('synthetic.bam',d)
            self.assertEqual(r['unique'],2);self.assertEqual(r['strict_reads'],3)
    def test_platform_layout_routes(self):
        for p in ['ILLUMINA','BGISEQ','DNBSEQ']:
            self.assertTrue(Q._supported_combo({'Platform':p,'LibraryLayout':'PAIRED'}))
            self.assertFalse(Q._supported_combo({'Platform':p,'LibraryLayout':'SINGLE'}))
        self.assertTrue(Q._supported_combo({'Platform':'ION_TORRENT','LibraryLayout':'SINGLE'}))
        self.assertFalse(Q._supported_combo({'Platform':'ION_TORRENT','LibraryLayout':'PAIRED'}))
        self.assertFalse(Q._supported_combo({'Platform':'OXFORD_NANOPORE','LibraryLayout':'SINGLE'}))
    def test_bacterial_challenge_requires_either_mate_at95(self):
        records=sam('candidate',99,100,nm=5)+sam('candidate',147,200,nm=6)
        def fake_run(*a,**kw):return Stream(records)
        with tempfile.TemporaryDirectory() as d,patch.object(A.subprocess,'Popen',side_effect=fake_run),patch.object(A,'run'):
            r=A.bacterial_challenge('artificial.fa','r1','r2',d,1,1)
            self.assertEqual(r,dict(ge1_95=1,both_95=0,ge1_98=0,both_98=0))
    def test_missing_assembly_not_silently_negative(self):
        rows=[dict(Accession='SYNTHETIC_A',Human_positive='YES',Assembly_check='FAILED',Quantification_unit='PAIRS',Unique_strict_units='1'),dict(Accession='SYNTHETIC_B',Human_positive='NO',Assembly_check='NOT_NEEDED_NEGATIVE',Quantification_unit='PAIRS',Unique_strict_units='0')]
        r=aggregate(rows)
        self.assertEqual(r['strict_positive_runs'],1);self.assertEqual(r['positive_runs_with_completed_challenge'],0)
        self.assertEqual(r['positive_runs_without_detected_bacterial_match'],0)
        with self.assertRaises(ValueError):aggregate(rows+rows)
    def test_master_commit_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);e=root/'evidence';e.mkdir()
            (e/'chromosome_counts.tsv').write_text('Chromosome\tStrict_fragments\tFraction\nchr1\t1\t1\n')
            s=MasterStore(root/'master','0.1.0rc1','2.3.0')
            row=dict(Accession='SYNTHETIC_A',Analysis_branch='PAIRED_SHORTREAD',Unique_strict_fragments=1,Unique_strict_units=1,Human_positive='YES',Assembly_check='PASS')
            s.commit(row,dict(cohort='ARTIFICIAL',notes=''),e)
            self.assertTrue(s.is_complete('SYNTHETIC_A'))
            old=s.results.read_bytes();s.commit(row,dict(cohort='ARTIFICIAL',notes=''),e)
            self.assertEqual(s.results.read_bytes(),old)
            s.chromosomes.write_text(s.chromosomes.read_text().splitlines()[0]+'\n')
            self.assertFalse(s.is_complete('SYNTHETIC_A'))
    def test_runtime_failed_jobs_requeue(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'state.sqlite';db.init_db(p)
            db.prepare_master_jobs(p,[('SYNTHETIC_A','TEST','ena',100)],'0.1.0rc1')
            job=db.claim_job(p,'ena','test');self.assertEqual(job['accession'],'SYNTHETIC_A')
            db.mark_analysis_error(p,'SYNTHETIC_A','ena','artificial failure')
            db.prepare_master_jobs(p,[('SYNTHETIC_A','TEST','ena',100)],'0.1.0rc1')
            self.assertIsNotNone(db.claim_job(p,'ena','test-again'))
    def test_half_open_archive_bounds(self):
        q=archive.range_query('ILLUMINA','2020-01-01','2020-02-01')
        self.assertIn('first_public<2020-02-01',q);self.assertNotIn('first_public<=',q)
        cfg=dict(platforms=list(archive.LAYOUT),start_date='2020-01-01',cutoff_date='2020-01-31')
        archive.validate_plan(archive.initial_plan(cfg),cfg)
    def test_selection_is_input_order_invariant(self):
        rows=[]
        for p in range(3):
            for b in range(3):
                rows.append(dict(run_accession=f'SYNTHETIC_{p}_{b}',canonical_bioproject=f'PROJECT_{p}',canonical_biosample=f'SAMPLE_{p}_{b}',primary_eligible=True,country_compare_status='MATCH',sample_country_iso3='TST',scientific_name='Artificial bacterium',instrument_model='Synthetic instrument',instrument_platform='ILLUMINA'))
        f=build_sampling_frame(rows,{'TST':('Synthetic subregion','Synthetic region')})
        g=build_sampling_frame(reversed(rows),{'TST':('Synthetic subregion','Synthetic region')})
        key=b'artificial-test-key-only-0000000' # synthetic test key, not a historical seed
        self.assertEqual(len(key),32)
        a=allocate_hamilton(f.capacities(),2)
        x=select_cohort(f,a,seed=key,cohort_label='TEST')['membership']
        y=select_cohort(g,a,seed=key,cohort_label='TEST')['membership']
        self.assertEqual(x,y);self.assertEqual(len(x),6)

class ExternalSamtoolsSmoke(unittest.TestCase):
    @unittest.skipUnless(shutil.which('samtools'),'samtools not installed; parser tests still run')
    def test_primary_filter_with_real_samtools(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);s=p/'artificial.sam';bam=p/'artificial.bam'
            s.write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n'+pair('normal',100)+sam('secondary',99|256,400)+sam('secondary',147|256,500)+sam('supplementary',99|2048,700)+sam('supplementary',147|2048,800))
            subprocess.run(['samtools','view','-b','-o',str(bam),str(s)],check=True,capture_output=True)
            r=A.strict_human_from_bam(bam,p/'out');self.assertEqual(r['unique'],1)

    @unittest.skipUnless(shutil.which('minimap2') and shutil.which('samtools'),'external mapping tools unavailable')
    def test_real_synthetic_exact_assembly_alignment(self):
        import random
        rng=random.Random('artificial-assembly-fixture-only')
        reference=''.join(rng.choice('ACGT') for _ in range(20000))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);fa=p/'artificial_bacterium.fa'
            fa.write_text('>artificial_scaffold\n'+reference+'\n')
            reads=[reference[1000:1150],reference[1250:1400].translate(str.maketrans('ACGT','TGCA'))[::-1]]
            for mate,seq in enumerate(reads,1):
                (p/f'R{mate}.fastq').write_text('@ARTIFICIAL_PAIR\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
            result=A.bacterial_challenge(fa,p/'R1.fastq',p/'R2.fastq',p,1,1)
            self.assertEqual(result,dict(ge1_95=1,both_95=1,ge1_98=1,both_98=1))

if __name__=='__main__':unittest.main()
