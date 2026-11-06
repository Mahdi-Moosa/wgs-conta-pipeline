"""Wholly artificial technical-failure and downstream-status controls."""
import io,json,subprocess,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
from wgs_conta import analysis as A,worker as W
from wgs_conta.summary import aggregate
from wgs_conta.master import MasterStore

RECORD='\t'.join(map(str,['ARTIFICIAL',99,'chr1',1,60,'100M','=',101,0,'A'*100,'I'*100,'NM:i:0']))+'\n'
class Proc:
    def __init__(self,text='',code=0):self.stdout=io.StringIO(text);self.returncode=code
    def communicate(self):return '', 'ARTIFICIAL_IO_ERROR' if self.returncode else ''
    def wait(self):return self.returncode

def processes(text='',code=0,mapper=0,sorter=0):
    def popen(args,**kwargs):
        if args[0]=='minimap2':return Proc(code=mapper)
        if args[1]=='sort':return Proc(code=sorter)
        return Proc(text,code)
    return popen

class AssemblyFailureTests(unittest.TestCase):
    def challenge(self,text='',code=0,**kw):
        with tempfile.TemporaryDirectory() as d,patch.object(A.subprocess,'Popen',side_effect=processes(text,code,**kw)):
            return A.bacterial_challenge('ARTIFICIAL.fa','A.fq','B.fq',d,1,1)
    def test_successful_evidence_unchanged(self):self.assertEqual(self.challenge(RECORD)['ge1_95'],1)
    def test_successful_empty_stream_is_legitimate_zero(self):self.assertEqual(self.challenge()['ge1_95'],0)
    def test_failed_empty_stream_raises(self):
        with self.assertRaises(RuntimeError):self.challenge(code=1)
    def test_partial_failed_stream_raises(self):
        with self.assertRaises(RuntimeError):self.challenge(RECORD,1)
    def test_malformed_stream_not_accepted(self):
        with self.assertRaises((ValueError,IndexError,RuntimeError)):self.challenge('MALFORMED\n',1)
    def test_missing_executable_not_accepted(self):
        with tempfile.TemporaryDirectory() as d,patch.object(A.subprocess,'Popen',side_effect=FileNotFoundError('ARTIFICIAL_MISSING_TOOL')):
            with self.assertRaises(FileNotFoundError):A.bacterial_challenge('A.fa','A.fq','B.fq',d,1,1)
    def test_mapper_and_sorter_pipe_failures_raise(self):
        for kw in [dict(mapper=1),dict(sorter=1)]:
            with self.subTest(**kw),self.assertRaises(RuntimeError):self.challenge(**kw)
    def test_real_tools_missing_input_raises(self):
        import shutil
        if not (shutil.which('samtools') and shutil.which('minimap2')):self.skipTest('external tools missing')
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):A.bacterial_challenge(Path(d)/'missing.fa',Path(d)/'missing_R1.fq',Path(d)/'missing_R2.fq',d,1,1)
    def test_worker_master_resume_and_aggregation_preserve_failed_state(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);s=SimpleNamespace(evidence_root=p/'evidence',human_index='ARTIFICIAL.mmi',threads_per_worker=1,sort_threads=1,kraken_db='ARTIFICIAL',kraken_threads=1)
            ac=SimpleNamespace(metadata=dict(Platform='ILLUMINA',LibraryLayout='PAIRED',BioSample='ARTIFICIAL_SAMPLE',BioProject='ARTIFICIAL_PROJECT'),r1=p/'R1.fq',r2=p/'R2.fq',source='ENA',fmt='FASTQ',bytes_downloaded=0,read_filter_policy='AS_PROVIDED',reject_pairs_removed=0)
            def strict(bam,e):
                (e/'chromosome_counts.tsv').write_text('Chromosome\tStrict_fragments\tFraction\nchr1\t1\t1\n')
                return dict(unique=1,qnames=e/'q.txt',nchr=1,distribution='SINGLE_CHROMOSOME',top_chr='chr1',top_n=1,top_frac=1,mapq30=2,id95=2,id98=2,id99=2,strict_pairs=1)
            def extracts(*args):
                x=p/'R1.fq';y=p/'R2.fq';x.write_text('ARTIFICIAL');y.write_text('ARTIFICIAL');return x,y
            stats=dict(pairs=10,bases=2000,mean_len=100,q20=100,q30=100,gc=50)
            with patch.object(W,'seqkit_pair_stats',return_value=stats),patch.object(W,'human_map'),patch.object(W,'strict_human_from_bam',side_effect=strict),patch.object(W,'extract_strict',side_effect=extracts),patch.object(W,'kraken_strict',return_value=(0,1)),patch.object(W,'_assembly_for_positive',return_value=(SimpleNamespace(accession='ARTIFICIAL_ASSEMBLY',source='TEST',lookup_status='COMPLETE'),'A.fa',None)),patch.object(A.subprocess,'Popen',side_effect=processes('',1)):
                import time
                row=W.analyze_paired(s,'ena','ARTIFICIAL','ARTIFICIAL',p,ac,time.time())
            self.assertEqual(row['Assembly_check'],'FAILED');self.assertEqual(row['Bact_ge1_ID95'],'NA')
            store=MasterStore(p/'master','0.1.0rc1','2.3.0');store.commit(row,dict(cohort='ARTIFICIAL'),p/'evidence/ARTIFICIAL')
            # Resume means strict processing was committed, not a successful assembly check.
            self.assertTrue(store.is_complete('ARTIFICIAL'))
            import csv
            with store.results.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(rows[0]['Assembly_check'],'FAILED')
            result=aggregate(rows);self.assertEqual(result['positive_runs_with_completed_challenge'],0)
            self.assertEqual(result['positive_runs_without_detected_bacterial_match'],0)
            self.assertEqual(result['assembly_status_counts']['FAILED'],1)
