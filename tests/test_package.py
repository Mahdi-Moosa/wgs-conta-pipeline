import tempfile,unittest,subprocess,sys
from pathlib import Path
from wgs_conta.config import load_config
from wgs_conta.preflight import validate_runtime
from wgs_conta.entry import main

class PackageTests(unittest.TestCase):
    def config(self,d):
        template=Path(__file__).resolve().parents[1]/'configs/pipeline.example.ini'
        p=Path(d)/'pipeline.ini';p.write_bytes(template.read_bytes());return p
    def test_relative_config_and_missing_reference(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.config(d);s=load_config(p)
            self.assertTrue(s.project_root.is_absolute())
            self.assertTrue(validate_runtime(s,False,p))
            with self.assertRaises(FileNotFoundError):validate_runtime(s,True,p)
    def test_invalid_config(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.config(d);p.write_text(p.read_text().replace('source_mode = both','source_mode = invalid'))
            with self.assertRaises(ValueError):load_config(p)
    def test_cli_help(self):
        for args in [['--help'],['run','--help'],['select','--help'],['classify','--help']]:
            p=subprocess.run([sys.executable,'-m','wgs_conta',*args],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
    def test_cli_missing_required(self):
        p=subprocess.run([sys.executable,'-m','wgs_conta','validate'],capture_output=True)
        self.assertEqual(p.returncode,2)
    def test_seed_creation_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.key';main(['new-seed','--output',str(p)])
            self.assertEqual(p.stat().st_size,32);self.assertEqual(p.stat().st_mode&0o777,0o600)
            with self.assertRaises(FileExistsError):main(['new-seed','--output',str(p)])

    def test_resume_identity_blocks_changed_reference(self):
        import json
        from wgs_conta.preflight import guard_resume_identity
        with tempfile.TemporaryDirectory() as d:
            p=self.config(d);s=load_config(p)
            marker=s.project_root/'master/release_scientific_identity.json'
            # Isolate this test's output root instead of the example's parent workspace.
            from dataclasses import replace
            s=replace(s,project_root=Path(d)/'isolated')
            marker=s.project_root/'master/release_scientific_identity.json';marker.parent.mkdir(parents=True)
            marker.write_text(json.dumps(dict(human_index_sha256='incorrect')))
            with self.assertRaises(ValueError):guard_resume_identity(s,p,dry_run=True)
    def test_eligibility_rejects_host_workflow_and_multiple_samples(self):
        from wgs_conta.metadata_rules import classify
        class NoGeo:
            def normalize(self,x):return {}
        # Synthetic canonical strings assembled at runtime, never real archive records.
        row=dict(study_accession='PRJNA'+str(10**20),sample_accession='SAMN'+str(10**20),instrument_platform='ILLUMINA',library_layout='PAIRED',library_strategy='WGS',library_source='GENOMIC',scientific_name='Artificial bacterium',instrument_model='Synthetic instrument',strain='ARTIFICIAL_ISOLATE',sample_title='isolate sequencing')
        self.assertTrue(classify(row,NoGeo())['primary_eligible'])
        self.assertFalse(classify({**row,'sample_title':'propagation in human cells'},NoGeo())['primary_eligible'])
        self.assertFalse(classify({**row,'sample_accession':row['sample_accession']+';'+row['sample_accession']},NoGeo())['primary_eligible'])
    def test_biological_read_indices_preserve_technical_middle_read(self):
        from wgs_conta import acquisition
        from unittest.mock import patch
        from types import SimpleNamespace
        with patch.object(acquisition.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='BIOLOGICAL, TECHNICAL, BIOLOGICAL',stderr='')):
            self.assertEqual(acquisition._sra_biological_read_indices('artificial.sra'),[1,3])
