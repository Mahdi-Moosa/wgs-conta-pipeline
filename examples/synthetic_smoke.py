"""Run only fabricated data through the unchanged strict-call predicate."""
import io,json,tempfile
from unittest.mock import patch
from wgs_conta.analysis import strict_human_from_bam
class Fake:
    def __init__(self,text):self.stdout=io.StringIO(text);self.returncode=0
    def communicate(self):return '', ''
def row(flag,pos):return '\t'.join(map(str,['SYNTHETIC_PAIR',flag,'chr1',pos,30,'100M','=',0,0,'A'*100,'I'*100,'NM:i:0']))+'\n'
with tempfile.TemporaryDirectory() as tmp:
    with patch('wgs_conta.analysis.subprocess.Popen',return_value=Fake(row(99,100)+row(147,200))):
        r=strict_human_from_bam('artificial.bam',tmp)
        assert r['unique']==1 and r['nchr']==1
print(json.dumps(dict(status='PASS',synthetic_strict_fragments=1,external_reference_used=False)))
