from __future__ import annotations
import csv,re,subprocess
from collections import Counter,defaultdict
from pathlib import Path
from .util import capture,run
CIGAR_RE=re.compile(r'(\d+)([MIDNSHP=X])')
def cigar_stats(cigar):
    q=a=r=0
    for n,op in CIGAR_RE.findall(cigar):
        n=int(n)
        if op in 'MI=X': q+=n
        if op in 'MID=X': a+=n
        if op in 'MDN=X': r+=n
    return q,a,r
def parse_nm(f):
    for x in f[11:]:
        if x.startswith('NM:i:'): return int(x.split(':')[2])
    return 0
def seqkit_pair_stats(r1,r2):
    txt=capture(['seqkit','stats','-T','-a',str(r1),str(r2)]); lines=[x for x in txt.splitlines() if x.strip()]
    if len(lines)<3: raise RuntimeError('Unexpected seqkit stats output')
    h=lines[0].split('\t'); rows=[dict(zip(h,x.split('\t'))) for x in lines[1:3]]
    def num(v): return float(str(v).replace(',',''))
    n=[int(num(x['num_seqs'])) for x in rows]
    if n[0]!=n[1]: raise RuntimeError(f'Unequal R1/R2 counts: {n}')
    bases=sum(int(num(x['sum_len'])) for x in rows)
    def weighted(col):
        vals=[]; weights=[]
        for x in rows:
            if col in x and x[col] not in ('','-'): vals.append(num(x[col])); weights.append(num(x['sum_len']))
        return sum(v*w for v,w in zip(vals,weights))/sum(weights) if vals else float('nan')
    return {'pairs':n[0],'bases':bases,'mean_len':bases/sum(n) if sum(n) else 0,'q20':weighted('Q20(%)'),'q30':weighted('Q30(%)'),'gc':weighted('GC(%)')}
def human_map(index,r1,r2,bam,log,threads,sort_threads):
    bam=Path(bam); log=Path(log); bam.parent.mkdir(parents=True,exist_ok=True); log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w') as lf:
        p1=subprocess.Popen(['minimap2','-ax','sr','-t',str(threads),str(index),str(r1),str(r2)],stdout=subprocess.PIPE,stderr=lf)
        p2=subprocess.Popen(['samtools','sort','-@',str(sort_threads),'-o',str(bam)],stdin=p1.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        p1.stdout.close(); _,err=p2.communicate(); rc1=p1.wait()
        if rc1!=0 or p2.returncode!=0: raise RuntimeError(f'human mapping failed minimap2={rc1} samtools={p2.returncode}: {(err or b"")[-1000:].decode(errors="ignore")}')
def strict_human_from_bam(bam,evidence):
    evidence=Path(evidence); passing=defaultdict(dict); mapq30=id95=id98=id99=0
    p=subprocess.Popen(['samtools','view','-F','0x904',str(bam)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    for line in p.stdout:
        f=line.rstrip().split('\t'); qn=f[0]; flag=int(f[1]); chrom=f[2]; pos=int(f[3]); mapq=int(f[4]); seq=f[9]
        if mapq<30: continue
        mapq30+=1; q,a,r=cigar_stats(f[5])
        if not seq or not a: continue
        cov=100*q/len(seq); ident=100*(a-parse_nm(f))/a
        if cov>=80 and ident>=95: id95+=1
        if cov>=80 and ident>=98: id98+=1
        if cov>=80 and ident>=99: id99+=1
        if cov<80 or ident<98 or not(flag&0x2): continue
        mate=1 if flag&0x40 else (2 if flag&0x80 else 0)
        if mate: passing[qn][mate]={'chrom':chrom,'pos':pos,'end':pos+r-1,'mapq':mapq,'cov':cov,'ident':ident}
    _,err=p.communicate()
    if p.returncode!=0: raise RuntimeError(f'samtools view failed: {err[-2000:]}')
    complete=0; unique={}
    for qn,m in passing.items():
        if 1 not in m or 2 not in m: continue
        r1,r2=m[1],m[2]
        if r1['chrom']!=r2['chrom']: continue
        complete+=1; key=(r1['chrom'],min(r1['pos'],r2['pos']),max(r1['end'],r2['end']))
        if key not in unique: unique[key]=(qn,r1,r2)
    evidence.mkdir(parents=True,exist_ok=True); qf=evidence/'strict_unique_qnames.txt'; qf.write_text(''.join(f'{v[0]}\n' for v in unique.values()))
    with (evidence/'strict_unique_fragments.tsv').open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t'); w.writerow(['QNAME','Chromosome','Left','Right','R1_MAPQ','R1_Coverage','R1_Identity','R2_MAPQ','R2_Coverage','R2_Identity'])
        for (chrom,left,right),(qn,r1,r2) in unique.items(): w.writerow([qn,chrom,left,right,r1['mapq'],f"{r1['cov']:.4f}",f"{r1['ident']:.4f}",r2['mapq'],f"{r2['cov']:.4f}",f"{r2['ident']:.4f}"])
    counts=Counter(k[0] for k in unique)
    with (evidence/'chromosome_counts.tsv').open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t'); w.writerow(['Chromosome','Strict_fragments','Fraction']); total=len(unique)
        for c,n in counts.most_common(): w.writerow([c,n,f'{n/total:.6f}'])
    if counts: top_chr,top_n=counts.most_common(1)[0]; top_frac=top_n/len(unique)
    else: top_chr,top_n,top_frac='NA',0,0.0
    dist='ZERO' if not counts else ('SINGLE_CHROMOSOME' if len(counts)==1 else 'MULTI_CHROMOSOME')
    return {'mapq30':mapq30,'id95':id95,'id98':id98,'id99':id99,'strict_pairs':complete,'unique':len(unique),'qnames':qf,'nchr':len(counts),'distribution':dist,'top_chr':top_chr,'top_n':top_n,'top_frac':top_frac}
def extract_strict(qnames,r1,r2,evidence):
    evidence=Path(evidence); u1=evidence/'strict_unique_R1.fastq'; u2=evidence/'strict_unique_R2.fastq'
    run(['seqkit','grep','-f',str(qnames),str(r1),'-o',str(u1)]); run(['seqkit','grep','-f',str(qnames),str(r2),'-o',str(u2)]); return u1,u2
def kraken_strict(db,r1,r2,evidence,threads):
    evidence=Path(evidence); out=evidence/'strict_kraken.out'; report=evidence/'strict_kraken.report'
    with out.open('w') as fh: run(['kraken2','--db',str(db),'--threads',str(threads),'--paired','--report',str(report),str(r1),str(r2)],stdout=fh,stderr=subprocess.DEVNULL)
    human=total=0
    for line in out.read_text().splitlines():
        f=line.split('\t'); total+=1
        if len(f)>=3 and f[0]=='C' and f[2]=='9606': human+=1
    return human,total
def bacterial_challenge(fasta,r1,r2,evidence,threads,sort_threads):
    evidence=Path(evidence); bam=evidence/'exact_bacterial.tmp.bam'; log=evidence/'bacterial_minimap2.log'
    with log.open('w') as lf:
        p1=subprocess.Popen(['minimap2','-ax','sr','-t',str(threads),str(fasta),str(r1),str(r2)],stdout=subprocess.PIPE,stderr=lf)
        p2=subprocess.Popen(['samtools','sort','-@',str(sort_threads),'-o',str(bam)],stdin=p1.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        p1.stdout.close(); _,err=p2.communicate(); rc1=p1.wait()
        if rc1!=0 or p2.returncode!=0: raise RuntimeError('bacterial mapping failed')
    h95=defaultdict(set); h98=defaultdict(set); p=subprocess.Popen(['samtools','view','-F','0x904',str(bam)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    for line in p.stdout:
        f=line.rstrip().split('\t'); qn=f[0]; flag=int(f[1]); seq=f[9]; q,a,_=cigar_stats(f[5])
        if not seq or not a: continue
        cov=100*q/len(seq); ident=100*(a-parse_nm(f))/a; mate=1 if flag&0x40 else (2 if flag&0x80 else 0)
        if cov>=80 and ident>=95: h95[qn].add(mate)
        if cov>=80 and ident>=98: h98[qn].add(mate)
    _,err=p.communicate()
    if p.returncode!=0: raise RuntimeError(f'bacterial samtools view failed: {err[-2000:]}')
    res={'ge1_95':len(h95),'both_95':sum(1 for x in h95.values() if 1 in x and 2 in x),'ge1_98':len(h98),'both_98':sum(1 for x in h98.values() if 1 in x and 2 in x)}
    with (evidence/'bacterial_competition.tsv').open('w',newline='') as fh: w=csv.DictWriter(fh,fieldnames=res.keys(),delimiter='\t'); w.writeheader(); w.writerow(res)
    bam.unlink(missing_ok=True); Path(str(bam)+'.bai').unlink(missing_ok=True); return res
