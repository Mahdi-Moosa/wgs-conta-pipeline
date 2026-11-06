from __future__ import annotations
import csv, subprocess
from collections import Counter
from pathlib import Path
from .analysis import cigar_stats, parse_nm
from .util import capture, run


def seqkit_single_stats(reads):
    txt=capture(['seqkit','stats','-T','-a',str(reads)])
    lines=[x for x in txt.splitlines() if x.strip()]
    if len(lines)<2: raise RuntimeError('Unexpected seqkit stats output for single-end reads')
    h=lines[0].split('\t'); row=dict(zip(h,lines[1].split('\t')))
    def num(v): return float(str(v).replace(',',''))
    n=int(num(row['num_seqs'])); bases=int(num(row['sum_len']))
    def val(col):
        return num(row[col]) if col in row and row[col] not in ('','-') else float('nan')
    return {'reads':n,'bases':bases,'mean_len':bases/n if n else 0,'q20':val('Q20(%)'),'q30':val('Q30(%)'),'gc':val('GC(%)')}


def human_map_single(index,reads,bam,log,threads,sort_threads):
    bam=Path(bam); log=Path(log); bam.parent.mkdir(parents=True,exist_ok=True); log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w') as lf:
        p1=subprocess.Popen(['minimap2','-ax','sr','-t',str(threads),str(index),str(reads)],stdout=subprocess.PIPE,stderr=lf)
        p2=subprocess.Popen(['samtools','sort','-@',str(sort_threads),'-o',str(bam)],stdin=p1.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        p1.stdout.close(); _,err=p2.communicate(); rc1=p1.wait()
        if rc1!=0 or p2.returncode!=0:
            raise RuntimeError(f'Ion human mapping failed minimap2={rc1} samtools={p2.returncode}: {(err or b"")[-1000:].decode(errors="ignore")}')


def strict_human_single_from_bam(bam,evidence):
    """Validated Ion branch: primary, MAPQ>=30, >=80% query coverage, >=98% identity,
    then coordinate-deduplicate passing single reads.
    """
    evidence=Path(evidence); passing=[]; mapq30=id95=id98=id99=0
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
        if cov<80 or ident<98: continue
        passing.append((qn,chrom,pos,pos+r-1,mapq,cov,ident,flag))
    _,err=p.communicate()
    if p.returncode!=0: raise RuntimeError(f'samtools view failed: {err[-2000:]}')

    # Same coordinate-dedup principle used by the validated Ion production branch.
    # Ignore strand in the coordinate key, matching the paired-fragment dedup philosophy.
    unique={}
    for rec in passing:
        qn,chrom,left,right,*_=rec
        key=(chrom,left,right)
        if key not in unique: unique[key]=rec

    evidence.mkdir(parents=True,exist_ok=True)
    qf=evidence/'strict_unique_qnames.txt'; qf.write_text(''.join(f'{v[0]}\n' for v in unique.values()))
    with (evidence/'strict_unique_reads.tsv').open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t'); w.writerow(['QNAME','Chromosome','Left','Right','MAPQ','Coverage','Identity'])
        for (_,_,_),(qn,chrom,left,right,mapq,cov,ident,flag) in unique.items():
            w.writerow([qn,chrom,left,right,mapq,f'{cov:.4f}',f'{ident:.4f}'])
    counts=Counter(k[0] for k in unique)
    with (evidence/'chromosome_counts.tsv').open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t'); w.writerow(['Chromosome','Strict_reads','Fraction']); total=len(unique)
        for c,n in counts.most_common(): w.writerow([c,n,f'{n/total:.6f}'])
    if counts: top_chr,top_n=counts.most_common(1)[0]; top_frac=top_n/len(unique)
    else: top_chr,top_n,top_frac='NA',0,0.0
    dist='ZERO' if not counts else ('SINGLE_CHROMOSOME' if len(counts)==1 else 'MULTI_CHROMOSOME')
    return {'mapq30':mapq30,'id95':id95,'id98':id98,'id99':id99,'strict_reads':len(passing),'unique':len(unique),'qnames':qf,'nchr':len(counts),'distribution':dist,'top_chr':top_chr,'top_n':top_n,'top_frac':top_frac}


def extract_strict_single(qnames,reads,evidence):
    evidence=Path(evidence); out=evidence/'strict_unique.fastq'
    run(['seqkit','grep','-f',str(qnames),str(reads),'-o',str(out)])
    return out


def kraken_strict_single(db,reads,evidence,threads):
    evidence=Path(evidence); out=evidence/'strict_kraken.out'; report=evidence/'strict_kraken.report'
    with out.open('w') as fh:
        run(['kraken2','--db',str(db),'--threads',str(threads),'--report',str(report),str(reads)],stdout=fh,stderr=subprocess.DEVNULL)
    human=total=0
    for line in out.read_text().splitlines():
        f=line.split('\t'); total+=1
        if len(f)>=3 and f[0]=='C' and f[2]=='9606': human+=1
    return human,total


def bacterial_challenge_single(fasta,reads,evidence,threads,sort_threads):
    evidence=Path(evidence); bam=evidence/'exact_bacterial.tmp.bam'; log=evidence/'bacterial_minimap2.log'
    with log.open('w') as lf:
        p1=subprocess.Popen(['minimap2','-ax','sr','-t',str(threads),str(fasta),str(reads)],stdout=subprocess.PIPE,stderr=lf)
        p2=subprocess.Popen(['samtools','sort','-@',str(sort_threads),'-o',str(bam)],stdin=p1.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        p1.stdout.close(); _,err=p2.communicate(); rc1=p1.wait()
        if rc1!=0 or p2.returncode!=0: raise RuntimeError('Ion bacterial mapping failed')
    h95=set(); h98=set()
    p=subprocess.Popen(['samtools','view','-F','0x904',str(bam)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    for line in p.stdout:
        f=line.rstrip().split('\t'); qn=f[0]; seq=f[9]; q,a,_=cigar_stats(f[5])
        if not seq or not a: continue
        cov=100*q/len(seq); ident=100*(a-parse_nm(f))/a
        if cov>=80 and ident>=95: h95.add(qn)
        if cov>=80 and ident>=98: h98.add(qn)
    _,err=p.communicate()
    if p.returncode!=0: raise RuntimeError(f'Ion bacterial samtools view failed: {err[-2000:]}')
    res={'ge1_95':len(h95),'both_95':'NA','ge1_98':len(h98),'both_98':'NA'}
    with (evidence/'bacterial_competition.tsv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=res.keys(),delimiter='\t'); w.writeheader(); w.writerow(res)
    bam.unlink(missing_ok=True); Path(str(bam)+'.bai').unlink(missing_ok=True)
    return res
