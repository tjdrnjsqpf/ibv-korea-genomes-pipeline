#!/bin/bash
# Rebuild each genome de novo, then place it against whichever reference actually fits.
#
# The deposited set was scaffolded onto one Indian isolate (NC_048213.1) with RagTag, which
# both invented sequence and hid divergent regions. Two things change here: the assembly is
# de novo and never N-padded into a single contig, and the reference is chosen per strain
# from a panel covering the lineages present in this collection (QX, QXE, KM91-like,
# Mass-like, GVI-1, B4-like, GI-13). A reference is only used to select and order contigs
# and to cross-check by read depth - it never contributes sequence.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
cd "$IBV_WORK"
RUN="$IBV_RUN"
T=${T:-16}

IN=05_reassembly
PANEL=$IN/ref/panel
SAMPLES=${SAMPLES:-$IN/samples_all.tsv}
MIN_CONTIG=500
mkdir -p $IN/{asm,consensus,logs,depth2}

while read -r strain prefix batch; do
  [ -z "$strain" ] && continue
  r1=$(find rawdata -path '*/_not_in_study/*' -prune -o \( -name "${prefix}_S*_L001_R1_001.fastq" -o -name "${prefix}_S*_L001_R1_001.fastq.gz" \) -print 2>/dev/null | sort | head -1)
  r2=$(find rawdata -path '*/_not_in_study/*' -prune -o \( -name "${prefix}_S*_L001_R2_001.fastq" -o -name "${prefix}_S*_L001_R2_001.fastq.gz" \) -print 2>/dev/null | sort | head -1)
  [ -n "$r1" ] && [ -n "$r2" ] || { echo "WAIT $strain"; continue; }
  [ -s "$IN/asm/$strain.v2.fa" ] && { echo "SKIP $strain"; continue; }
  echo "[$(date +%H:%M:%S)] $strain"

  $RUN fastp -i "$r1" -I "$r2" -o $IN/asm/${strain}_1.fq.gz -O $IN/asm/${strain}_2.fq.gz \
      --detect_adapter_for_pe -w 8 -j $IN/logs/$strain.fastp.json -h /dev/null \
      2> $IN/logs/$strain.fastp.log || { echo "FAIL fastp $strain"; continue; }

  rm -rf $IN/asm/$strain.megahit
  $RUN megahit -1 $IN/asm/${strain}_1.fq.gz -2 $IN/asm/${strain}_2.fq.gz \
      -o $IN/asm/$strain.megahit -t $T --min-contig-len $MIN_CONTIG \
      > $IN/logs/$strain.megahit.log 2>&1 || { echo "FAIL megahit $strain"; continue; }
  cp $IN/asm/$strain.megahit/final.contigs.fa $IN/asm/$strain.raw.fa

  # whichever panel member the contigs align to best is this strain's lineage reference
  best=""; bestcov=0
  for ref in $PANEL/*.fasta; do
    cov=$($RUN minimap2 -c -x asm20 --secondary=no "$ref" $IN/asm/$strain.raw.fa 2>/dev/null \
          | awk '{for(i=$8;i<$9;i++) s[i]=1} END{print length(s)+0}')
    [ "${cov:-0}" -gt "$bestcov" ] && { bestcov=$cov; best=$ref; }
  done
  [ -n "$best" ] || { echo "FAIL noref $strain"; continue; }
  echo "$strain	$(basename $best .fasta)	$bestcov" >> $IN/ref_choice.tsv

  # keep contigs that explain a real stretch of that reference and add span the longer
  # contigs have not already covered - no padding, no forced single contig
  $RUN minimap2 -c -x asm20 --secondary=no "$best" $IN/asm/$strain.raw.fa 2>/dev/null \
    > $IN/asm/$strain.paf
  $RUN python - "$IN/asm/$strain.raw.fa" "$IN/asm/$strain.paf" "$IN/asm/$strain.v2.fa" <<'PY'
import sys
fa,paf,out=sys.argv[1:4]
seq={};k=None
for l in open(fa):
    if l.startswith('>'): k=l[1:].split()[0]; seq[k]=[]
    else: seq[k].append(l.strip())
seq={k:''.join(v) for k,v in seq.items()}
hits={}
for l in open(paf):
    p=l.split('\t'); hits.setdefault(p[0],[]).append((int(p[7]),int(p[8])))
covered=set(); kept=[]
for n in sorted(seq,key=lambda n:-len(seq[n])):
    span=set()
    for a,b in hits.get(n,[]): span|=set(range(a,b))
    if len(span)<1000: continue
    if len(span-covered)>=500 or not kept:
        kept.append(n); covered|=span
with open(out,'w') as fh:
    for i,n in enumerate(kept,1):
        fh.write('>contig%d len=%d\n'%(i,len(seq[n])))
        s=seq[n]
        for j in range(0,len(s),70): fh.write(s[j:j+70]+'\n')
PY

  # read depth on the rebuilt genome - any region without reads is still suspect
  $RUN bwa-mem2 index $IN/asm/$strain.v2.fa 2>>$IN/logs/$strain.map.log
  $RUN bwa-mem2 mem -t $T $IN/asm/$strain.v2.fa $IN/asm/${strain}_1.fq.gz $IN/asm/${strain}_2.fq.gz \
      2>>$IN/logs/$strain.map.log | $RUN samtools sort -@ 4 -o $IN/consensus/$strain.v2.bam -
  $RUN samtools index $IN/consensus/$strain.v2.bam
  $RUN samtools depth -a -Q 20 $IN/consensus/$strain.v2.bam | gzip > $IN/depth2/$strain.perbase.gz

  rm -f $IN/asm/${strain}_1.fq.gz $IN/asm/${strain}_2.fq.gz $IN/asm/$strain.v2.fa.[0-9a-z]*
  rm -rf $IN/asm/$strain.megahit
done < $SAMPLES

echo "[$(date +%H:%M:%S)] done"
