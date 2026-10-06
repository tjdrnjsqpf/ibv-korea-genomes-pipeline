#!/bin/bash
# Bring K2 back into the final set, assembled from its two libraries pooled.
# Every step is the one the pipeline applies to the other genomes.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
R="$IBV_WORK"; W=$R/38_k2; F=$R/final
RUN="$IBV_RUN"; QC="$IBV_QC_RUN"
PANEL=$R/05_reassembly/ref/panel

echo "[1] reference choice: greatest aligned span over the panel"
best=""; bestcov=0
for ref in $PANEL/*.fasta; do
  cov=$($RUN minimap2 -c -x asm20 --secondary=no "$ref" $W/both.megahit/final.contigs.fa 2>/dev/null \
        | awk '{for(i=$8;i<$9;i++) s[i]=1} END{print length(s)+0}')
  echo "    $(basename $ref .fasta)  $cov"
  [ "${cov:-0}" -gt "$bestcov" ] && { bestcov=$cov; best=$ref; }
done
bn=$(basename $best .fasta); echo "    chosen: $bn ($bestcov)"
grep -v -P "^K2\t" $F/05_reassembly/ref_choice.tsv > $W/ref_choice.tmp
printf "K2\t%s\t%s\n" "$bn" "$bestcov" >> $W/ref_choice.tmp; cp $W/ref_choice.tmp $F/05_reassembly/ref_choice.tsv

echo "[2] contigs placed on that reference (>= 1000 bp aligned)"
$RUN minimap2 -c -x asm20 --secondary=no $best $W/both.megahit/final.contigs.fa 2>/dev/null | awk '$9-$8>=1000{print "    "$1, $2" bp", "ref "$8"-"$9}'
seq=$(grep -v ">" $W/both.genome.fa | tr -d "\n"); printf ">contig1 len=%d\n%s\n" ${#seq} "$seq" > $F/05_reassembly/asm/K2.v2.fa

echo "[3] depth from the pooled reads (bwa-mem2, samtools depth -a -Q 20)"
cp $F/05_reassembly/asm/K2.v2.fa $W/K2.fa; $RUN bwa-mem2 index $W/K2.fa 2> $W/map.log
$RUN bwa-mem2 mem -t 8 $W/K2.fa $W/both_1.fq.gz $W/both_2.fq.gz 2>> $W/map.log | $RUN samtools sort -@ 4 -o $W/K2.v2.bam -
$RUN samtools index $W/K2.v2.bam
rm -f $F/05_reassembly/depth2/K2.perbase.gz
$RUN samtools depth -a -Q 20 $W/K2.v2.bam | gzip > $F/05_reassembly/depth2/K2.perbase.gz
zcat $F/05_reassembly/depth2/K2.perbase.gz | awk '{d[NR]=$3; if($3==0) z++} END{n=asort(d); printf "    %d positions, median depth %d, minimum %d, positions with no read %d\n", n, d[int(n/2)], d[1], z}'

echo "[4] library composition: kraken2 on the raw reads of both libraries"
cat $R/rawdata/Vaccine/IBV-K2-1_S3_L001_R1_001.fastq.gz $R/rawdata/Vaccine/IBV-K2-2_S7_L001_R1_001.fastq.gz > $W/raw_1.fq.gz
cat $R/rawdata/Vaccine/IBV-K2-1_S3_L001_R2_001.fastq.gz $R/rawdata/Vaccine/IBV-K2-2_S7_L001_R2_001.fastq.gz > $W/raw_2.fq.gz
$QC kraken2 --db "$KRAKEN2_DB" --threads 8 --paired --gzip-compressed \
    --report $W/K2.pooled.report --output /dev/null $W/raw_1.fq.gz $W/raw_2.fq.gz > $W/kraken.log 2>&1
rm -f $W/raw_1.fq.gz $W/raw_2.fq.gz
tail -2 $W/kraken.log
echo done
