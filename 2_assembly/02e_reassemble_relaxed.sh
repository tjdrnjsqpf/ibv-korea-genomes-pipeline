#!/bin/bash
# Try to close the four real gaps from the reads alone: re-assemble with singleton
# k-mers kept and a finer k ladder, then see whether any contig spans the gap.
R="$IBV_WORK"
IN=$R/05_reassembly; W=$R/33_gapclose; mkdir -p $W
RUN="$IBV_RUN"
MM2="${IBV_BIN:+$IBV_BIN/}minimap2"
for spec in 19ND70:MT978193.1:25503:25565 20ND039CT:MT978193.1:25500:25544 Kr-D15-15-CE7:PZ330336.1:13423:13711 20L215F:PQ387070.1:8178:9070; do
  IFS=: read s ref a b <<< "$spec"
  raw1=$(ls $R/rawdata/*/${s}_*R1*.fastq.gz 2>/dev/null | head -1)
  [ -z "$raw1" ] && raw1=$(ls $R/rawdata/*/${s}*R1*.fastq.gz 2>/dev/null | head -1)
  if [ -z "$raw1" ]; then echo "$s: raw reads not found"; continue; fi
  raw2=${raw1/R1/R2}
  r1=$W/${s}_1.fq.gz; r2=$W/${s}_2.fq.gz
  $RUN fastp -i $raw1 -I $raw2 -o $r1 -O $r2 --detect_adapter_for_pe -w 8 -j $W/$s.fastp.json -h /dev/null > $W/$s.fastp.log 2>&1
  rm -rf $W/$s.megahit
  $RUN megahit -1 $r1 -2 $r2 -o $W/$s.megahit -t 8 --min-count 1 --k-min 21 --k-max 141 --k-step 10 \
       --min-contig-len 300 > $W/$s.log 2>&1 || { echo "$s: megahit failed"; continue; }
  echo "=== $s  gap $ref:$a-$b  (relaxed assembly: $(grep -c '>' $W/$s.megahit/final.contigs.fa) contigs) ==="
  $MM2 -x asm20 --secondary=no $R/31_fixable/$ref.fa $W/$s.megahit/final.contigs.fa 2>/dev/null \
    | awk -v a=$a -v b=$b '$8<b+300 && $9>a-300 {s=($8<=a-20 && $9>=b+20)?"SPANS THE GAP":"touches"; printf "   %-14s len=%-6d ref %6d-%-6d  %s\n",$1,$2,$8,$9,s}'
done
echo done
