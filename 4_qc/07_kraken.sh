#!/bin/bash
# What else is in each sample besides IBV.
#
# These are clinical metagenomes, so the fraction of reads that are actually IBV - and what
# the rest are - belongs in the technical validation. It also gives an independent read on
# co-infection: a sample carrying two IBV strains still classifies as one taxon here, but
# unusually low IBV fractions are worth a second look.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
cd "$IBV_WORK"
RUN="$IBV_QC_RUN"
DB="$KRAKEN2_DB"        # PlusPF-16GB here; core_nt would need 344 GB
T=${T:-8}
K=08_kraken; mkdir -p $K/reports

while read -r strain prefix batch; do
  [ -z "$strain" ] && continue
  [ -s "$K/reports/$strain.report" ] && continue
  r1=$(find rawdata -path '*/_not_in_study/*' -prune -o \( -name "${prefix}_S*_L001_R1_001.fastq" -o -name "${prefix}_S*_L001_R1_001.fastq.gz" \) -print 2>/dev/null | sort | head -1)
  r2=$(find rawdata -path '*/_not_in_study/*' -prune -o \( -name "${prefix}_S*_L001_R2_001.fastq" -o -name "${prefix}_S*_L001_R2_001.fastq.gz" \) -print 2>/dev/null | sort | head -1)
  [ -n "$r1" ] && [ -n "$r2" ] || { echo "SKIP $strain"; continue; }
  echo "[$(date +%H:%M:%S)] $strain"
  gz=""; [[ "$r1" == *.gz ]] && gz="--gzip-compressed"
  $RUN kraken2 --db $DB --threads $T --paired $gz \
      --report $K/reports/$strain.report --output /dev/null \
      "$r1" "$r2" >> $K/kraken.log 2>&1
done < 05_reassembly/samples_all.tsv
echo "[$(date +%H:%M:%S)] kraken2 done"
