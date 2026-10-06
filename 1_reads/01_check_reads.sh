#!/bin/bash
# Check what actually landed before anything reads it.
#
# rsync --partial leaves a half-written file under its real name if a transfer is cut off,
# and a truncated fastq assembles without complaining, so every file is checked for a
# clean gzip stream / whole record count and R1-R2 agreement before it is trusted.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
cd "$IBV_WORK"
IN=05_reassembly
OUT=$IN/inventory.tsv

printf "genome\tread_prefix\tbatch\tstatus\tr1_reads\tr2_reads\tr1_file\n" > $OUT
ok=0; bad=0; missing=0
while read -r sid prefix batch; do
  [ -z "$sid" ] && continue
  r1=$(find rawdata -name "${prefix}_S*_L001_R1_001.fastq" -o -name "${prefix}_S*_L001_R1_001.fastq.gz" 2>/dev/null | sort | head -1)
  r2=$(find rawdata -name "${prefix}_S*_L001_R2_001.fastq" -o -name "${prefix}_S*_L001_R2_001.fastq.gz" 2>/dev/null | sort | head -1)
  if [ -z "$r1" ] || [ -z "$r2" ]; then
    printf "%s\t%s\t%s\tMISSING\t\t\t\n" "$sid" "$prefix" "$batch" >> $OUT
    missing=$((missing+1)); continue
  fi
  # a gzip stream that does not decompress cleanly means the transfer was cut off
  count() {
    if [[ "$1" == *.gz ]]; then gzip -t "$1" 2>/dev/null || { echo BAD; return; }
                                echo $(( $(gzip -dc "$1" | wc -l) / 4 ))
    else echo $(( $(wc -l < "$1") / 4 )); fi
  }
  a=$(count "$r1"); b=$(count "$r2")
  if [ "$a" = "BAD" ] || [ "$b" = "BAD" ] || [ "$a" != "$b" ] || [ "$a" -eq 0 ] 2>/dev/null; then
    printf "%s\t%s\t%s\tCORRUPT\t%s\t%s\t%s\n" "$sid" "$prefix" "$batch" "$a" "$b" "$(basename "$r1")" >> $OUT
    bad=$((bad+1))
  else
    printf "%s\t%s\t%s\tOK\t%s\t%s\t%s\n" "$sid" "$prefix" "$batch" "$a" "$b" "$r1" >> $OUT
    ok=$((ok+1))
  fi
done < $IN/samples_all.tsv

echo "ok $ok / corrupt $bad / missing $missing  (total $((ok+bad+missing)))"
echo
echo "=== by batch ==="
awk -F'\t' 'NR>1{c[$3"\t"$4]++} END{for(k in c) print k"\t"c[k]}' $OUT | sort | awk -F'\t' '{printf "  %-18s %-8s %s\n",$1,$2,$3}'
[ $bad -gt 0 ] && { echo; echo "=== corrupt files, re-transfer needed ==="; awk -F'\t' '$4=="CORRUPT"{print "  "$1"\t"$7}' $OUT; }
echo
echo "-> $OUT"
