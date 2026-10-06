#!/bin/bash
# Run each genome through table2asn, the validator GenBank actually applies on submission.
#
# It reads the .fsa/.tbl pair and reports what would block or flag the record: internal
# stops, bad start/stop codons, coordinates past the sequence end, short CDS. The .val
# output is what decides which of the flagged genomes are real problems rather than
# ordinary strain variation.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
cd "$IBV_WORK"
PY="$IBV_RUN python"   # system python3 predates f-strings
T2A="$IBV_TABLE2ASN"
V=10_validation
rm -rf $V; mkdir -p $V/in $V/out

# a minimal submission template; table2asn refuses to start without one
cat > $V/template.sbt <<'SBT'
Submit-block ::= {
  contact {
    contact {
      name name { last "Lee", first "Sung-Gwon" },
      affil std { affil "Chonnam National University", country "South Korea" }
    }
  },
  cit {
    authors {
      names std { { name name { last "Lee", first "Sung-Gwon" } } },
      affil std { affil "Chonnam National University", country "South Korea" }
    }
  },
  subtype new
}
SBT

# the annotator already emits oriented fasta matching the table coordinates
cp 09_annotation/fsa/*.fsa $V/in/
cp 09_annotation/tbl/*.tbl $V/in/

echo "[$(date +%H:%M:%S)] validating $(ls $V/in/*.fsa 2>/dev/null | wc -l) files"
$T2A -indir $V/in -outdir $V/out -t $V/template.sbt -V v -M n -Z > $V/table2asn.log 2>&1
echo "[$(date +%H:%M:%S)] done · .val $(ls $V/out/*.val 2>/dev/null | wc -l) files"
