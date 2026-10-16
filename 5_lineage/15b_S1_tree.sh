#!/bin/bash
# S1 tree, used only to assign a lineage to each genome.
#
# The framework is the 32 prototypes of the harmonised classification plus every prototype of the
# four genotypes described since it (GVII-1, GVIII, GIX-1 and GX-1 - seven sequences, because one
# reference per genotype understates identity by several points). Columns with a base in fewer than
# 90% of sequences are dropped; -cons is not used, as it would keep columns present only in the
# full-genome prototypes.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
W="$IBV_WORK/41_phylo_allproto"; cd "$W"
MM=$HOME/bin/micromamba
T=${T:-8}                # match this to the cores actually free on the machine

echo "[$(date +%H:%M)] framework alignment"
$MM run -n ibv_phylo mafft --auto --thread $T 00_data/framework_S1.fa > 01_align/framework.aln 2> logs/mafft_fw.log

echo "[$(date +%H:%M)] add our 183"
$MM run -n ibv_phylo mafft --add 00_data/ours_S1.fa --keeplength --thread $T 01_align/framework.aln > 01_align/s1_all.aln 2> logs/mafft_add.log

echo "[$(date +%H:%M)] trimAl -gt 0.9"
"${IBV_BIN:+$IBV_BIN/}trimal" -in 01_align/s1_all.aln -out 01_align/s1_all.trimal.fa -gt 0.9 > logs/trimal.log 2>&1
$MM run -n ibv_phylo python -c "
n=0;L=0
for l in open(\"01_align/s1_all.trimal.fa\"):
    if l.startswith(\">\"): n+=1; L=0
    else: L+=len(l.strip())
print(\"trimmed alignment:\",n,\"sequences x\",L,\"columns\")"

echo "[$(date +%H:%M)] IQ-TREE"
rm -f 02_tree/s1.*
$MM run -n ibv_phylo iqtree3 -s 01_align/s1_all.trimal.fa --prefix 02_tree/s1 -m MFP -B 1000 --alrt 1000 -T $T --seqtype DNA --redo > logs/iqtree.log 2>&1
grep -E "Best-fit model|Total wall-clock" 02_tree/s1.iqtree logs/iqtree.log | head -3
echo "[$(date +%H:%M)] DONE"
