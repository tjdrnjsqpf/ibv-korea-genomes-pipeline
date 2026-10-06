#!/bin/bash
# Tree on S1 columns only: keep columns with a base in at least 90% of sequences.
# (The earlier run also passed -cons 60, which forces 60% of all columns to be kept
#  and so retained ~15 kb of columns present only in two full-genome prototypes.)
N="$IBV_WORK/27_phylo_final"; cd "$N"
"${IBV_BIN:+$IBV_BIN/}trimal" -in 01_align/s1_all.aln -out 01_align/s1_all.trimal.fa -gt 0.9 > logs/trimal.log 2>&1
python3 -c "
n=0;L=0
for l in open('01_align/s1_all.trimal.fa'):
    if l.startswith('>'): n+=1; L=0
    else: L+=len(l.strip())
print('trimmed alignment:',n,'sequences x',L,'columns')"
rm -f 02_tree/s1.*
$IBV_PHYLO_RUN iqtree3 -s 01_align/s1_all.trimal.fa --prefix 02_tree/s1 -m MFP -B 1000 --alrt 1000 -T 8 --seqtype DNA --redo > logs/iqtree.log 2>&1
$IBV_PHYLO_RUN python assign_lineage.py
grep -E "Best-fit model|Total wall-clock" 02_tree/s1.iqtree logs/iqtree.log | head -3
echo "DONE $(date +%H:%M)"
