#!/usr/bin/env python3
"""Identity of each genome's S1 to the prototype of the lineage it was assigned.

A tree can only place a genome among the lineages its framework contains, and support values do
not reveal a missing one: four genomes of this collection sat on a 97/100 node inside GII-1 while
being only 72.4-72.8% identical to that prototype, because the genotype they belong to had not yet
been described. This script reports the identity so that such a case cannot pass unnoticed. Where a
lineage has more than one prototype the best match is used, since a single reference understates
identity by several points.

The comparison is meaningful within a lineage, not against a fixed cut-off: lineages differ in how
tightly their members sit around a prototype.

usage: 15c_lineage_identity.py <alignment.fasta> <taxa.tsv> <lineage_assignment.tsv> [out.tsv]
"""
import csv, sys, collections

aln_path, taxa_path, assign_path = sys.argv[1:4]
out_path = sys.argv[4] if len(sys.argv) > 4 else "-"

aln = {}
for blk in open(aln_path).read().split(">")[1:]:
    head, *body = blk.splitlines()
    aln[head.split()[0]] = "".join(body).upper()

proto = collections.defaultdict(list)
for r in csv.DictReader(open(taxa_path), delimiter="\t"):
    if r["set"].startswith("lineage prototype") and r["label"]:
        proto[r["label"]].append(r["taxon"])

def identity(a, b):
    """percent identity over the columns where both have a base"""
    n = m = 0
    for x, y in zip(aln[a], aln[b]):
        if x == "-" or y == "-": continue
        n += 1; m += (x == y)
    return 100 * m / n if n else 0.0

rows = []
for r in csv.DictReader(open(assign_path), delimiter="\t"):
    s, lin = r["sample"], r["lineage"]
    refs = [p for p in proto.get(lin, []) if p in aln]
    if lin == "unresolved" or not refs or s not in aln:
        rows.append((s, lin, "", ""))
        continue
    best = max(refs, key=lambda p: identity(s, p))
    rows.append((s, lin, "%.1f" % identity(s, best), best))

fh = sys.stdout if out_path == "-" else open(out_path, "w", newline="")
w = csv.writer(fh, delimiter="\t")
w.writerow(["sample", "lineage", "identity_pct", "closest_prototype"])
w.writerows(rows)
if fh is not sys.stdout: fh.close()

vals = collections.defaultdict(list)
for s, lin, v, _ in rows:
    if v: vals[lin].append(float(v))
for lin, v in sorted(vals.items()):
    print("%-8s n=%-4d %5.1f%% - %5.1f%%  median %5.1f%%"
          % (lin, len(v), min(v), max(v), sorted(v)[len(v) // 2]), file=sys.stderr)
