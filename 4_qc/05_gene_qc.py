#!/usr/bin/env python
"""Judge each CDS against the dataset's own consensus length, not a single GenBank record.

Comparing to whichever panel reference a strain matched made the verdict depend on how well
that record happened to be annotated: 3c is given as 282 nt in the B4-like reference and 333
in the QXE one, a 51 nt spread for the same gene, so 14 perfectly ordinary genomes were
flagged purely for having matched the short one. The reference length is a poor standard.

The dataset's own median per gene is a better one, and for a data descriptor it is also the
right thing to report - it says how each genome sits relative to the collection it belongs to.
"""
import sys, os, glob, csv, statistics, collections

base = sys.argv[1]
TBL  = f"{base}/09_annotation/tbl"
FSA  = f"{base}/09_annotation/fsa"
OUT  = f"{base}/09_annotation"
# 1b is folded into the joined 1ab wherever the frameshift could be placed, so the
# checklist carries 1ab and falls back to a standalone 1b only when there is no join
GENES = ["1a", "1ab", "S", "3a", "3b", "3c", "M", "5a", "5b", "N"]
PRODUCT = {"1a": "1a polyprotein", "1ab": "1ab polyprotein", "1b": "1b polyprotein", "S": "spike glycoprotein",
           "3a": "3a protein", "3b": "3b protein", "3c": "envelope protein",
           "M": "membrane glycoprotein", "5a": "5a protein", "5b": "5b protein",
           "N": "nucleocapsid protein"}
BY_PRODUCT = {v: k for k, v in PRODUCT.items()}

def read_tbl(p):
    """(gene -> nt length, partial?) from a five-column feature table."""
    out, coords = {}, []
    lines = open(p).read().splitlines()
    for i, l in enumerate(lines):
        if not l.startswith("\t\t\tproduct\t"): continue
        prod = l.split("\t")[-1]
        g = BY_PRODUCT.get(prod)
        if g is None: continue
        segs, j = [], i - 1
        while j >= 0 and not lines[j].startswith("\t"):
            f = lines[j].split("\t")
            if len(f) >= 2 and f[0] and f[1]:
                segs.insert(0, (f[0], f[1]))
                if len(f) >= 3 and f[2] == "CDS": break
            j -= 1
        if not segs: continue
        nt = sum(abs(int(b.strip("<>")) - int(a.strip("<>"))) + 1 for a, b in segs)
        partial = any("<" in a or ">" in b for a, b in segs)
        out[g] = (nt, partial)
    return out

per = {}
for p in sorted(glob.glob(f"{TBL}/*.tbl")):
    g = read_tbl(p)
    if "1ab" not in g and "1b" in g: g["1ab"] = g.pop("1b")
    per[os.path.basename(p)[:-4]] = g

# the dataset consensus, taken from complete copies only
median = {}
for g in GENES:
    vals = [v[g][0] for v in per.values() if g in v and not v[g][1]]
    if vals: median[g] = int(statistics.median(vals))
print("=== dataset reference length per gene ===")
print("%-5s %8s %8s %8s %8s" % ("gene", "median", "min", "max", "n"))
for g in GENES:
    vals = sorted(v[g][0] for v in per.values() if g in v and not v[g][1])
    if vals: print("%-5s %8d %8d %8d %8d" % (g, median[g], vals[0], vals[-1], len(vals)))

rows, flags = [], collections.Counter()
for sid, genes in sorted(per.items()):
    miss, odd, partial = [], [], []
    for g in GENES:
        if g not in genes: miss.append(g); continue
        nt, part = genes[g]
        if part: partial.append(g); continue
        tol = max(0.10 * median[g], 45)
        if abs(nt - median[g]) > tol:
            odd.append("%s:%+d" % (g, nt - median[g]))
    status = ("complete" if not (miss or odd or partial) else
              "gene absent" if miss else "length outlier" if odd else "partial CDS")
    flags[status] += 1
    rows.append(dict(strain=sid, cds=len(genes), missing=";".join(miss),
                     length_outlier=";".join(odd), partial=";".join(partial), status=status))

with open(f"{OUT}/gene_qc.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

print("\nverdict:", " / ".join(f"{k} {v}" for k, v in flags.most_common()))
out = [r for r in rows if r["length_outlier"]]
print("\n=== length outliers: %d ===" % len(out))
for r in out:
    print("  %-14s %s" % (r["strain"], r["length_outlier"]))
c = collections.Counter(g for r in rows for x in r["missing"].split(";") if x for g in [x])
if c: print("\ngenes absent:", dict(c.most_common()))
