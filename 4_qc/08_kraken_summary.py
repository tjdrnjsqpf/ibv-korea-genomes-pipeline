#!/usr/bin/env python
"""Per-sample read composition: how much of each library is actually IBV.

These libraries came from clinical tissue, so the IBV fraction and what the rest of the
reads are belong in the technical validation - and a low IBV fraction flags a genome that
was assembled from thin evidence.
"""
import sys, os, glob, csv
base = sys.argv[1]
rows = []
for p in sorted(glob.glob(f"{base}/08_kraken/reports/*.report")):
    sid = os.path.basename(p)[:-len(".report")]
    tot = unclass = ibv = corona = host = bact = 0
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6: continue
        pct, clade, direct, rank, taxid, name = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5].strip()
        if rank == "U": unclass = clade
        if rank == "R" and name == "root": tot = clade + unclass
        if name == "Infectious bronchitis virus": ibv = clade
        if name in ("Gammacoronavirus", "Orthocoronavirinae", "Coronaviridae"): corona = max(corona, clade)
        if name in ("Gallus gallus", "Gallus"): host = max(host, clade)
        if name == "Bacteria": bact = clade
    total = tot or (unclass + ibv)
    pctof = lambda x: round(100 * x / total, 2) if total else 0.0
    rows.append(dict(strain=sid, total_reads=total, IBV=ibv, pct_IBV=pctof(ibv),
                     pct_corona=pctof(corona), pct_host=pctof(host),
                     pct_bacteria=pctof(bact), pct_unclassified=pctof(unclass)))
rows.sort(key=lambda r: r["pct_IBV"])
with open(f"{base}/08_kraken/kraken_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
v = sorted(r["pct_IBV"] for r in rows)
print("%d strains · IBV read fraction  min %.1f%% / median %.1f%% / max %.1f%%"
      % (len(rows), v[0], v[len(v)//2], v[-1]))
print("below 5%% IBV: %d" % sum(1 for x in v if x < 5))
print("\n=== 12 lowest IBV fractions ===")
print("%-14s %11s %8s %8s %8s %8s" % ("strain","reads","IBV%","host%","bacteria%","unclassified%"))
for r in rows[:12]:
    print("%-14s %11d %8.2f %8.2f %8.2f %8.2f" % (r["strain"], r["total_reads"],
          r["pct_IBV"], r["pct_host"], r["pct_bacteria"], r["pct_unclassified"]))
