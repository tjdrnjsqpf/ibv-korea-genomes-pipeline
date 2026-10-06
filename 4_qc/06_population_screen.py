#!/usr/bin/env python
"""Classify every genome by how its within-sample variation is distributed.

Three shapes come out of the minor-allele spectrum. A two-strain mixture pushes the
differing sites to one shared fraction well away from zero and spreads them along the whole
genome. Sequencing noise and quasispecies pile up just above the 5% floor and decay. A
clonal sample has almost no qualifying sites at all.
"""
import sys, os, glob, subprocess, collections, csv
base = sys.argv[1]
RUN = os.environ.get("IBV_RUN", "").split()   # optional prefix, e.g. "micromamba run -n ibv_asm"

def profile(sid):
    ref = f"{base}/05_reassembly/asm/{sid}.v2.fa"
    bam = f"{base}/05_reassembly/consensus/{sid}.v2.bam"
    mp = subprocess.run(RUN + ["bcftools", "mpileup", "-f", ref, "-d", "5000",
                               "-a", "AD", "-q", "20", "-Q", "20", "-Ou", bam], capture_output=True)
    cl = subprocess.run(RUN + ["bcftools", "call", "-m", "-A", "-Ov"],
                        input=mp.stdout, capture_output=True)
    fr, pos, L = [], [], 0
    for line in cl.stdout.decode(errors="replace").splitlines():
        if line.startswith("##contig=<ID="): L += int(line.split("length=")[1].rstrip(">"))
        if line.startswith("#"): continue
        f = line.split("\t")
        if f[4] in (".", "<*>"): continue
        try: ad = f[9].split(":")[f[8].split(":").index("AD")]
        except (ValueError, IndexError): continue
        c = [int(x) for x in ad.split(",") if x.isdigit()]
        dp = sum(c)
        if dp < 50 or len(c) < 2: continue
        a = max(c[1:]) / dp
        if 0.05 <= a <= 0.95: fr.append(a); pos.append((f[0], int(f[1])))
    if not fr: return dict(n=0, peak=0, spread=0.0, low_frac=0.0)
    h = collections.Counter(min(int(a * 20), 19) for a in fr)
    win = {(c, p // 1000) for c, p in pos}
    return dict(n=len(fr), peak=max(h, key=lambda k: h[k]) * 5,
                spread=round(100 * len(win) / max(1, L // 1000), 1),
                low_frac=round(100 * sum(1 for a in fr if a < 0.15) / len(fr), 1))

rows = []
for p in sorted(glob.glob(f"{base}/05_reassembly/asm/*.v2.fa")):
    sid = os.path.basename(p)[:-len(".v2.fa")]
    if not os.path.exists(f"{base}/05_reassembly/consensus/{sid}.v2.bam"): continue
    r = profile(sid); r["strain"] = sid
    r["call"] = ("mixed" if r["n"] >= 200 and r["peak"] >= 25 and r["low_frac"] < 45
                 else "low-level diversity" if r["n"] >= 200 else "clonal")
    rows.append(r)
rows.sort(key=lambda r: -r["n"])
with open(f"{base}/05_reassembly/mixclass.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["strain", "n", "peak", "spread", "low_frac", "call"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(" / ".join(f"{k} {v}" for k, v in collections.Counter(r["call"] for r in rows).items()))
print("\n%-14s %8s %9s %8s %9s  %s" % ("strain", "sites", "peak", "spread%", "low%", "verdict"))
for r in rows[:18]:
    print("%-14s %8d %8d%% %7.0f%% %8.0f%%  %s" % (r["strain"], r["n"], r["peak"],
                                                   r["spread"], r["low_frac"], r["call"]))
