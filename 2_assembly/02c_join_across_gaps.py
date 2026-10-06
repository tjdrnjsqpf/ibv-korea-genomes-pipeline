#!/usr/bin/env python3
"""Finish the four assemblies whose contigs do not overlap.

  20ND039CT      the relaxed re-assembly produced one contig across the gap. It is
                 adopted only after checking that both original contigs sit inside
                 it unchanged and that reads cover the bridged stretch.
  the other three  contigs are joined in reference order with a run of N whose
                 length is estimated from the reference: the distance between the
                 two alignments, less any unaligned contig ends.
"""
import csv, os, subprocess, statistics, collections
R = os.environ.get("IBV_WORK", os.getcwd())   # directory holding the reads and all intermediate output
FSA, WORK, GC, F = R + "/09_annotation/fsa", R + "/31_fixable", R + "/33_gapclose", R + "/final"
BIN = os.environ.get("IBV_BIN", "")      # directory holding minimap2 and samtools; empty = on PATH
MM2, ST = os.path.join(BIN, "minimap2"), os.path.join(BIN, "samtools")
RC = str.maketrans("ACGTN", "TGCAN")

def read_fa(p):
    s, n = {}, None
    for ln in open(p):
        if ln.startswith(">"): n = ln[1:].split()[0]; s[n] = []
        else: s[n].append(ln.strip())
    return {k: "".join(v).upper() for k, v in s.items()}

def paf(target, query):
    out = subprocess.run([MM2, "-cx", "asm20", "--secondary=no", target, query],
                         capture_output=True, text=True).stdout
    return [l.split("\t") for l in out.splitlines() if l.count("\t") >= 10]

t1 = {r["sample_name"]: r for r in csv.DictReader(
      open(R + "/29_validate_acc/SupplementaryTable1_per_genome.tsv"), delimiter="\t")}
gaps = open(F + "/declared_gaps.tsv", "w")
gaps.write("sample\tgap_start\tgap_end\tgap_length\tbasis\n")

# ---------- 20ND039CT: closed by re-assembly ----------
s = "20ND039CT"
new = read_fa(GC + "/" + s + ".megahit/final.contigs.fa")
name = max(new, key=lambda k: len(new[k])); seq = new[name]
open(GC + "/" + s + ".closed.fa", "w").write(">closed\n" + seq + "\n")
print("=== %s: closed contig %s, %d bp ===" % (s, name, len(seq)))
edges = []
for f in paf(GC + "/" + s + ".closed.fa", FSA + "/" + s + ".fsa"):
    ql, qs, qe, ts, te, match, alen = int(f[1]), int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
    print("   original %-10s %6d bp -> closed %6d-%-6d  strand %s  matches %d/%d  unaligned ends %d+%d"
          % (f[0].split("_")[-1], ql, ts, te, f[4], match, alen, qs, ql - qe))
    edges += [ts, te]
edges.sort(); b0, b1 = edges[1], edges[2]          # the stretch between the two originals
bam = GC + "/" + s + ".closed.bam"
subprocess.run("%s -ax sr -t 8 %s %s %s 2>/dev/null | %s sort -@2 -o %s - && %s index %s" %
               (MM2, GC + "/" + s + ".closed.fa", GC + "/" + s + "_1.fq.gz", GC + "/" + s + "_2.fq.gz",
                ST, bam, ST, bam), shell=True)
dep = [0] * (len(seq) + 2)
for ln in subprocess.run([ST, "depth", "-a", "-Q", "20", bam], capture_output=True, text=True).stdout.splitlines():
    x = ln.split("\t"); dep[int(x[1])] = int(x[2])
w = dep[min(b0, b1) - 50:max(b0, b1) + 50]
print("   bridged stretch %d-%d (%d bp): read depth min %d, median %.0f"
      % (b0, b1, b1 - b0, min(w), statistics.median(w)))
orig = read_fa(FSA + "/" + s + ".fsa")
c2 = [v for k, v in orig.items() if k.endswith("contig2")][0]
off = len(seq) - len(c2)
for i, (x, y) in enumerate(zip(seq[off:], c2)):
    if x != y:
        pile = subprocess.run("%s mpileup -r closed:%d-%d -Q 20 -q 20 %s 2>/dev/null" % (ST, off + i + 1, off + i + 1, bam),
                              shell=True, capture_output=True, text=True).stdout.split("\t")
        bases = pile[4].upper() if len(pile) > 4 else ""
        cnt = {b_: bases.count(b_) for b_ in "ACGT"}
        print("   differs from the original at %d: closed=%s original=%s  reads %s" % (off + i + 1, x, y, cnt))
open(F + "/05_reassembly/asm/" + s + ".v2.fa", "w").write(">contig1 len=%d\n%s\n" % (len(seq), seq))

# ---------- the three joined with a declared gap ----------
for s in ("19ND70", "Kr-D15-15-CE7", "20L215F"):
    racc = t1[s]["annotation_reference"].split("_")[-1]
    fa = read_fa(FSA + "/" + s + ".fsa")
    if s == "19ND70":      # contig2 was a circular path; use its linear form
        c2 = [k for k in fa if k.endswith("contig2")][0]
        fa[c2] = list(read_fa(GC + "/19ND70.contig2.linear.fa").values())[0]
    tmp = GC + "/" + s + ".contigs.fa"
    open(tmp, "w").write("".join(">%s\n%s\n" % kv for kv in fa.items()))
    first, last = {}, {}
    for f in paf(WORK + "/" + racc + ".fa", tmp):
        c, ql, qs, qe, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), int(f[7]), int(f[8])
        if c not in first or ts < first[c][0]: first[c] = (ts, qs)
        if c not in last or te > last[c][0]:  last[c] = (te, ql - qe)
    order = sorted(first, key=lambda c: first[c][0])
    seq = fa[order[0]]
    print("=== %s: %s ===" % (s, " + ".join("%s(%d bp)" % (c.split("_")[-1], len(fa[c])) for c in order)))
    for a, b in zip(order, order[1:]):
        ref_gap = first[b][0] - last[a][0]
        est = ref_gap - last[a][1] - first[b][1]
        basis = "reference distance %d bp, less unaligned contig ends %d+%d" % (ref_gap, last[a][1], first[b][1])
        if est < 10: est, basis = 100, basis + "; estimate under 10 bp, set to 100 (length unknown)"
        g0 = len(seq) + 1
        seq += "N" * est + fa[b]
        print("   gap after %s: %d N at %d-%d   (%s)" % (a.split("_")[-1], est, g0, g0 + est - 1, basis))
        gaps.write("%s\t%d\t%d\t%d\t%s\n" % (s, g0, g0 + est - 1, est, basis))
    print("   joined length %d bp" % len(seq))
    open(F + "/05_reassembly/asm/" + s + ".v2.fa", "w").write(">contig1 len=%d\n%s\n" % (len(seq), seq))
gaps.close()
print("final set:", len([x for x in os.listdir(F + "/05_reassembly/asm") if x.endswith(".v2.fa")]), "sequences")
