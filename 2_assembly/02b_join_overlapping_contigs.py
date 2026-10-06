#!/usr/bin/env python3
"""Join neighbouring contigs that share an exact terminal overlap, then check the
join against the reads.

Contigs are ordered by their main alignment block on the reference. A pair is
joined only if the end of one equals the start of the next, base for base, over
the assembler's final k (141 nt); one copy of the overlap is removed and nothing
else is changed. The reads are then mapped to the joined sequence: once the
duplicate is gone they map uniquely, so depth across the junction is a direct
test of whether the data support the join.
"""
import csv, glob, os, subprocess, statistics, collections
R = os.environ.get("IBV_WORK", os.getcwd())   # directory holding the reads and all intermediate output
S, FSA, WORK, OUT = R + "/29_validate_acc", R + "/09_annotation/fsa", R + "/31_fixable", R + "/32_merged"
BIN = os.environ.get("IBV_BIN", "")      # directory holding minimap2 and samtools; empty = on PATH
MM2, ST = os.path.join(BIN, "minimap2"), os.path.join(BIN, "samtools")
K = 141
os.makedirs(OUT, exist_ok=True)

def read_fasta(p):
    s, n = {}, None
    for ln in open(p):
        if ln.startswith(">"): n = ln[1:].split()[0]; s[n] = []
        else: s[n].append(ln.strip())
    return {k: "".join(v).upper() for k, v in s.items()}

def find_reads(s):
    for p in (s, s.replace("-", "_"), s.replace("_", "-")):
        for pat in ("/rawdata/*/%s_*R1*.fastq.gz", "/rawdata/*/%s*R1*.fastq.gz", "/rawdata/%s*R1*.fastq.gz"):
            h = sorted(glob.glob(R + pat % p))
            if h:
                r2 = h[0].replace("R1", "R2"); return h[0], (r2 if os.path.exists(r2) else "")
    return None, None

t1 = list(csv.DictReader(open(S + "/SupplementaryTable1_per_genome.tsv"), delimiter="\t"))
print("%-16s %7s %8s %8s  %-10s %-34s %s" % ("isolate", "contigs", "sum bp", "joined", "junction", "depth at junction (min / median)", "verdict"), flush=True)
for r in sorted([x for x in t1 if int(x["contigs"]) > 1], key=lambda x: -int(x["median_depth"])):
    s, racc = r["sample_name"], r["annotation_reference"].split("_")[-1]
    if s == "K2": continue
    fa = read_fasta(FSA + "/" + s + ".fsa")
    paf = subprocess.run([MM2, "-x", "asm20", "--secondary=no", WORK + "/" + racc + ".fa", FSA + "/" + s + ".fsa"],
                         capture_output=True, text=True).stdout
    best = {}
    for ln in paf.splitlines():
        f = ln.split("\t")
        if len(f) >= 9 and (f[0] not in best or int(f[8]) - int(f[7]) > best[f[0]][1]):
            best[f[0]] = (int(f[7]), int(f[8]) - int(f[7]))
    order = sorted(best, key=lambda c: best[c][0])
    if not all(fa[order[i]][-K:] == fa[order[i + 1]][:K] for i in range(len(order) - 1)):
        continue                      # not an exact-overlap case; handled separately
    seq, junc = fa[order[0]], []
    for c in order[1:]:
        junc.append(len(seq) - K // 2)          # centre of the shared k-mer
        seq += fa[c][K:]
    open(OUT + "/" + s + ".fa", "w").write(">" + s + "\n" + seq + "\n")
    r1, r2 = find_reads(s)
    bam = OUT + "/" + s + ".bam"
    subprocess.run("%s -ax sr -t 8 %s %s %s 2>/dev/null | %s sort -@2 -o %s - && %s index %s"
                   % (MM2, OUT + "/" + s + ".fa", r1, r2, ST, bam, ST, bam), shell=True)
    dep = [0] * (len(seq) + 2)
    for ln in subprocess.run([ST, "depth", "-a", "-Q", "20", bam], capture_output=True, text=True).stdout.splitlines():
        f = ln.split("\t"); dep[int(f[1])] = int(f[2])
    for j in junc:
        w = dep[j - 100:j + 100]
        ok = min(w) >= 10
        print("%-16s %7s %8d %8d  %-10d %-34s %s" % (s, r["contigs"], sum(len(fa[c]) for c in order), len(seq), j,
              "%d / %.0f" % (min(w), statistics.median(w)), "supported" if ok else "WEAK: under 10 reads"), flush=True)
print("done", flush=True)
