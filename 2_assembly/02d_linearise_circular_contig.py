#!/usr/bin/env python3
"""19ND70 contig2 is a circular path: its two halves sit on the reference in the
opposite order to the contig, and its first and last 141 nt are the same. Cut the
circle where it does not correspond to the genome and write the linear order."""
import os, subprocess
R = os.environ.get("IBV_WORK", os.getcwd())   # directory holding the reads and all intermediate output
BIN = os.environ.get("IBV_BIN", "")      # directory holding minimap2; empty = on PATH
MM2 = os.path.join(BIN, "minimap2")
def read_fa(p):
    s, n = {}, None
    for ln in open(p):
        if ln.startswith(">"): n = ln[1:].split()[0]; s[n] = []
        else: s[n].append(ln.strip())
    return {k: "".join(v).upper() for k, v in s.items()}
fa = read_fa(R + "/09_annotation/fsa/19ND70.fsa")
c2 = [k for k in fa if k.endswith("contig2")][0]; q = fa[c2]
print("contig2 %d bp; first 141 == last 141: %s" % (len(q), q[:141] == q[-141:]))
# base-level alignment blocks on the reference
open("/tmp/_c2.fa", "w").write(">c2\n" + q + "\n")
out = subprocess.run([MM2, "-cx", "asm20", "--secondary=no", R + "/31_fixable/PQ387070.1.fa", "/tmp/_c2.fa"],
                     capture_output=True, text=True).stdout
blocks = sorted((int(f[7]), int(f[8]), int(f[2]), int(f[3]), int(f[9]), int(f[10]))
                for f in (l.split("\t") for l in out.splitlines()))
for b in blocks: print("   ref %d-%d  query %d-%d  identity %.1f%%" % (b[0], b[1], b[2], b[3], 100.0 * b[4] / b[5]))
(x_rs, x_re, x_qs, x_qe, _, _), (y_rs, y_re, y_qs, y_qe, _, _) = blocks
# linear order: the upstream half to the end of the contig, then, skipping the
# 141 nt the circle repeats, the downstream half to the end of its alignment
linear = q[x_qs:] + q[141:y_qe]
print("linear contig2: %d bp  (query %d-end, then 141-%d); dropped %d nt of junction between the halves"
      % (len(linear), x_qs, y_qe, x_qs - y_qe))
open(R + "/33_gapclose/19ND70.contig2.linear.fa", "w").write(">contig2_linear\n" + linear + "\n")
out = subprocess.run([MM2, "-cx", "asm20", "--secondary=no", R + "/31_fixable/PQ387070.1.fa",
                      R + "/33_gapclose/19ND70.contig2.linear.fa"], capture_output=True, text=True).stdout
for l in out.splitlines():
    f = l.split("\t"); print("   after: ref %s-%s  query %s-%s of %s  identity %.1f%%" % (f[7], f[8], f[2], f[3], f[1], 100.0 * int(f[9]) / int(f[10])))
