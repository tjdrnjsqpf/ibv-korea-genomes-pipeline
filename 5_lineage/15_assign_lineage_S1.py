#!/usr/bin/env python3
"""Assign a lineage and a sub-lineage to each deposited genome from the S1 tree.

For each of our genomes, walk from its tip towards the root. The first node whose
descendant reference tips all carry the same label, and which is supported
(SH-aLRT >= 80 and UFBoot >= 95), gives that label. Lineage comes from the
Valastro prototypes, sub-lineage from the labelled Korean context. A genome that
reaches the root without such a node is reported as unresolved rather than forced.
"""
import csv, os, re, sys

D = os.path.join(os.environ.get("IBV_WORK", os.getcwd()), "41_phylo_allproto")
MIN_ALRT, MIN_UFB = 80.0, 95.0

# ---------- newick ----------
class Node:
    __slots__ = ("children", "parent", "name", "alrt", "ufb")
    def __init__(self):
        self.children = []; self.parent = None
        self.name = None; self.alrt = None; self.ufb = None

def parse(nwk):
    nwk = nwk.strip().rstrip(";")
    root = Node(); cur = root; tok = ""
    def flush(n, t):
        t = t.strip()
        if not t: return
        lab = t.split(":")[0]
        if n.children:                       # internal: label is support
            m = re.match(r"^([\d.]+)/([\d.]+)$", lab)
            if m: n.alrt, n.ufb = float(m.group(1)), float(m.group(2))
        elif lab:
            n.name = lab
    for ch in nwk:
        if ch == "(":
            nn = Node(); nn.parent = cur; cur.children.append(nn); cur = nn
        elif ch == ",":
            flush(cur, tok); tok = ""
            nn = Node(); nn.parent = cur.parent; cur.parent.children.append(nn); cur = nn
        elif ch == ")":
            flush(cur, tok); tok = ""
            cur = cur.parent
        else:
            tok += ch
    flush(cur, tok)
    return root

def tips(n, acc=None):
    acc = [] if acc is None else acc
    if not n.children: acc.append(n)
    for c in n.children: tips(c, acc)
    return acc

# ---------- labels ----------
proto = {r["taxon"]: r["lineage"] for r in
         csv.DictReader(open(D + "/00_data/valastro_labels.tsv"), delimiter="\t")}
# four genotypes described after Valastro 2016 (GVII-1, GVIII, GIX-1, GX-1)
proto.update({r["taxon"]: r["lineage"] for r in
              csv.DictReader(open(D + "/00_data/new_proto_labels.tsv"), delimiter="\t")})
sub = {r["taxon"]: r["lineage"] for r in
       csv.DictReader(open(D + "/00_data/s1_lineage_labels.tsv"), delimiter="\t")}
ours = {l[1:].strip() for l in open(D + "/00_data/ours_S1.fa") if l.startswith(">")}

def sub_of(k):
    return sub.get(k) or sub.get(k.split(".")[0])

root = parse(open(D + "/02_tree/s1.treefile").read())
all_tips = tips(root)
byname = {t.name: t for t in all_tips if t.name}
print("tips in tree: %d, ours: %d" % (len(all_tips), sum(1 for t in all_tips if t.name in ours)),
      file=sys.stderr)

def climb(tip, labeller):
    """nearest supported ancestor whose reference descendants agree on a label"""
    n = tip
    while n.parent is not None:
        n = n.parent
        labs = {labeller(t.name) for t in tips(n) if t.name and t.name not in ours}
        labs.discard(None)
        if len(labs) == 1:
            ok = (n.alrt is None or n.alrt >= MIN_ALRT) and (n.ufb is None or n.ufb >= MIN_UFB)
            sup = "%.0f/%.0f" % (n.alrt or 0, n.ufb or 0)
            if ok:
                return labs.pop(), sup, len(tips(n))
        elif len(labs) > 1:
            return None, "", 0
    return None, "", 0

rows = []
for s in sorted(ours):
    t = byname.get(s)
    if t is None:
        rows.append((s, "no S1", "", "", "", "")); continue
    lin, lsup, _ = climb(t, lambda k: proto.get(k))
    sl,  ssup, _ = climb(t, sub_of)
    rows.append((s, lin or "unresolved", lsup, sl or "unresolved", ssup, ""))

with open(D + "/lineage_from_tree.tsv", "w") as f:
    f.write("sample\tlineage\tlineage_support\tsublineage\tsublineage_support\tnote\n")
    for r in rows: f.write("\t".join(r) + "\n")

import collections
for col, nm in ((1, "lineage"), (3, "sub-lineage")):
    print("\n=== %s ===" % nm)
    for k, v in collections.Counter(r[col] for r in rows).most_common():
        print("   %-16s %3d" % (k, v))
