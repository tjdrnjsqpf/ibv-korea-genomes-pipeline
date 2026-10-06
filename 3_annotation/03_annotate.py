#!/usr/bin/env python
"""Annotate each rebuilt genome against its own lineage reference, then judge the result.

Three things the previous pass got wrong, all of them in the checking rather than the data:

  1a was reported missing in 28 genomes because several panel references annotate only a
  joined 1ab feature. 1a and 1b are now taken from that join's first and second segments
  whenever they are not annotated separately.

  Start codons were judged at exactly the projected position, so a single indel upstream
  moved the frame and 127 genomes came back as codon errors. The start is now sought in a
  window around the projection, preferring the ATG whose ORF length lands closest to the
  reference.

  A flat 10% length tolerance is far too tight on the short ORFs - 3a is 174 bp, so six
  codons of ordinary strain variation tripped it. Tolerance is now whichever is larger,
  10% or 45 bp, and a real premature stop is reported separately.
"""
import sys, os, re, glob, csv, subprocess, collections

base = sys.argv[1]
ASM   = f"{base}/05_reassembly/asm"
REFGB = f"{base}/09_annotation/refgb"
PANEL = f"{base}/05_reassembly/ref/panel"
OUT   = f"{base}/09_annotation"
RUN   = os.environ.get("IBV_RUN", "").split()   # optional prefix, e.g. "micromamba run -n ibv_asm"
os.makedirs(f"{OUT}/tbl", exist_ok=True)
os.makedirs(f"{OUT}/fsa", exist_ok=True)
STOPS = {"TAA", "TAG", "TGA"}
SLIPPERY = "TTTAAAC"   # the -1 frameshift site of the replicase
START_WINDOW = 30      # nt either side of the projected start to look for the real ATG

GENES = ["1a", "1b", "S", "3a", "3b", "3c", "M", "5a", "5b", "N"]
PATTERN = {"1a": r"^(orf ?)?1a( polyprotein)?$|^polyprotein 1a$",
           "1b": r"^(orf ?)?1b( polyprotein)?$|^polyprotein 1b$",
           "S":  r"spike", "3a": r"^(orf ?)?3a( protein)?$", "3b": r"^(orf ?)?3b( protein)?$",
           "3c": r"^(orf ?)?3c( protein)?$|envelope", "M": r"membrane",
           "5a": r"^(orf ?)?5a( protein)?$", "5b": r"^(orf ?)?5b( protein)?$",
           "N":  r"nucleocapsid|^n protein$"}
PRODUCT = {"1a": "1a polyprotein", "1b": "1b polyprotein", "S": "spike glycoprotein",
           "3a": "3a protein", "3b": "3b protein", "3c": "envelope protein",
           "M": "membrane glycoprotein", "5a": "5a protein", "5b": "5b protein",
           "N": "nucleocapsid protein"}

def read_fa(p):
    d, k = {}, None
    for line in open(p):
        if line.startswith(">"): k = line[1:].split()[0]; d[k] = []
        else: d[k].append(line.strip().upper())
    return {k: "".join(v) for k, v in d.items()}

def parse_cds(p):
    feats, cur = [], None
    for line in open(p):
        if re.match(r"^     CDS  ", line):
            if cur: feats.append(cur)
            cur = {"loc": line[21:].strip(), "product": None, "gene": None}
        elif cur is not None and line.startswith(" " * 21):
            t = line.strip()
            # some records label the polyprotein only with /gene, so fall back to it
            if t.startswith("/product="): cur["product"] = t.split("=", 1)[1].strip('"')
            elif t.startswith("/gene=") and not cur.get("gene"):
                cur["gene"] = t.split("=", 1)[1].strip('"')
            elif cur["product"] is None and not t.startswith("/"): cur["loc"] += t
        elif cur is not None and re.match(r"^     \S", line):
            feats.append(cur); cur = None
    if cur: feats.append(cur)
    out = []
    for f in feats:
        seg = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", f["loc"])]
        name = (f["product"] or f.get("gene") or "").strip().lower()
        if seg: out.append(dict(product=name, seg=seg))
    return out

def canonical(gbfile):
    genes, feats = {}, parse_cds(gbfile)
    for f in feats:
        for g in GENES:
            if g in genes: continue
            if re.search(PATTERN[g], f["product"]):
                genes[g] = dict(start=f["seg"][0][0],
                                expected=sum(b - a + 1 for a, b in f["seg"]))
                break
    # references that annotate only the joined 1ab still define 1a and 1b through it
    joined = next((f for f in feats if len(f["seg"]) > 1 and re.search(r"1ab", f["product"])), None)
    if joined:
        if "1a" not in genes:
            a, b = joined["seg"][0]
            genes["1a"] = dict(start=a, expected=b - a + 1)
        if "1b" not in genes:
            a = joined["seg"][1][0]
            genes["1b"] = dict(start=a, expected=sum(y - x + 1 for x, y in joined["seg"][1:]))
    # one panel record (PZ330336.1) labels the 1a reading frame "ORF1ab" as a single
    # interval and gives 1b separately. A one-segment "1ab" of 1a length is 1a.
    if "1a" not in genes:
        mis = next((f for f in feats if len(f["seg"]) == 1 and re.search(r"1ab", f["product"])
                    and f["seg"][0][1] - f["seg"][0][0] + 1 < 13000), None)
        if mis:
            a, b = mis["seg"][0]
            genes["1a"] = dict(start=a, expected=b - a + 1)
    return genes

def rc(s): return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def align(q, r):
    open("/tmp/_r.fa", "w").write(">ref\n%s\n" % r)
    open("/tmp/_q.fa", "w").write(">q\n%s\n" % q)
    return subprocess.run(RUN + ["minimap2", "-c", "--cs", "-x", "asm20", "--secondary=no",
                                 "/tmp/_r.fa", "/tmp/_q.fa"], capture_output=True, text=True).stdout

def project(qseq, refseq):
    paf = align(qseq, refseq)
    plus  = sum(int(l.split("\t")[9]) for l in paf.splitlines() if l.split("\t")[4] == "+")
    minus = sum(int(l.split("\t")[9]) for l in paf.splitlines() if l.split("\t")[4] == "-")
    if minus > plus:
        qseq = rc(qseq); paf = align(qseq, refseq)
    r2q = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if f[4] != "+": continue
        cs = next((x[5:] for x in f if x.startswith("cs:Z:")), "")
        if not cs: continue
        qp, rp = int(f[2]), int(f[7])
        for op, val in re.findall(r"([:=*+-])([0-9A-Za-z]+)", cs):
            if op in ":=":
                n = int(val) if op == ":" else len(val)
                for i in range(n): r2q[rp + i] = qp + i
                qp += n; rp += n
            elif op == "*": r2q[rp] = qp; qp += 1; rp += 1
            elif op == "+": qp += len(val)
            elif op == "-": rp += len(val)
    return r2q, qseq

def orf_from(seq, s):
    i = s
    while i + 3 <= len(seq) and seq[i:i+3] not in STOPS: i += 3
    return i + 3, (i + 3 <= len(seq))

def best_start(seq, s0, expected, require_atg=True):
    """The projected position can be off by an indel; take the nearby ATG that yields the
    ORF closest to the reference length.

    1b is the exception: it has no start codon of its own, being reached from 1a through a
    -1 ribosomal frameshift, so there the projected position is used as given."""
    if not require_atg:
        end, complete = orf_from(seq, s0)
        return (abs((end - s0) - expected), s0, end, complete)
    best = None
    for off in sorted(range(-START_WINDOW, START_WINDOW + 1), key=abs):
        s = s0 + off
        if s < 0 or s + 3 > len(seq) or seq[s:s+3] != "ATG": continue
        end, complete = orf_from(seq, s)
        score = abs((end - s) - expected)
        if best is None or score < best[0]: best = (score, s, end, complete)
        if score == 0: break
    return best


def slip_site(seq, a_start, b):
    """Find where the -1 frameshift actually sits in this genome, and where 1b then ends.

    Projecting the reference join boundary is not reliable: an indel upstream moves it, and
    a boundary landing after 1a's own stop codon drags that stop into the joined product -
    which is what table2asn reported as internal stops. The boundary is therefore chosen on
    the sequence itself. A candidate has to keep the first segment in frame and stop-free,
    and reading on from it has to give a second segment of about 1b's expected length. The
    second segment's frame is set by the candidate, not by the projected 1b start, since
    that projection is exactly what may be off. Ties go to the candidate nearest it.
    """
    proj, exp = b["start"], b["expected"]
    lo, hi = max(a_start + 3, proj - 300), min(proj + 300, len(seq) - 3)
    best = None
    for x in range(lo, hi + 1):
        if (x - a_start + 1) % 3: continue
        seg1 = seq[a_start - 1:x]
        # every codon of the first segment is translated inside the join, including the
        # last one - so 1a's own stop codon must not fall anywhere in it
        if any(seg1[i:i+3] in STOPS for i in range(0, len(seg1) - 2, 3)): continue
        end, complete = orf_from(seq, x - 1)
        if not complete: continue
        n2 = end - (x - 1)
        if abs(n2 - exp) > max(0.10 * exp, 60): continue
        # the frameshift happens on the slippery sequence UUUAAAC, and the join belongs on
        # its last base. Several panel references annotate 1b from a point 27 nt further
        # on, which also reads through without a stop but mistranslates the nine codons in
        # between - so a candidate on the heptamer always wins over one that is merely
        # nearer the projection.
        on_heptamer = seq[x - 7:x] == SLIPPERY
        d = (0 if on_heptamer else 1, abs(x - proj))
        if best is None or d < best[0]: best = (d, x, end)
    if best is None: return proj, b["end"], b["complete"]
    return best[1], best[2], True

def partial_replicase(seq, b):
    """A genome that begins inside 1a has no 1a start to anchor the join. If 1b begins
    on the slippery sequence and the frame upstream of it is open back to the first
    bases, the replicase can still be given as 5' partial features.
    Returns (codon_start, join position, end of 1a) or None."""
    x = b["start"]
    if seq[x - 7:x] != SLIPPERY: return None
    for cs in (1, 2, 3):
        if (x - cs + 1) % 3: continue
        seg1 = seq[cs - 1:x]
        if any(seg1[i:i+3] in STOPS for i in range(0, len(seg1) - 2, 3)): return None
        a_end, complete = orf_from(seq, cs - 1)
        if not complete or a_end - x > 90: return None
        return cs, x, a_end
    return None

refchoice = dict(l.split("\t")[:2] for l in open(f"{base}/05_reassembly/ref_choice.tsv"))
cache = {}
def get_ref(name):
    if name not in cache:
        acc = name.split("_", 1)[1]
        cache[name] = (list(read_fa(f"{PANEL}/{name}.fasta").values())[0],
                       canonical(f"{REFGB}/{acc}.gb"))
    return cache[name]

rows = []
for p in sorted(glob.glob(f"{ASM}/*.v2.fa")):
    sid = os.path.basename(p)[:-len(".v2.fa")]
    rname = refchoice.get(sid)
    if not rname: continue
    refseq, genes = get_ref(rname)
    contigs, found = read_fa(p), {}
    for cname in list(contigs):
        r2q, oriented = project(contigs[cname], refseq)
        contigs[cname] = oriented
        if not r2q: continue
        for g in GENES:
            if g in found or g not in genes: continue
            exp = genes[g]["expected"]
            s0 = r2q.get(genes[g]["start"] - 1)
            if s0 is None:
                # the reference start can fall in a stretch that does not align (the
                # first codons of M vary between lineages). Extrapolate from the first
                # aligned position inside the gene; the ATG search below then settles it.
                st = genes[g]["start"] - 1
                k = next((p_ for p_ in range(st, st + min(300, exp)) if p_ in r2q), None)
                if k is None: continue
                s0 = r2q[k] - (k - st)
                if s0 < 0: continue
            hit = best_start(oriented, s0, exp, require_atg=(g != "1b"))
            if hit is None: continue
            _, s, end, complete = hit
            # an ORF that hits a stop within a few codons is a projection artefact, not a gene
            if end - s < 60: continue
            found[g] = dict(contig=cname, start=s + 1, end=end, nt=end - s,
                            expected=exp, complete=complete)

    miss, offlen, premature, runoff = [], [], [], []
    for g in GENES:
        v = found.get(g)
        if v is None: miss.append(g); continue
        tol = max(0.10 * v["expected"], 45)
        if not v["complete"]: runoff.append(g)
        elif v["nt"] < v["expected"] - tol: premature.append(g)
        elif abs(v["nt"] - v["expected"]) > tol: offlen.append(g)
    status = ("complete" if not (miss or offlen or premature or runoff) else
              "gene absent" if miss else "premature stop" if premature else
              "at contig end" if runoff else "length variant")
    # coordinates below refer to the oriented sequence, so that sequence - not the original
    # contig - is what gets written out and validated
    names = {c: (sid if len(contigs) == 1 else f"{sid}_{c}") for c in contigs}
    with open(f"{OUT}/fsa/{sid}.fsa", "w") as fh:
        for cname, seq in contigs.items():
            fh.write(">%s [organism=Infectious bronchitis virus] [isolate=%s] [gcode=1]\n"
                     % (names[cname], sid))
            for i in range(0, len(seq), 70): fh.write(seq[i:i+70] + "\n")
    with open(f"{OUT}/tbl/{sid}.tbl", "w") as fh:
        for cname, seq in contigs.items():
            fh.write(">Feature %s\n" % names[cname])
            a, b = found.get("1a"), found.get("1b")
            # GenBank annotates the replicase as a single 1ab joined across the -1 frameshift,
            # with 1a as its own CDS. Emitting a bare 1b instead leaves a CDS with no start
            # codon, which is what table2asn was objecting to.
            joined = (a and b and a["contig"] == cname and b["contig"] == cname)
            part = None
            if not a and b and b["contig"] == cname:
                part = partial_replicase(seq, b)
                if part:
                    cs, x, a_end = part
                    bend = min(b["end"], len(seq))
                    ben = "%d" % bend if b["complete"] and bend == b["end"] else ">%d" % bend
                    fh.write("<1\t%d\tCDS\n%d\t%s\n\t\t\tproduct\t1ab polyprotein\n"
                             "\t\t\tribosomal_slippage\n\t\t\tcodon_start\t%d\n" % (x, x, ben, cs))
                    fh.write("<1\t%d\tCDS\n\t\t\tproduct\t1a polyprotein\n"
                             "\t\t\tcodon_start\t%d\n" % (a_end, cs))
            for g in GENES:
                v = found.get(g)
                if not v or v["contig"] != cname: continue
                if g == "1b" and (joined or part): continue
                end = min(v["end"], len(seq))
                st = "%d" % v["start"] if seq[v["start"]-1:v["start"]+2] == "ATG" else "<%d" % v["start"]
                en = "%d" % end if v["complete"] and end == v["end"] else ">%d" % end
                if g == "1a" and joined:
                    slip, b_end, b_ok = slip_site(seq, v["start"], b)
                    bend = min(b_end, len(seq))
                    ben = "%d" % bend if b_ok and bend == b_end else ">%d" % bend
                    fh.write("%s\t%d\tCDS\n%d\t%s\n\t\t\tproduct\t1ab polyprotein\n"
                             "\t\t\tribosomal_slippage\n" % (st, slip, slip, ben))
                fh.write("%s\t%s\tCDS\n\t\t\tproduct\t%s\n" % (st, en, PRODUCT[g]))
    rows.append(dict(strain=sid, reference=rname, contigs=len(contigs), cds=len(found),
                     missing=";".join(miss), premature_stop=";".join(premature),
                     length_variant=";".join(offlen), at_contig_end=";".join(runoff),
                     status=status))

with open(f"{OUT}/annotation_qc.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print("%d genomes" % len(rows))
print("verdict:", " / ".join(f"{k} {v}" for k, v in collections.Counter(r["status"] for r in rows).most_common()))
print("CDS annotated:", dict(sorted(collections.Counter(r["cds"] for r in rows).items())))
for lbl, key in (("gene absent", "missing"), ("premature stop", "premature_stop"),
                 ("length variant", "length_variant"), ("at contig end", "at_contig_end")):
    c = collections.Counter(g for r in rows for g in r[key].split(";") if g)
    if c: print(f"{lbl}:", dict(c.most_common()))
