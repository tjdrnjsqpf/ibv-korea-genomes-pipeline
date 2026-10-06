#!/usr/bin/env python
"""Assemble everything NCBI will ask for into one directory.

Three submissions have to agree with each other - BioProject, then BioSample per isolate,
then SRA for the reads and GenBank for the assemblies - and they are linked by sample name,
so every table here is keyed on the GenBank Sequence_ID already in use. The per-genome QC
columns travel with it, since a data descriptor has to state what each record's evidence is.
"""
import sys, os, glob, csv, gzip, hashlib, statistics, collections

base = sys.argv[1]
OUT  = f"{base}/11_submission"
for d in ("genomes", "annotation", "metadata"): os.makedirs(f"{OUT}/{d}", exist_ok=True)

def load(p, key=None):
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    return {r[key]: r for r in rows} if key else rows

bios  = load(f"{base}/05_reassembly/biosample_attributes.tsv", "sample_name")
samp  = [l.split("\t") for l in open(f"{base}/05_reassembly/samples_all.tsv").read().splitlines() if l.strip()]
geneq = load(f"{base}/09_annotation/gene_qc.tsv", "strain")
mixc  = load(f"{base}/05_reassembly/mixclass.tsv", "strain")
krak  = load(f"{base}/08_kraken/kraken_summary.tsv", "strain")
refch = dict(l.split("\t")[:2] for l in open(f"{base}/05_reassembly/ref_choice.tsv"))

def md5(p, chunk=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""): h.update(b)
    return h.hexdigest()

manifest, sra = [], []
for sid, prefix, batch in samp:
    fsa = f"{base}/09_annotation/fsa/{sid}.fsa"
    tbl = f"{base}/09_annotation/tbl/{sid}.tbl"
    if not os.path.exists(fsa): continue
    for src, dst in ((fsa, f"{OUT}/genomes/{sid}.fsa"), (tbl, f"{OUT}/annotation/{sid}.tbl")):
        if os.path.exists(src):
            open(dst, "w").write(open(src).read())

    seqs, k = {}, None
    for line in open(fsa):
        if line.startswith(">"): k = line[1:].split()[0]; seqs[k] = []
        else: seqs[k].append(line.strip())
    total = sum(len(("".join(v))) for v in seqs.values())
    ns    = sum("".join(v).count("N") for v in seqs.values())

    dp = []
    dpath = f"{base}/05_reassembly/depth2/{sid}.perbase.gz"
    if os.path.exists(dpath):
        dp = [int(l.rstrip("\n").split("\t")[2]) for l in gzip.open(dpath, "rt")]

    r1 = sorted(glob.glob(f"{base}/rawdata/**/{prefix}_S*_L001_R1_001.fastq*", recursive=True))
    r2 = sorted(glob.glob(f"{base}/rawdata/**/{prefix}_S*_L001_R2_001.fastq*", recursive=True))
    r1 = [f for f in r1 if "_not_in_study" not in f]
    r2 = [f for f in r2 if "_not_in_study" not in f]
    if r1 and r2:
        sra.append(dict(sample_name=sid, library_ID=f"{sid}_lib",
                        title=f"Illumina MiniSeq sequencing of infectious bronchitis virus {sid}",
                        library_strategy="WGS", library_source="VIRAL RNA",
                        library_selection="RANDOM", library_layout="paired",
                        platform="ILLUMINA", instrument_model="Illumina MiniSeq",
                        design_description=("Total RNA from clinical tissue, random-primed cDNA, "
                                            "Illumina library, paired-end sequencing"),
                        filetype="fastq", filename=os.path.basename(r1[0]),
                        filename2=os.path.basename(r2[0]), sequencing_batch=batch))

    g, m, kr = geneq.get(sid, {}), mixc.get(sid, {}), krak.get(sid, {})
    manifest.append(dict(
        sample_name=sid, contigs=len(seqs), total_bp=total, ambiguous_bp=ns,
        median_depth=int(statistics.median(dp)) if dp else "",
        zero_cov_bp=sum(1 for x in dp if x == 0) if dp else "",
        pct_IBV_reads=kr.get("pct_IBV", ""), lineage_reference=refch.get(sid, ""),
        cds_annotated=g.get("cds", ""), annotation_status=g.get("status", ""),
        cds_missing=g.get("missing", ""), cds_length_outlier=g.get("length_outlier", ""),
        population=m.get("call", ""), sequencing_batch=batch,
        assembly_md5=md5(f"{OUT}/genomes/{sid}.fsa")))

def write(path, rows):
    if not rows: return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("  %-34s %4d rows" % (os.path.basename(path), len(rows)))

print("written:")
write(f"{OUT}/metadata/biosample_attributes.tsv", [bios[s] for s, _, _ in samp if s in bios])
write(f"{OUT}/metadata/sra_metadata.tsv", sra)
write(f"{OUT}/MANIFEST.tsv", manifest)

print("\nBioSample %d / SRA %d / genomes %d" % (len(bios), len(sra), len(manifest)))
print("annotation verdict:", " / ".join(f"{k} {v}" for k, v in
      collections.Counter(r["annotation_status"] for r in manifest).most_common()))
miss = [r["sample_name"] for r in manifest if not r["pct_IBV_reads"]]
if miss: print("no kraken result:", miss[:5])
