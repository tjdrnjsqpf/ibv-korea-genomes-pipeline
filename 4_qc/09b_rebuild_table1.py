#!/usr/bin/env python3
"""Rebuild Supplementary Table 1 for the final 182 genomes.

Every assembly-derived column is recomputed from the final sequences; columns
that describe the library or the sample are carried over. Values that were in
Korean in the working table are written in English.
"""
import csv, gzip, hashlib, json, os, statistics
R = os.environ.get("IBV_WORK", os.getcwd())   # directory holding the reads and all intermediate output
F = R + "/final"
STATUS = {"complete": "complete", "gene absent": "gene absent", "length outlier": "length outlier",
          "partial CDS": "partial CDS", "length variant": "length outlier", "premature stop": "premature stop",
          "at contig end": "partial CDS"}
POP = {"clonal": "clonal", "low-level diversity": "low-level diversity", "mixed": "mixed"}
MERGED141 = set("20ND27 20ND42 21ND079CT 20ND30 22ND370CT 20ND71 Kr-Q35-06-CE 20L124CL 19ND16 "
                "21L034CL 20AD042TraCE1 15D44CT-CE2 22ND259CT Kr-342-03".split())
old = list(csv.DictReader(open(R + "/29_validate_acc/SupplementaryTable1_per_genome.tsv"), delimiter="\t"))
gq = {r["strain"]: r for r in csv.DictReader(open(F + "/09_annotation/gene_qc.tsv"), delimiter="\t")}
gaps = {}
for r in csv.DictReader(open(F + "/declared_gaps.tsv"), delimiter="\t"):
    gaps[r["sample"]] = r
cols = list(old[0].keys())
i = cols.index("ambiguous_bp")
cols[i:i+1] = ["ambiguous_bp", "assembly_route", "declared_gap"]
# lineage and sub-lineage from the S1 tree of the final genomes, next to the prior label
tree = {r["sample"]: r for r in csv.DictReader(open(R + "/27_phylo_final/lineage_from_tree.tsv"), delimiter="\t")}
j = cols.index("prior_sublineage")
cols[j:j+1] = ["lineage_S1", "lineage_S1_support", "sublineage_S1", "sublineage_S1_support", "prior_sublineage"]
out = []
for r in old:
    s = r["sample_name"]
    seq = "".join(l.strip() for l in open(F + "/09_annotation/fsa/" + s + ".fsa") if not l.startswith(">")).upper()
    dp = [int(l.split("\t")[2]) for l in gzip.open(F + "/05_reassembly/depth2/" + s + ".perbase.gz", "rt")]
    assert len(dp) == len(seq), (s, len(dp), len(seq))
    g = gq[s]
    if s in MERGED141:        route = "contigs joined at exact 141-nt overlap"
    elif s == "20ND039CT":    route = "contigs joined at 2-nt overlap, confirmed by re-assembly"
    elif s in gaps:           route = "contigs joined across declared gap"
    else:                     route = "single contig"
    if s == "K2":             route = "single contig, assembled from two pooled libraries"
    gp = gaps.get(s)
    r = dict(r)
    r.update(contigs=1, total_bp=len(seq), ambiguous_bp=seq.count("N"), assembly_route=route,
             declared_gap=("%s-%s (%s)" % (gp["gap_start"], gp["gap_end"],
                           "length unknown" if "unknown" in gp["basis"] else "estimated %s nt" % gp["gap_length"])) if gp else "",
             median_depth=int(statistics.median(dp)),
             zero_cov_bp=sum(1 for i, x in enumerate(dp) if x == 0 and seq[i] != "N"),
             cds_annotated=g["cds"], annotation_status=STATUS[g["status"]],
             cds_missing=g["missing"], cds_length_outlier=g["length_outlier"],
             population=POP.get(r["population"], r["population"]),
             assembly_md5=hashlib.md5(seq.encode()).hexdigest())
    if s == "K2":   # library composition and population recomputed on the pooled reads
        k = json.load(open(R + "/38_k2/K2.metrics.json"))
        r["pct_IBV_reads"] = k["pct_IBV"]; r["population"] = POP[k["call"]]
        r["sra_run"] = "SRR40942476; SRR40981834"   # both libraries were used
    x = tree[s]
    r.update(lineage_S1=x["lineage"], lineage_S1_support=x["lineage_support"],
             sublineage_S1=x["sublineage"], sublineage_S1_support=x["sublineage_support"])
    out.append(r)
with open(F + "/SupplementaryTable1_per_genome.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(out)

import collections
L = [int(r["total_bp"]) for r in out]; D = [int(r["median_depth"]) for r in out]
print("genomes %d" % len(out))
print("length median %d, range %d-%d" % (statistics.median(L), min(L), max(L)))
print("depth median %d, min %d; over 100x: %d" % (statistics.median(D), min(D), sum(1 for d in D if d > 100)))
print("with N: %d genomes; total N %d" % (sum(1 for r in out if int(r["ambiguous_bp"]) > 0), sum(int(r["ambiguous_bp"]) for r in out)))
z = [int(r["zero_cov_bp"]) for r in out]
print("genomes with any uncovered base: %d; largest %d bp" % (sum(1 for x in z if x > 0), max(z)))
for k in ("assembly_route", "annotation_status", "population", "cds_annotated"):
    print(k, dict(collections.Counter(r[k] for r in out).most_common()))
print("vaccines:", [r["sample_name"] for r in out if r["sample_name"] in ("H120", "AVR1-08", "K40-09", "QXE", "K2")])
