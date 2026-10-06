#!/usr/bin/env python
"""Source data for every figure panel, one CSV each, from the supplementary tables.
usage: 16_figure_source_data.py <SciData_Submission directory>"""
import sys, csv, re, json, collections
base = sys.argv[1]; OUT = f"{base}/02_figures/source_data"
t1 = list(csv.DictReader(open(f"{base}/03_tables/SupplementaryTable1_per_genome.tsv"), delimiter="\t"))
t2 = {r["sample_name"]: r for r in csv.DictReader(open(f"{base}/03_tables/SupplementaryTable2_sample_metadata.tsv"), delimiter="\t")}
gq = {r["strain"]: r for r in csv.DictReader(open(f"{base}/05_metadata/gene_qc_final.tsv"), delimiter="\t")}
prov = json.load(open(f"{base}/06_scripts/skorea-provinces-geo.json"))
muni = json.load(open(f"{base}/06_scripts/skorea-municipalities-2013-geo.json"))
# a few name_eng values in the 2013 file lack the usual hyphenation; Cheonan, Seongnam and Yongin
# are split into districts there, which we report at city level
CNAME = {"Dangjinsi": "Dangjin-si", "Pocheonsi": "Pocheon-si", "Hwaseongsi": "Hwaseong-si", "Gwangju": "Gwangju-si",
         "Cheonansidongnamgu": "Cheonan-si", "Cheonansiseobukgu": "Cheonan-si", "Seongnamsijungwongu": "Seongnam-si",
         "Seongnamsisujeonggu": "Seongnam-si", "Seongnamsibundanggu": "Seongnam-si", "Yonginsigiheunggu": "Yongin-si",
         "Yonginsisujigu": "Yongin-si", "Yonginsicheoingu": "Yongin-si"}

def write(name, header, rows):
    with open(f"{OUT}/{name}", "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print("%-36s %d rows" % (name, len(rows)))

GROUP = {"QX-like I": "QX-like", "QX-like II": "QX-like", "QX-like III": "QX-like", "QX-like IV": "QX-like",
         "QXE": "QX-like", "B4-like": "B4-like", "KM91-like": "KM91-like", "K40/09-like": "K40/09-like", "Mass-like": "Mass-like"}
grp = lambda r: GROUP.get(r["prior_sublineage"], "Other or not assigned")
year = lambda s: (re.match(r"(19|20)\d{2}", t2[s]["collection_date"]) or [None])[0]

# Figure 1a
ORDER = ["QX-like", "B4-like", "KM91-like", "K40/09-like", "Mass-like", "Other or not assigned"]
cnt = collections.Counter((year(r["sample_name"]), grp(r)) for r in t1 if year(r["sample_name"]))
years = sorted({y for y, _ in cnt})
write("Fig1a_genomes_per_year_by_sublineage_group.csv", ["year"] + ORDER + ["total"],
      [[y] + [cnt[(y, g)] for g in ORDER] + [sum(cnt[(y, g)] for g in ORDER)] for y in years])
write("Fig1a_per_genome.csv", ["sample_name", "collection_year", "prior_sublineage", "group"],
      [[r["sample_name"], year(r["sample_name"]) or "", r["prior_sublineage"], grp(r)] for r in t1])

# Figure 1b
def rings(g): return [p[0] for p in (g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]])]
def inside(lon, lat, ring):
    ok = False
    for i in range(len(ring) - 1):
        x0, y0 = ring[i]; x1, y1 = ring[i + 1]
        if (y0 > lat) != (y1 > lat) and lon < (x1 - x0) * (lat - y0) / (y1 - y0) + x0: ok = not ok
    return ok
NAMES = {f["properties"]["NAME_1"] for f in prov["features"]}
per = []
for r in t1:
    m = t2[r["sample_name"]]; hit, how, county = "", "", ""
    if re.match(r"\d", m["lat_lon"]):
        la, _, lo, _ = m["lat_lon"].split()
        for f in prov["features"]:
            if any(inside(float(lo), float(la), ring) for ring in rings(f["geometry"])): hit, how = f["properties"]["NAME_1"], "county coordinate"; break
        for f in muni["features"]:
            if any(inside(float(lo), float(la), ring) for ring in rings(f["geometry"])): county = CNAME.get(f["properties"]["name_eng"], f["properties"]["name_eng"]); break
    elif ": " in m["geo_loc_name"] and m["geo_loc_name"].split(": ")[1] in NAMES:
        hit, how = m["geo_loc_name"].split(": ")[1], "province recorded directly"
    if not hit: how = "vaccine strain" if m["isolation_source"] == "commercial vaccine" else "no recorded site"
    per.append([r["sample_name"], m["geo_loc_name"], m["lat_lon"], county, hit, how])
write("Fig1b_per_genome_province.csv", ["sample_name", "geo_loc_name", "lat_lon", "county_2013", "province", "placed_by"], per)
pc = collections.Counter(p[4] for p in per if p[4])
write("Fig1b_genomes_per_province.csv", ["province", "genomes"], [[k, v] for k, v in sorted(pc.items(), key=lambda x: -x[1])])
cc = collections.Counter((p[3], p[4]) for p in per if p[3])
write("Fig1b_genomes_per_county.csv", ["county_2013", "province", "genomes"],
      [[k[0], k[1], v] for k, v in sorted(cc.items(), key=lambda x: (-x[1], x[0]))])

# Figure 2
write("Fig2a_depth_per_genome.csv", ["rank", "sample_name", "median_depth"],
      [[i + 1, r["sample_name"], r["median_depth"]] for i, r in enumerate(sorted(t1, key=lambda r: int(r["median_depth"])))])
write("Fig2b_length_per_genome.csv", ["sample_name", "total_bp", "ambiguous_bp", "assembly_route", "declared_gap"],
      [[r["sample_name"], r["total_bp"], r["ambiguous_bp"], r["assembly_route"], r["declared_gap"]] for r in t1])
write("Fig2c_ibv_reads_per_genome.csv", ["sample_name", "pct_IBV_reads"], [[r["sample_name"], r["pct_IBV_reads"]] for r in t1])
GENES = ["1a", "1ab", "S", "3a", "3b", "3c", "M", "5a", "5b", "N"]
rows, per = [], []
for g in GENES:
    c = collections.Counter()
    for r in t1:
        q = gq[r["sample_name"]]
        st = ("not annotated" if g in (q["missing"] or "").split(";") else "partial" if g in (q["partial"] or "").split(";")
              else "length outlier" if any(x.split(":")[0] == g for x in (q["length_outlier"] or "").split(";") if x) else "complete")
        c[st] += 1; per.append([r["sample_name"], g, st])
    rows.append([g, c["complete"], c["length outlier"], c["partial"], c["not annotated"]])
write("Fig2d_cds_recovery_by_gene.csv", ["gene", "complete", "length_outlier", "partial", "not_annotated"], rows)
write("Fig2d_per_genome_gene_status.csv", ["sample_name", "gene", "status"], per)
