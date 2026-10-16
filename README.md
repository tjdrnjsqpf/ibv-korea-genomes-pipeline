# Analysis pipeline for 183 infectious bronchitis virus genomes from the Republic of Korea

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23179043.svg)](https://doi.org/10.5281/zenodo.23179043)

The analysis behind the data descriptor *Whole-genome sequences of 183 infectious bronchitis
virus isolates from national poultry surveillance in the Republic of Korea, 1986-2022*: from
the deposited reads to the genomes, the per-genome quality record and the figures.

The data are at NCBI under **BioProject PRJNA1537317** - BioSample SAMN63679324-SAMN63679506,
SRA study SRP741337 (184 runs; the vaccine strain K2 has two) and GenBank PV902964-PV903161.
The dataset is defined by the accession list in Supplementary Table 1, not by the accession
ranges, which also contain unrelated records.

## Pipeline

Stages run in order. Everything after `2_assembly/02_assemble.sh` reads only that stage's
output, so any later stage can be re-run on its own.

### 1_reads

| Script | What it does |
|---|---|
| `01_check_reads.sh` | Verifies every read file before anything reads it: gzip integrity, record counts, R1/R2 agreement. A truncated fastq assembles without complaint, so this runs first. |

### 2_assembly

| Script | What it does |
|---|---|
| `02_assemble.sh` | fastp then MEGAHIT, placing each assembly against a seven-genome reference panel and keeping the contigs the reads support. No reference scaffolding and no N-padding: a reference only selects and orders contigs, and never contributes sequence. |
| `02b_join_overlapping_contigs.py` | Joins neighbouring contigs sharing an exact 141-nt terminal overlap, the assembler's final k-mer left on both sides where the graph branched, and checks each junction against the reads. |
| `02d_linearise_circular_contig.py` | One contig was assembled as a circular path (identical first and last 141 nt, halves in reverse order on the reference); it is cut and written in linear order. |
| `02e_reassemble_relaxed.sh` | Re-assembles the non-overlapping cases with `--min-count 1 --k-min 21 --k-step 10` to test whether the reads close the junction. |
| `02c_join_across_gaps.py` | Adopts a junction closed by re-assembly, and joins the rest with a declared run of N whose length is estimated from the reference. |
| `02f_K2_pooled_libraries.sh` | Pools the two libraries of the vaccine strain K2, neither of which assembles alone. |

### 3_annotation

| Script | What it does |
|---|---|
| `03_annotate.py` | Projects the CDS coordinates of each isolate's annotation reference onto its assembly, refines start codons on the sequence itself, and locates the -1 ribosomal frameshift per genome on the slippery sequence UUUAAAC. |
| `04_validate_table2asn.sh` | Validates the annotation with NCBI table2asn. |

### 4_qc

| Script | What it does |
|---|---|
| `05_gene_qc.py` | Judges each CDS against the dataset's own median length per gene, because reference records differ in annotation convention by as much as 51 nt for the same gene. |
| `06_population_screen.py` | Minor-allele frequency spectrum per genome: clonal, low-level diversity, or mixed population. |
| `07_kraken.sh`, `08_kraken_summary.py` | Read classification: what each library held besides IBV. |
| `09_build_manifest.py`, `09b_rebuild_table1.py` | Assemble the per-genome record that becomes Supplementary Table 1. |

### 5_lineage

| Script | What it does |
|---|---|
| `15b_S1_tree.sh` | Aligns the first 1,700 nt of spike onto the reference framework (the 32 prototypes of the harmonised classification, every prototype of the four genotypes described since it, and 103 labelled Korean and other Asian sequences), trims with trimAl `-gt 0.9` and infers a maximum-likelihood tree with IQ-TREE. |
| `15_assign_lineage_S1.py` | Reads lineage and sub-lineage off that tree: the nearest ancestral node whose reference descendants all carry one label, with SH-aLRT at least 80 % and ultrafast bootstrap at least 95 %; otherwise unresolved. |
| `15c_lineage_identity.py` | Identity of each genome to the prototype of the lineage it was assigned. A tree cannot place a genome in a lineage its framework lacks, and support does not reveal the gap, so this is reported alongside the label. |

### 6_figure_data

| Script | What it does |
|---|---|
| `16_figure_source_data.py` | Derives the values behind every panel of Figures 1 and 2 from the per-genome record and the sample metadata: genomes per year by sub-lineage, genomes per county and per province, read depth, assembly length by finishing route, library composition, and coding-sequence status per gene. One file per panel; these are the figure source data published with the descriptor. |

The figures themselves were drawn by hand from these files, so no plotting code is released.

`samples_all.tsv` maps each sample to its read-file prefix and sequencing run; every script
keys on it.

## Running it

```bash
micromamba env create -f environment.yml      # or conda/mamba
micromamba activate ibv_asm
bash data/fetch_maps.sh                       # boundary files for Figure 1b
```

`table2asn` is not on conda-forge or bioconda; take the binary from
<https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/>. Kraken2 needs a database;
the analysis used PlusPF-16GB (2026-06-26), downloaded separately.

Settings live in [config.sh](config.sh), which every shell script sources and which the Python
scripts read from the environment. Nothing in the repository refers to a particular machine.

| Variable | Meaning |
|---|---|
| `IBV_WORK` | Directory holding the reads and all intermediate output. Scripts `cd` here first. Defaults to the current directory. |
| `IBV_RUN` | Command prefix if the tools are not on `PATH`, for example `micromamba run -n ibv_asm`. Empty by default. |
| `IBV_QC_RUN`, `IBV_PHYLO_RUN` | The same for Kraken2 and for MAFFT, trimAl and IQ-TREE, which lived in their own environments here. Both default to `IBV_RUN`. |
| `IBV_BIN` | Directory holding `minimap2`, `samtools` and `trimal` if they are not on `PATH`. |
| `IBV_TABLE2ASN` | The `table2asn` binary. Defaults to `table2asn`. |
| `KRAKEN2_DB` | Kraken2 database directory. |

The boundary files for the map are third-party data under their own terms and are not
redistributed here; `data/fetch_maps.sh` fetches one and says where to get the other.

Reads are not in this repository; download them from SRA, for example

```bash
fasterq-dump --split-files SRR40942385
```

## Licence and citation

MIT, see [LICENSE](LICENSE). If you use these scripts, please cite the data descriptor, and the
code as <https://doi.org/10.5281/zenodo.23179043> (all versions) or
<https://doi.org/10.5281/zenodo.23179044> (v1.0.0). See [CITATION.cff](CITATION.cff).
