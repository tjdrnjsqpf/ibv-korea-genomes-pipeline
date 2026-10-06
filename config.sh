# Settings every script reads. Export these, or edit the defaults here.
#
# Nothing else in this repository refers to a particular machine: the analysis ran on a shared
# Linux server with the tools in micromamba environments, but any layout works.

# Directory holding the reads and all intermediate output. Scripts cd here first.
IBV_WORK="${IBV_WORK:-$PWD}"

# Optional command prefix, if the tools are not on PATH; leave empty if they are.
#   e.g. IBV_RUN="micromamba run -n ibv_asm"
IBV_RUN="${IBV_RUN:-}"
IBV_QC_RUN="${IBV_QC_RUN:-$IBV_RUN}"          # Kraken2 lived in its own environment here
IBV_PHYLO_RUN="${IBV_PHYLO_RUN:-$IBV_RUN}"    # MAFFT, trimAl and IQ-TREE likewise

# Directory holding minimap2, samtools and trimal, if they are not on PATH.
IBV_BIN="${IBV_BIN:-}"

# NCBI table2asn (not on conda-forge/bioconda; see the README).
IBV_TABLE2ASN="${IBV_TABLE2ASN:-table2asn}"

# Kraken2 database. The analysis used PlusPF-16GB (2026-06-26).
KRAKEN2_DB="${KRAKEN2_DB:-}"

export IBV_WORK IBV_RUN IBV_QC_RUN IBV_PHYLO_RUN IBV_BIN IBV_TABLE2ASN KRAKEN2_DB
