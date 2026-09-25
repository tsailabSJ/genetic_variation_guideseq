#!/bin/bash
# Genome-wide association test (Methods: "Genetic data analysis").
#
# PROVENANCE: original. data/original_run.sh is the run.sh from association_tests/,
# verbatim. The commands below are the same, extended with the GUIDE-seq-2 read-count
# phenotype (which was run separately) and pointed at the copies of the inputs kept in
# data/ so the step is self-contained.
#
# plink_glm.py (copied unmodified from HemTools/bin) runs one PLINK v2 --glm linear
# regression per off-target site, in parallel over the phenotype columns. For each site it
#   1. parses the column name "<chr>_<start>_<end>" back into coordinates,
#   2. writes a temporary BED of that window PADDED BY 100 bp ON EACH SIDE, and
#   3. runs  plink2 --vcf <vcf> --extract bed0 <tmp.bed> --pheno <pheno.tsv>
#                   --pheno-name <site> --glm --output-chr 26 --freq counts --out <jid>/<site>
#
# The +/-100 bp padding is essential -- it is why the VCF is named "flank100". Most tested
# variants lie in the flanks, not inside the 23 bp protospacer: in data/example.glm.linear
# all three reported variants are outside it.
#
# Additive genotype coding (PLINK's default TEST=ADD), 95 donors, no covariates.
# BETA is the effect size and P the association p value used throughout Fig 4.
#
# These three output folders are the ones behind Supplementary Table 4 (verified: their
# P and BETA reproduce the master table; the other phenotype runs in association_tests/
# -- raw counts, normalized indel, primer-SNP-masked indel -- do not).
#
# Requires plink2 and plink_glm.py on $PATH.  Runtime: ~1 h for all three on 64 cores.

set -euo pipefail
cd "$(dirname "$0")"
VCF=${1:-data/donor95.1115.flank100.final.vcf}
export PATH="$PWD:$PATH"

plink_glm.py -f data/OT_LOG2_count_matrix.pheno.tsv      -v "$VCF" -j OT_LOG2_count_matrix_PLINK_GLM_results
plink_glm.py -f data/rhampseq_indel_freq.pheno.tsv       -v "$VCF" -j rhampseq_indel_freq_PLINK_GLM_results
plink_glm.py -f data/rhampseq_integration_freq.pheno.tsv -v "$VCF" -j rhampseq_integration_freq_PLINK_GLM_results
