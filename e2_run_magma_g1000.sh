#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH --time 30:00
#SBATCH --mem=7G
# #SBATCH --mail-type=BEGIN,END
# #SBATCH --mail-user=karl.heilbron@charite.de

# Parse arguments
MAINDIR=$1

# Set variables
MAGMA=/home/heilbron/software/magma/magma
SNPS_MAPPED_TO_GENES="${MAINDIR}"/snps_mapped_to_genes_filtered.genes.annot
GWAS_PVALUES="${MAINDIR}"/gwas_pvalues.tsv
BED_PATH=/home/heilbron/projects/pops/data/g1000_eur/g1000_eur
OUT_PREFIX="${MAINDIR}"/magma

# Print arguments
echo "Using BED files: ${BED_PATH}"
echo "Using MAGMA binary: ${MAGMA}"
echo "Using SNP map: ${SNPS_MAPPED_TO_GENES}"
echo "Using GWAS P values: ${GWAS_PVALUES}"
echo -e "Using outfile prefix: ${OUT_PREFIX}\n"

# Run MAGMA
"${MAGMA}" \
  --bfile "${BED_PATH}" \
  --gene-annot "${SNPS_MAPPED_TO_GENES}" \
  --pval "${GWAS_PVALUES}" ncol=N \
  --gene-model snp-wise=mean \
  --out "${OUT_PREFIX}"





