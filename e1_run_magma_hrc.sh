#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH --time 30:00
#SBATCH --mem=7G

# Parse arguments
CHR=$1
MAINDIR=$2

# Set variables
MAGMA=/projects/0/prjs0817/software2/magma/magma
HRC_DIR=/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/
SNPS_MAPPED_TO_GENES="${MAINDIR}"/snps_mapped_to_genes_filtered.genes.annot
GWAS_PVALUES="${MAINDIR}"/gwas_pvalues.tsv
BED_PATH="${HRC_DIR}"/HRC.r1-1.EGA.GRCh37.chr"${CHR}".impute.plink.EUR
OUTDIR="${MAINDIR}"/magma/
OUT_PREFIX="${OUTDIR}"/chr"${CHR}"

# Print arguments
echo "Using chromosome: ${CHR}"
echo "Using BED files: ${BED_PATH}"
echo "Using MAGMA binary: ${MAGMA}"
echo "Using SNP map: ${SNPS_MAPPED_TO_GENES}"
echo "Using GWAS P values: ${GWAS_PVALUES}"
echo -e "Using outfile prefix: ${OUT_PREFIX}\n"

# Make output directory
mkdir -p "${OUTDIR}"

# Run MAGMA
"${MAGMA}" \
  --bfile "${BED_PATH}" \
  --gene-annot "${SNPS_MAPPED_TO_GENES}" \
  --pval "${GWAS_PVALUES}" ncol=N \
  --gene-model snp-wise=mean \
  --out "${OUT_PREFIX}"





