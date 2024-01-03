#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=rome
#SBATCH --time 25:00
#SBATCH --mem=7G

# Parse arguments
CHR=$1
MAINDIR=$2
POPULATION=$3

# Set variables
MAGMA=/projects/0/prjs0817/software/magma/magma
SNPS_MAPPED_TO_GENES="${MAINDIR}"/snps_mapped_to_genes_filtered.genes.annot
GWAS_PVALUES="${MAINDIR}"/gwas_pvalues.tsv
OUTDIR="${MAINDIR}"/magma/
OUT_PREFIX="${OUTDIR}"/chr"${CHR}"

# Set population-dependent variables: EUR
if [ "${POPULATION}" = "eur" ]
then
  HRC_DIR=/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/
  BED_PATH="${HRC_DIR}"/HRC.r1-1.EGA.GRCh37.chr"${CHR}".impute.plink.EUR
fi

# Set population-dependent variables: EAS
if [ "${POPULATION}" = "eas" ]
then
  HRC_DIR=/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/
  BED_PATH="${HRC_DIR}"/HRC.r1-1.EGA.GRCh37.chr"${CHR}".impute.plink.EAS
fi

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





