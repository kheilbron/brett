#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=rome
#SBATCH --time 59:00
#SBATCH --mem=27G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=karl.heilbron@charite.de


# Parse arguments
MAINDIR=$1

# Set variables
pops_path=/projects/0/prjs0817/repos/pops/pops2.py
gene_locs=/projects/0/prjs0817/projects/pops/data/gene_locations.tsv
feature_prefix=/projects/0/prjs0817/projects/pops/data/features_munged/features
magma_prefix="${MAINDIR}"/magma
control_path=/projects/0/prjs0817/repos/pops/example/data/utils/features_jul17_control.txt
pops_prefix="${MAINDIR}"/pops

# Print arguments
echo "Using POPS path: ${pops_path}"
echo "Using gene location file: ${gene_locs}"
echo "Using feature prefix: ${feature_prefix}"
echo "Using MAGMA prefix: ${magma_prefix}"
echo "Using control features path: ${control_path}"
echo -e "Using output prefix: ${pops_prefix}\n"

# Run POPS
python "${pops_path}" \
  --gene_annot_path "${gene_locs}" \
  --feature_mat_prefix "${feature_prefix}" \
  --num_feature_chunks 12 \
  --magma_prefix "${magma_prefix}" \
  --control_features_path "${control_path}" \
  --out_prefix "${pops_prefix}" \
  --verbose








