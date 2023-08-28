
# c_munge_pops_features.sh


# Load modules
module load Python/3.10.4-GCCcore-11.3.0

# Create an output directory, navigate there
OUTDIR=/home/heilbron/projects/pops/data/features_munged
mkdir "$OUTDIR"
cd "$OUTDIR"

# Arguments
SCRIPT=/home/heilbron/repos/pops/munge_feature_directory.py
GENE_LOC_FILE=/home/heilbron/projects/pops/data/gene_locations.tsv
SPLIT_FEATURES_DIR=/home/heilbron/projects/pops/data/features_split
OUTFILE_PREFIX=features

# Run
python "$SCRIPT" \
  --gene_annot_path "$GENE_LOC_FILE" \
  --feature_dir "$SPLIT_FEATURES_DIR" \
  --save_prefix "$OUTFILE_PREFIX" \
  --nan_policy raise

