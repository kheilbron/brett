#!/bin/bash

# Set up directories and file names
main_dir=/home/heilbron/projects/pops/data
in_dir="${main_dir}"/features_raw
out_dir="${main_dir}"/features_split
in_file="${in_dir}"/PoPS.features.txt.gz
unzip_file="${in_dir}"/PoPS.features.txt
n_cols_per_out_file=1000

# Unzip the input file
gzip -dk ${in_file}

# Find number of columns in the input file
n_in_cols=$(awk 'NR==1{print NF}' "$unzip_file")

# Split the feature file into 58 files of up to 1k columns each
for((i=1, start=1, end=$n_cols_per_out_file; i < n_in_cols/n_cols_per_out_file + 2; i++, start+=n_cols_per_out_file, end+=n_cols_per_out_file)); do
  echo "Starting on file ${i}/58"
  cut -f$start-$end "${unzip_file}" > "${outdir}"/features."${i}"
done

# Delete the unzipped file
rm "${unzip_file}"
echo "Done"

