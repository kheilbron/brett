
#-------------------------------------------------------------------------------------------
# a_split_feature_file.R
#   Split the large PoPS feature file of 58k gene-level features into 58 files of up to
#   1000 features each. Ensure that Ensembl ID is the first column of each file.
#-------------------------------------------------------------------------------------------

# Modules required to run R:
#   module load 2022
#   module load R/4.2.1-foss-2022a

# Set arguments
message("Set arguments")
infile <- "/home/heilbron/projects/pops/data/features_raw/PoPS.features.txt.gz"
outdir <- "/home/heilbron/projects/pops/data/features_split"
n_cols_per_out_file <- 1000

# Unzip the input file
message("Unzip the input file")
unzip_file <- sub( pattern=".gz$", replacement="", x=infile )
system( paste0( "gzip -dk ", infile ) )

# Find number of columns in the input file
message("Find number of columns in the input file")
cmd <- paste0( "head -n1 ", unzip_file, " | awk '{print NF}'" )
n_in_cols <- as.integer( system( cmd, intern=TRUE ) )

# Split the feature file into 58 files of up to 1k columns each
# Ensure that each file contains column 1 (Ensembl ID, using "cut -f 1,")
dir.create( path=outdir, showWarnings=FALSE )
for( i in 1:58 ){
  outfile <- paste0( outdir, "/features.", i )
  if( !file.exists(outfile) ){
    message( "Creating sub-file: ", i, "/58" )
    start <- 1 + (i-1)*n_cols_per_out_file
    end   <- i * n_cols_per_out_file
    cmd   <- paste0( "cut -f 1,", start, "-", end, " ", unzip_file, " > ", outfile )
    system(cmd)
    #cmd <- paste( "sbatch bash name_of_script.sh", unzip_file, start, end, i, outdir )
  }
}

# Delete the unzipped file
message("Delete the unzipped file")
file.remove(unzip_file)
message("Done")


