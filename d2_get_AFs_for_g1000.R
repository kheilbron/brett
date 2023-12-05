
#   d_get_AFs_for_g1000.R


# Before firing up R
# module load PLINK/2.00a3.6-GCC-11.3.0

# Make scratch directories to work in
scratch_dir <- file.path( "/scratch-shared", Sys.getenv("USER"), "high_quality_g1000" )
dir.create( path=scratch_dir,  showWarnings=FALSE, recursive=TRUE )

# # If output exists, skip
common_outfile   <- file.path( scratch_dir, "g1000_eur_snps_mac_ge_10.tsv" )
if( file.exists(common_outfile) ){
  stop("Output already exists")
}

# Run plink2 to get EAF of HRC SNPs
bfile     <- "/home/heilbron/projects/pops/data/g1000_eur/g1000_eur"
plink_out <- file.path( scratch_dir, "plink" )
plink_cmd <- paste( "plink2 --bfile", bfile, "--freq --out", plink_out )
system(plink_cmd)

# Read in AF file, calculate allele count
library(data.table)
af <- fread( paste0( plink_out, ".afreq" ) )
names(af) <- c( "chr", "snp", "ref", "alt", "af", "n" )
af$ac <- round( af$af * af$n )

# Read in .bim file
bim <- fread( paste0( bfile, ".bim" ) )
names(bim) <- c( "chr", "snp", "cm", "bp", "alt", "ref" )

# Merge
gcols_bim <- c( "snp", "chr", "bp", "ref", "alt" )
gcols_af  <- c( "af", "ac" )
m <- cbind( bim[ , ..gcols_bim ],
            af[  , ..gcols_af  ] )

# Remove MAC < 10
# Remove chrX
common_df <- m[ m$ac >= 10 & m$chr != 23 , ]

# Write
fwrite( x=common_df, file=common_outfile, sep="\t" )

# Move back to my home directory
common_home_dir_file   <- "/home/heilbron/projects/pops/data/g1000_eur_snps_mac_ge_10.tsv"
file.copy( from=common_outfile, to=common_home_dir_file )






