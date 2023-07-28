
#   d_get_AFs_for_HRC.R


# Before firing up R
# module load PLINK/2.00a3.6-GCC-11.3.0

# Make scratch directories to work in
scratch_dir <- file.path( "/scratch-local", Sys.getenv("USER"), "high_quality_HRC" )
plink_dir   <- file.path( scratch_dir, "plink_freqs" )
merged_dir  <- file.path( scratch_dir, "merged_files" )
dir.create( path=plink_dir,  showWarnings=FALSE, recursive=TRUE )
dir.create( path=merged_dir, showWarnings=FALSE, recursive=TRUE )

# Function: make a .bim file including AF and AC
bim_w_af_and_ac <- function(chromosome){
  
  # If output exists, skip
  merged_file <- file.path( merged_dir, paste0( chromosome, ".tsv" ) )
  if( file.exists(merged_file) ){
    message( "Output already exists for chromosome ", chromosome, ", skipping" )
    return()
  }
  
  # Run plink2 to get EAF of HRC SNPs
  bfile     <- paste0( "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/HRC.r1-1.EGA.GRCh37.chr", 
                       chromosome, ".impute.plink.EUR" )
  plink_out <- file.path( plink_dir, chromosome )
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
  
  # Remove MAC <10
  m2 <- m[ m$ac >= 10 , ]
  
  # Write
  fwrite( x=m2, file=merged_file, sep="\t" )
}

# Loop through chromosomes making improved .bim files
for( chromosome in 22:1 ){
  message("\n\n\n")
  message( "Starting chromosome: ", chromosome )
  bim_w_af_and_ac(chromosome)
}

# Concatenate files
setwd(merged_dir)
merged_files <- list.files( path=merged_dir, full.names=FALSE )
cat_files    <- paste( merged_files, collapse=" " )
outfile <- file.path( scratch_dir, "hrc_eur_snps_mac_ge_10.tsv" )
system( paste0( "head -n1 ", merged_files[1], " > ", outfile ) )
system( paste0( 'cat ', cat_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', outfile ) )

# Move back to my home directory
home_dir_file <- "/home/heilbron/projects/pops/data/hrc_eur_snps_mac_ge_10.tsv"
file.copy( from=outfile, to=home_dir_file )






