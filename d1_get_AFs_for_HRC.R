
#   d_get_AFs_for_HRC.R


# Before firing up R
# module load PLINK/2.00a3.6-GCC-11.3.0

# Make scratch directories to work in
scratch_dir <- file.path( "/scratch-shared", Sys.getenv("USER"), "high_quality_HRC" )
plink_eur_dir   <- file.path( scratch_dir,"eur", "plink_freqs" )
rare_eur_dir    <- file.path( scratch_dir,"eur" ,"rare_snp_files" )
common_eur_dir  <- file.path( scratch_dir, "eur", "common_snp_files" )

plink_eas_dir   <- file.path( scratch_dir,"eas", "plink_freqs" )
rare_eas_dir    <- file.path( scratch_dir,"eas" ,"rare_snp_files" )
common_eas_dir  <- file.path( scratch_dir, "eas", "common_snp_files" )


dir.create( path=plink_eur_dir,  showWarnings=FALSE, recursive=TRUE )
dir.create( path=rare_eur_dir,   showWarnings=FALSE, recursive=TRUE )
dir.create( path=common_eur_dir, showWarnings=FALSE, recursive=TRUE )


dir.create( path=plink_eas_dir,  showWarnings=FALSE, recursive=TRUE )
dir.create( path=rare_eas_dir,   showWarnings=FALSE, recursive=TRUE )
dir.create( path=common_eas_dir, showWarnings=FALSE, recursive=TRUE )


# Function: make a .bim file including AF and AC
bim_w_af_and_ac <- function( chromosome, population ){
  
  # If output exists, skip
  rare_file   <- file.path( scratch_dir,population,"rare_snp_files",  paste0( chromosome, ".tsv" ) )
  common_file <- file.path( scratch_dir,population, "common_snp_files",  paste0( chromosome, ".tsv" ) )
  
  
  if( all( file.exists( rare_file, common_file ) ) ){
    message( "Output already exists for chromosome ", chromosome, ", skipping" )
    return()
  }
  
  # Run plink2 to get EAF of HRC SNPs
  if( population == "eur" ){  
    bfile     <- paste0( "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/HRC.r1-1.EGA.GRCh37.chr",
                         chromosome, ".impute.plink.EUR" )
    
  }else if(population == "eas"){
    bfile     <- paste0( "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/HRC.r1-1.EGA.GRCh37.chr",
                         chromosome, ".impute.plink.EAS" )
    
  }else{
    stop("population must be either eur or eas.")
  }
  
  plink_out <- file.path(scratch_dir,population,"plink_freqs",chromosome )
  
  
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
  # Remove MAF < 1%
  rare_df   <- m[ m$ac >= 10 , ]
  common_df <- m[ m$af >= 0.01 , ]
  
  # Write
  fwrite( x=rare_df,   file=rare_file,   sep="\t" )
  fwrite( x=common_df, file=common_file, sep="\t" )
}


# Loop through chromosomes making improved .bim files
for(chromosome in 1:22 ){
  message("\n\n\n")
  #  message( "Starting chromosome: ", chromosome, ", Population: European")
  #  bim_w_af_and_ac(chromosome, population = "eur") 
  message( "Starting chromosome: ", chromosome, ", Population: East-Asian")
  bim_w_af_and_ac(chromosome, population = "eas")
}


# Prepare to concatenate files
rare_eur_snp_files   <- list.files( path=rare_eur_dir,   full.names=FALSE )
common_eur_snp_files <- list.files( path=common_eur_dir, full.names=FALSE )
rare_eas_snp_files   <- list.files( path=rare_eas_dir,   full.names=FALSE )
common_eas_snp_files <- list.files( path=common_eas_dir, full.names=FALSE )


cat_rare_eur_files   <- paste( rare_eur_snp_files,   collapse=" " )
cat_common_eur_files <- paste( common_eur_snp_files, collapse=" " )
rare_eur_outfile     <- file.path( scratch_dir, "hrc_eur_snps_mac_ge_10.tsv" )
common_eur_outfile   <- file.path( scratch_dir, "hrc_eur_snps_maf_ge_0.01.tsv" )

cat_rare_eas_files   <- paste( rare_eas_snp_files,   collapse=" " )
cat_common_eas_files <- paste( common_eas_snp_files, collapse=" " )
rare_eas_outfile     <- file.path( scratch_dir, "hrc_eas_snps_mac_ge_10.tsv" )
common_eas_outfile   <- file.path( scratch_dir, "hrc_eas_snps_maf_ge_0.01.tsv" )


# Concatenate rare files
setwd(rare_eur_dir)
system( paste0( "head -n1 ", rare_eur_snp_files[1], " > ", rare_eur_outfile ) )
system( paste0( 'cat ', cat_rare_eur_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', rare_eur_outfile ) )

setwd(rare_eas_dir)
system( paste0( "head -n1 ", rare_eas_snp_files[1], " > ", rare_eas_outfile ) )
system( paste0( 'cat ', cat_rare_eas_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', rare_eas_outfile ) )


# Concatenate common files
setwd(common_eur_dir)
system( paste0( "head -n1 ", common_eur_snp_files[1], " > ", common_eur_outfile ) )
system( paste0( 'cat ', cat_common_eur_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', common_eur_outfile ) )


setwd(common_eas_dir)
system( paste0( "head -n1 ", common_eas_snp_files[1], " > ", common_eas_outfile ) )
system( paste0( 'cat ', cat_common_eas_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', common_eas_outfile ) )



# Move back to my home directory

rare_eur_home_dir_file   <- "/projects/0/prjs0817/projects/pops/analyses/scz/w3/eur/hrc_eur_snps_mac_ge_10.tsv"
common_eur_home_dir_file <- "/projects/0/prjs0817/projects/pops/analyses/scz/w3/eur/hrc_eur_snps_maf_ge_0.01.tsv"

rare_eas_home_dir_file   <- "/projects/0/prjs0817/projects/pops/analyses/scz/w3/eas/hrc_eas_snps_mac_ge_10.tsv"
common_eas_home_dir_file <- "/projects/0/prjs0817/projects/pops/analyses/scz/w3/eas/hrc_eas_snps_maf_ge_0.01.tsv"

file.copy( from=rare_eur_outfile,   to=rare_eur_home_dir_file )
file.copy( from=common_eur_outfile, to=common_eur_home_dir_file )


file.copy( from=rare_eas_outfile,   to=rare_eas_home_dir_file )
file.copy( from=common_eas_outfile, to=common_eas_home_dir_file )






