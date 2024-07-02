
#   d_get_AFs_for_HRC.R


# Before firing up R
# module load PLINK/2.00a3.6-GCC-11.3.0


#-------------------------------------------------------------------------------
#   bim_w_af_and_ac: Make a .bim file including AF and AC
#-------------------------------------------------------------------------------

bim_w_af_and_ac <- function( chromosome, population, rare.or.common ){
  
  # Assign directories
  scratch_dir <- file.path( "/scratch-shared", Sys.getenv("USER"), "high_quality_HRC" )
  plink_dir   <- file.path( scratch_dir, population, rare.or.common, "plink_freqs" )
  snp_dir     <- file.path( scratch_dir, population, rare.or.common, "snp_files" )
  
  # If output exists, skip
  bim_outfile <- file.path( snp_dir, paste0( chromosome, ".tsv" ) )
  if( file.exists(bim_outfile) ){
    message( "Output already exists for chromosome ", chromosome, ", skipping" )
    return()
  }
  
  # Assign the correct BED file based on population
  if( population == "eur" ){  
    bfile <- file.path( "/gpfs/work5/0/pgcdac/imputation_references/",
                        "HRC.r1-1_merged_EUR_EAS_panel/HRC.r1-1_consistent_size/pop_EUR",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.EUR.16860" ) )
  }else if(population == "eas"){
    bfile <- file.path( "/gpfs/work5/0/pgcdac/imputation_references/",
                        "HRC.r1-1_merged_EUR_EAS_panel/HRC.r1-1_consistent_size/pop_EAS",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.EAS.538" ) )
  }else if(population == "afr"){
    bfile <- file.path( "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/",
                        "HRC_reference.r1-1/pop_AFR/",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.AFR" ) )
  }else if(population == "eur0.92_eas0.08"){
    bfile <- file.path( "/gpfs/work5/0/pgcdac/imputation_references/",
                        "HRC.r1-1_merged_EUR_EAS_panel/BIP",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.combined.EUR.6529.EAS.538" ) )
  }else if(population == "eur0.89_eas0.11"){
    bfile <- file.path( "/gpfs/work5/0/pgcdac/imputation_references/",
                        "HRC.r1-1_merged_EUR_EAS_panel/PD",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.combined.EUR.4322.EAS.538" ) )
  }else if(population == "eur0.80_eas0.20"){
    bfile <- file.path( "/gpfs/work5/0/pgcdac/imputation_references/",
                        "HRC.r1-1_merged_EUR_EAS_panel/SCZ",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.combined.EUR.2191.EAS.538" ) )
  }else{
    stop("population must be either eur, eas, afr, eur0.80_eas0.20, or eur0.89_eas0.11")
  }
  
  # Run plink2 to get EAF of HRC SNPs
  setwd(plink_dir)
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
  
  # Remove SNPs that are too rare
  #   Rare:   MAC < 5
  #   Common: MAF < 1%
  if( rare.or.common == "rare" ){
    out_df   <- m[ m$ac >= 5 , ]
  }else if( rare.or.common == "common" ){
    out_df   <- m[ m$af >= 0.01 , ]
  }else{
    stop("rare.or.common must be: rare or common")
  }
  
  # Write
  fwrite( x=out_df, file=bim_outfile, sep="\t" )
}


#-------------------------------------------------------------------------------
#   get_high_quality_snps: Write out a file of high-quality reference panel SNPs
#-------------------------------------------------------------------------------

# Turn it all into a function
get_high_quality_snps <- function( population, rare.or.common ){
  
  # Check inputs: population
  if( population != "eur" & population != "eas" & population != "afr" & 
      population != "eur0.80_eas0.20" & population != "eur0.89_eas0.11" & 
      population != "eur0.92_eas0.08" ){
    stop( "population must be one of: eur, eas, afr, eur0.80_eas0.20, ",
          "eur0.89_eas0.11, eur0.92_eas0.08" )
  }
  
  # Check inputs: rare.or.common
  if( rare.or.common != "rare" & rare.or.common != "common" ){
    stop("rare.or.common must be either: rare or common")
  }
  
  # Make scratch directories to work in
  scratch_dir <- file.path( "/scratch-shared", Sys.getenv("USER"), "high_quality_HRC" )
  plink_dir   <- file.path( scratch_dir, population, rare.or.common, "plink_freqs" )
  snp_dir     <- file.path( scratch_dir, population, rare.or.common, "snp_files" )
  dir.create( path=plink_dir, showWarnings=FALSE, recursive=TRUE )
  dir.create( path=snp_dir,   showWarnings=FALSE, recursive=TRUE )
  
  # Loop through chromosomes making improved .bim files
  for( chromosome in 1:22 ){
    message("\n\n\n")
    message( "Starting chromosome: ", chromosome )
    bim_w_af_and_ac( chromosome     = chromosome, 
                     population     = population, 
                     rare.or.common = rare.or.common )
  }
  
  # Prepare to concatenate files
  snp_files <- list.files( path=snp_dir,   full.names=FALSE )
  cat_files <- paste( snp_files, collapse=" " )
  if( rare.or.common == "rare" ){
    suffix    <- paste0( "hrc_", population, "_snps_mac_ge_5.tsv" )
  }else if( rare.or.common == "common" ){
    suffix    <- paste0( "hrc_", population, "_snps_maf_ge_0.01.tsv" )
  }else{
    stop("rare.or.common must be either: rare or common")
  }
  outfile   <- file.path( scratch_dir,                               suffix )
  home_file <- file.path( "/projects/0/prjs0817/projects/pops/data", suffix )
  
  # Concatenate rare files
  setwd(snp_dir)
  system( paste0( "head -n1 ", snp_files[1], " > ", outfile ) )
  system( paste0( 'cat ', cat_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', outfile ) )
  
  # Move back to the shared project space
  file.copy( from=outfile, to=home_file )
}


#-------------------------------------------------------------------------------
#   Run
#-------------------------------------------------------------------------------

# EUR
get_high_quality_snps( population     = "eur", 
                       rare.or.common = "rare" )
get_high_quality_snps( population     = "eur", 
                       rare.or.common = "common" )

# EAS
get_high_quality_snps( population     = "eas", 
                       rare.or.common = "rare" )
get_high_quality_snps( population     = "eas", 
                       rare.or.common = "common" )

# AFR
get_high_quality_snps( population     = "afr", 
                       rare.or.common = "rare" )
get_high_quality_snps( population     = "afr", 
                       rare.or.common = "common" )

# eur0.89_eas0.11
get_high_quality_snps( population     = "eur0.89_eas0.11", 
                       rare.or.common = "rare" )
get_high_quality_snps( population     = "eur0.89_eas0.11", 
                       rare.or.common = "common" )

# eur0.92_eas0.08
get_high_quality_snps( population     = "eur0.92_eas0.08", 
                       rare.or.common = "rare" )
get_high_quality_snps( population     = "eur0.92_eas0.08", 
                       rare.or.common = "common" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------


