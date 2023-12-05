
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
    bfile <- file.path( "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/",
                        "HRC_reference.r1-1/pop_EUR/",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.EUR" ) )
  }else if(population == "eas"){
    bfile <- file.path( "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/",
                        "HRC_reference.r1-1/pop_EAS/",
                        paste0( "HRC.r1-1.EGA.GRCh37.chr", chromosome, 
                                ".impute.plink.EAS" ) )
  }else{
    stop("population must be either eur or eas.")
  }
  
  # Run plink2 to get EAF of HRC SNPs
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
  #   Rare:   MAC < 10
  #   Common: MAF < 1%
  if( rare.or.common == "rare" ){
    out_df   <- m[ m$ac >= 10 , ]
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
  if( population != "eur" & population != "eas" ){
    stop("population must be either: eur or eas")
  }
  
  # Check inputs: rare.or.common
  if( rare.or.common != "rare" & rare.or.common != "common" ){
    stop("rare.or.common must be either: rare or common")
  }
  
  # Make scratch directories to work in
  scratch_dir <- file.path( "/scratch-shared", Sys.getenv("USER"), "high_quality_HRC" )
  plink_dir   <- file.path( scratch_dir, population, rare.or.common, "plink_freqs" )
  snp_dir     <- file.path( scratch_dir, population, rare.or.common, "snp_files" )
  dir.create( path=plink_dir,  showWarnings=FALSE, recursive=TRUE )
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
    suffix    <- paste0( "hrc_", population, "_snps_mac_ge_10.tsv" )
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

# Run: EUR, rare
get_high_quality_snps( population     = "eur", 
                       rare.or.common = "rare" )

# Run: EUR, common
get_high_quality_snps( population     = "eur", 
                       rare.or.common = "common" )

# Run: EAS, rare
get_high_quality_snps( population     = "eas", 
                       rare.or.common = "rare" )

# Run: EAS, common
get_high_quality_snps( population     = "eas", 
                       rare.or.common = "common" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------


