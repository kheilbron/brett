
pops_ma <- function( outdir, 
                     pops_dir_pop1,
                     pops_dir_pop2,
                     loci_dir_ma ){
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/brett/z_brett.R")
  message2("Load libraries and sources")
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(dplyr) )
  suppressPackageStartupMessages( library(rmarkdown) )
  
  # Create an output directory
  message2("Create an output directory")
  dir.create( path=outdir, showWarnings=FALSE, recursive=TRUE )
  
  # Copy from loci_dir_ma: the GWAS file
  message2("Copy from loci_dir_ma: the GWAS file")
  file.copy( from = file.path( loci_dir_ma, "gwas_sumstats.tsv" ), 
             to   = file.path( outdir,      "gwas_sumstats.tsv" ) )
  
  # Read in POPS for both populations
  message2("Read in POPS for both populations")
  pops1 <- fread( file.path( pops_dir_pop1, "pops.preds" ) )
  pops2 <- fread( file.path( pops_dir_pop2, "pops.preds" ) )
  
  # Read in MAGMA for both populations (just to get per-gene N)
  message2("Read in MAGMA for both populations (just to get per-gene N)")
  mag1 <- fread( file.path( pops_dir_pop1, "magma.genes.out" ) )
  mag2 <- fread( file.path( pops_dir_pop2, "magma.genes.out" ) )
  
  # Add MAGMA N to POPS data
  message2("Add MAGMA N to POPS data")
  pops1$N1 <- mag1$N[ match( pops1$ENSGID, mag1$GENE ) ]
  pops2$N2 <- mag2$N[ match( pops2$ENSGID, mag2$GENE ) ]
  
  # Mean impute missing N values
  message2("Mean impute missing N values")
  pops1$N1[ is.na(pops1$N1) ] <- mean( pops1$N1, na.rm=TRUE )
  pops2$N2[ is.na(pops2$N2) ] <- mean( pops2$N2, na.rm=TRUE )
  
  # Merge POPS datasets
  message2("Merge POPS datasets")
  names(pops1)[ names(pops1) == "PoPS_Score" ] <- "pops1"
  names(pops2)[ names(pops2) == "PoPS_Score" ] <- "pops2"
  pops3 <- full_join( x = pops1[ , c( "ENSGID", "pops1", "N1" ) ], 
                      y = pops2[ , c( "ENSGID", "pops2", "N2" ) ] )
  
  # Merge MAGMA datasets
  message2("Merge MAGMA datasets")
  names(mag1)[ names(mag1) == "ZSTAT" ] <- "ZSTAT1"
  names(mag2)[ names(mag2) == "ZSTAT" ] <- "ZSTAT2"
  names(mag1)[ names(mag1) == "N" ]     <- "N1"
  names(mag2)[ names(mag2) == "N" ]     <- "N2"
  mag1$P <- mag2$P <- mag1$NSNPS <- mag2$NSNPS <- mag1$NPARAM <- mag2$NPARAM <- NULL
  mag3   <- full_join( x = mag1, 
                       y = mag2 )
  mag3$N1[     is.na(mag3$N1)     ] <- 0
  mag3$N2[     is.na(mag3$N2)     ] <- 0
  mag3$ZSTAT1[ is.na(mag3$ZSTAT1) ] <- 0
  mag3$ZSTAT2[ is.na(mag3$ZSTAT2) ] <- 0
  
  # Meta-analysis function
  n_based_ma <- function( n, z ){
    w     <- sqrt(n)
    num   <- sum(z*w)
    denom <- sqrt( sum(w^2) )
    return( num / denom )
  }
  
  # Meta-analyze POPS
  message2("Meta-analyze POPS")
  pops3$PoPS_Score <- apply( X=pops3[,-1], MARGIN=1, FUN=function(x){
    n <- c( x["N1"],    x["N2"] )
    z <- c( x["pops1"], x["pops2"] )
    n_based_ma( n=n, z=z )
  })
  
  # Meta-analyze MAGMA
  message2("Meta-analyze MAGMA")
  mag3$N     <- mag3$N1 + mag3$N2
  mag3$ZSTAT <- apply( X=mag3[,-1], MARGIN=1, FUN=function(x){
    n <- c( x["N1"],     x["N2"] )
    z <- c( x["ZSTAT1"], x["ZSTAT2"] )
    n_based_ma( n=n, z=z )
  })
  mag3$P <- z_to_p( z=mag3$ZSTAT )
  
  # Write the new pops.preds file into the output directory
  message2("Write the new pops.preds file into the output directory")
  fwrite( x=pops3, file=file.path( outdir, "pops.preds" ), sep="\t" )
  
  # Write the new magma.genes.out file into the output directory
  message2("Write the new magma.genes.out file into the output directory")
  fwrite( x=mag3, file=file.path( outdir, "magma.genes.out" ), sep="\t" )
  
  # Run pops_plots
  message2("Run pops_plots")
  pops_plots( maindir  = outdir, 
              loci_dir = loci_dir_ma )
  
  # Run magma_plots
  message2("Run magma_plots")
  magma_plots( maindir  = outdir, 
               loci_dir = loci_dir_ma )
  
  # Render an HTML report
  message2("Render an HTML report")
  rmd_file  <- file.path( outdir, "report.Rmd" )
  html_file <- file.path( outdir, "report.html" )
  file.copy( from = "/projects/0/prjs0817/repos/brett/g_brett_template.Rmd",
             to   = rmd_file, overwrite=TRUE )
  args <- list( maindir    = outdir,
                ld_panel   = "hrc",
                gw_file    = file.path( outdir, "gwas_sumstats.tsv" ),
                loci_dir   = loci_dir_ma,
                check_args = FALSE )
  render( input       = rmd_file, 
          params      = args, 
          output_file = html_file )
  
  # Done
  message2("Done")
}

