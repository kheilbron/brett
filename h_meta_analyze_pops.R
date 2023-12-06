
pops_ma <- function( outdir, 
                     maindir1,
                     maindir2 ){
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/brett/z_brett.R")
  message2("Load libraries and sources")
  library(data.table)
  library(dplyr)
  library(rmarkdown)
  
  # Create an output directory
  message2("Create an output directory")
  dir.create( path=outdir, showWarnings=FALSE, recursive=TRUE )
  
  # Copy from maindir1: the peaks file, GWAS file, MAGMA file
  message2("Copy from maindir1: the peaks file, GWAS file, MAGMA file")
  file.copy( from = file.path( maindir1, "peaks.tsv" ), 
             to   = file.path( outdir,   "peaks.tsv" ) )
  file.copy( from = file.path( maindir1, "gwas_pvalues.tsv" ), 
             to   = file.path( outdir,   "gwas_pvalues.tsv" ) )
  file.copy( from = file.path( maindir1, "magma.genes.out" ), 
             to   = file.path( outdir,   "magma.genes.out" ) )
  
  # Read in POPS for both populations
  message2("Read in POPS for both populations")
  pops1 <- fread( file.path( maindir1, "pops.preds" ) )
  pops2 <- fread( file.path( maindir2, "pops.preds" ) )
  
  # Read in MAGMA for both populations (just to get per-gene N)
  message2("Read in MAGMA for both populations (just to get per-gene N)")
  mag1 <- fread( file.path( maindir1, "magma.genes.out" ) )
  mag2 <- fread( file.path( maindir2, "magma.genes.out" ) )
  
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
  m <- full_join( x = pops1[ , c( "ENSGID", "pops1", "N1" ) ], 
                  y = pops2[ , c( "ENSGID", "pops2", "N2" ) ] )
  
  # Meta-analysis function
  n_based_ma <- function( n, z ){
    w     <- sqrt(n)
    num   <- sum(z*w)
    denom <- sqrt( sum(w^2) )
    return( num / denom )
  }
  
  # Meta-analyze
  message2("Meta-analyze")
  m$PoPS_Score <- apply( X=m[,-1], MARGIN=1, FUN=function(x){
    n <- c( x["N1"],    x["N2"] )
    z <- c( x["pops1"], x["pops2"] )
    n_based_ma( n=n, z=z )
  })
  
  # Write the new pops.preds file into the output directory
  message2("Write the new pops.preds file into the output directory")
  fwrite( x=m, file=file.path( outdir, "pops.preds" ), sep="\t" )
  
  # Run pops_plots to get plots
  message2("Run pops_plots to get plots")
  pops_plots( maindir = outdir, 
              z.or.p  = "z" )
  
  # Render an HTML report
  message2("Render an HTML report")
  rmd_file  <- file.path( outdir, "report.Rmd" )
  html_file <- file.path( outdir, "report.html" )
  file.copy( from = "/projects/0/prjs0817/repos/brett/g_brett_template.Rmd",
             to   = rmd_file, overwrite=TRUE )
  args <- list( maindir    = outdir,
                ld.panel   = "",
                gw.file    = file.path( maindir1, "gwas_pvalues.tsv" ),
                chr.bp.col = "",
                chr.col    = "",
                bp.col     = "",
                a1.col     = "",
                a2.col     = "",
                p.col      = "",
                eaf.col    = "",
                n1.col     = "",
                n0.col     = "",
                n.col      = "",
                n          = "",
                z.or.p     = "",
                check.args = "" )
  render( input       = rmd_file, 
          params      = args, 
          output_file = html_file )
  
  # Done
  message2("Done")
}

