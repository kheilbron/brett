
# z_brett.R
#   Author: Karl Heilbron
#   Created: July 25, 2023


#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# ensure_finished_jobs: Make sure all cluster jobs finish running before proceeding
# mem_used:             Memory used by an R session
# message2:             Just like message, but with the date and a gap pasted in
# message_header:       Nicely-formatted text to break up your log files into 
#                       readable chunks
# n_eff:                Get the effective sample size of a case-control GWAS
# p_to_z:               Convert a P value into a z-score
# z_to_p:               Convert a z-score into a P value


# ensure_finished_jobs: Make sure all cluster jobs finish running before proceeding
ensure_finished_jobs <- function(identifier){
  
  external.call <- paste0( "squeue | grep ", identifier, " | wc -l" )
  running.jobs  <- as.numeric( system( external.call, intern=TRUE ) ) 
  total_sleep   <- 0
  while( running.jobs > 0){
    message( "Waited for ", total_sleep, " minutes, there are still ",
             running.jobs, " jobs running")
    Sys.sleep(60)
    total_sleep <- total_sleep + 1
    running.jobs <- as.numeric( system( external.call, intern=TRUE ) )
  }
}

# mem_used: Memory used by an R session
mem_used <- function(){
  gc <- gc()
  max_used <- sum( gc[,6] )
  paste0( max_used, 'Mb of memory used in this R session' )
}

# message2: Just like message, but with the date and a gap pasted in
message2 <- function(...) message(date(), "     ", ...)

# message_header: Nicely-formatted text to break up your log files into readable chunks
message_header <- function(...){
  message( "\n\n#--------------------------------------------------" )
  message( "#   ", ...)
  message( "#--------------------------------------------------" )
}

# n_eff: Get the effective sample size of a case-control GWAS
n_eff <- function( n.cases, n.controls ){
  n.eff <- 4 / ( (1/n.cases) + (1/n.controls) )
  round(n.eff)
}

# p_to_z: Convert a P value into a z-score
p_to_z <- function( p, direction=NULL, limit=.Machine$double.xmin, log.p=FALSE ){
  
  # Set P value lower limit to avoid Inf/-Inf
  if( !is.null( limit ) )  p[ which( p < limit ) ] <- limit
  
  # Get z
  if(log.p){
    z <- -qnorm( p - log(2), log.p=TRUE )
  }else{
    z <- -qnorm(p/2)
  }
  
  # Correct sign, return
  if ( !is.null( direction) )  z <-  z * sign(direction)
  z
}

# z_to_p: Convert a z-score into a P value
z_to_p <- function( z, log.p=FALSE ){
  if(log.p){
    log(2) + pnorm( -abs(z), log.p=TRUE )
  }else{
    2*pnorm( -abs(z) )
  }
}  


#-------------------------------------------------------------------------------
#   check_arguments
#-------------------------------------------------------------------------------

check_arguments <- function( ld_panel   = NULL, 
                             population = NULL,
                             gw_file    = NULL ){
  
  # ld_panel must be either "hrc" or "g1000"
  if( ld_panel != "hrc" & ld_panel != "g1000" ) stop("ld_panel must be either 'hrc' or 'g1000'")
  
  # population must be either "eur" or "eas"
  if( population != "eur" & population != "eas" ) stop("population must be either 'eur' or 'eas'")
  
  # Does the GWAS file exist?
  if( !file.exists(gw_file) )  stop("GWAS file does not exist")
}


#-------------------------------------------------------------------------------
#   format_gwas_and_snp_loc_files
#-------------------------------------------------------------------------------

format_gwas_and_snp_loc_files <- function( maindir    = "/projects/0/prjs0817/projects/pops/analyses/pd",
                                           ld_panel   = "hrc",
                                           population = "eur",
                                           gw_file    = "/projects/0/prjs0817/projects/pops/analyses/pd/meta5_raw.tab.gz" ){
  
  #-------------------------------------------------------------------------------
  #   Input descriptions
  #-------------------------------------------------------------------------------
  
  #   maindir:    Main directory in which to store results
  #   ld_panel:   Which LD reference panel should be used? Options are either: 
  #               "hrc" or "g1000".
  #   population: Which populations should be used? Options are either: 
  #               "eur" or "eas".
  #   gw_file:    GWAS file name
  
  
  #-------------------------------------------------------------------------------
  #   Read in GWAS and HRC, format columns
  #-------------------------------------------------------------------------------
  
  # Load libraries and sources
  library(data.table)
  
  # Read in reference panel SNPs
  if( ld_panel == "hrc" ){                                  ### HRC
    rare.or.common.snps <- "common"
    if( rare.or.common.snps == "common"){                   ##  Common
      if( population == "eur" ){                            #   EUR
        message2("Read in EUR HRC SNPs with MAF >= 1%")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eur_snps_maf_ge_0.01.tsv")
      }else if( population == "eas" ){                      #   EAS
        message2("Read in EAS HRC SNPs with MAF >= 1%")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eas_snps_maf_ge_0.01.tsv")
      }else{
        stop("population must be either: eur or eas")
      }
    }else if( rare.or.common.snps == "rare"){               ##  Rare
      if( population == "eur" ){                            #   EUR
        message2("Read in EUR HRC SNPs with MAC >= 10")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eas_snps_mac_ge_10.tsv")
      }else if( population == "eas" ){                      #   EAS
        message2("Read in EAS HRC SNPs with MAC >= 10")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eas_snps_mac_ge_10.tsv")
      }else{
        stop("population must be either: eur or eas")
      }
    }else{
      stop("rare.or.common must be 'rare' or 'common'")
    }
  }else if( ld_panel == "g1000" ){                          ### 1000G
    message2("Read in 1000 Genomes SNPs with EUR MAC >= 10")
    hrc <- fread("/projects/0/prjs0817/projects/pops/data/g1000_eur_snps_mac_ge_10.tsv")
  }else{
    stop("ld_panel must be either 'hrc' or 'g1000'")
  }
  
  # Read in GWAS
  message2("Read in GWAS")
  gw <- fread(gw_file)
  
  # Re-name GWAS columns
  message2("Re-name GWAS columns")
  names(gw)[ names(gw) == "A1" ] <- "a1"
  names(gw)[ names(gw) == "A2" ] <- "a2"
  names(gw)[ names(gw) == "p"  ] <- "P"
  
  
  #-------------------------------------------------------------------------------
  #   Subset to shared SNPs
  #-------------------------------------------------------------------------------
  
  # Subset GWAS and reference panel to shared SNPs based on chromosome and position
  message2("Subset GWAS and reference panel to shared SNPs")
  snps_both <- intersect( gw$SNP, hrc$snp )
  hrc2 <- hrc[ match( snps_both, hrc$snp ) , ]
  gw2  <- gw[  match( snps_both, gw$SNP  ) , ]
  message2( "Of the ", NROW(gw), " GWAS SNPs, ", NROW(gw2), 
           " (", round( 100*NROW(gw2)/NROW(gw), 2 ), "%) were found in the reference panel" )
  message2( "Of the ", NROW(hrc), " reference panel SNPs, ", NROW(hrc2), 
           " (", round( 100*NROW(hrc2)/NROW(hrc), 2 ), "%) were found in the GWAS" )
  
  
  #-------------------------------------------------------------------------------
  #   Format and write outputs
  #-------------------------------------------------------------------------------
  
  # Dump GWAS P values
  # Column names/order: SNP, P, N
  message2("Dump GWAS P values")
  gw_outfile <- file.path( maindir, "gwas_pvalues.tsv" )
  gw_out <- gw2[ , c( "SNP", "P", "N" ) ]
  fwrite( x=gw_out, file=gw_outfile, sep="\t" )
  
  # Dump a file of SNP locations
  # No header, but column order: SNP, CHR, BP
  message2("Dump a file of SNP locations")
  snp_loc_outfile <- file.path( maindir, "snp_locations.tsv" )
  snp_loc_out <- gw2[ , c( "SNP", "chr", "bp" ) ]
  fwrite( x=snp_loc_out, file=snp_loc_outfile, sep="\t", col.names=FALSE )
  
  
  #-------------------------------------------------------------------------------
  #   Done
  #-------------------------------------------------------------------------------
}


#-------------------------------------------------------------------------------
#   map_snps_to_genes
#-------------------------------------------------------------------------------

map_snps_to_genes <- function(maindir){
  
  # Make input file paths
  snp_loc_file   <- file.path( maindir, "snp_locations.tsv" )
  outfile_prefix <- file.path( maindir, "snps_mapped_to_genes" )
  gene_loc_file  <- "/projects/0/prjs0817/projects/pops/data/gene_locations.tsv"
  
  # Run
  cmd <- paste( "/projects/0/prjs0817/software/magma/magma",
                "--annotate",
                "--snp-loc", snp_loc_file, 
                "--gene-loc", gene_loc_file, 
                "--out", outfile_prefix )
  system(cmd)
}


#-------------------------------------------------------------------------------
#   rm_genes_without_enough_snps
#-------------------------------------------------------------------------------

rm_genes_without_enough_snps <- function(maindir){
  
  # Parse the number of SNPs mapped to each gene
  map_file <- file.path( maindir, "snps_mapped_to_genes.genes.annot" )
  map   <- readLines(map_file)
  nsnps <- as.integer( system( paste0( "awk '{print NF}' ", map_file ), intern=TRUE ) ) - 2
  
  # Subset and write to file
  map2  <- map[ nsnps >= 3 ]
  outfile <- sub( pattern     = "snps_mapped_to_genes.genes.annot",
                  replacement = "snps_mapped_to_genes_filtered.genes.annot", 
                  x           = map_file )
  writeLines( text=map2, con=outfile )
}


#-------------------------------------------------------------------------------
#   run_magma
#-------------------------------------------------------------------------------

run_magma <- function( maindir, ld_panel, population, do.local=FALSE ){
  
  # Create a job identifier
  job.id <- paste0( "m", sample( x=1:999, size=1 ) )
  
  # If using HRC, run each chromosome separately
  # If using 1000 Genomes, run all at once
  if( ld_panel == "hrc" ){
    
    # HRC: loop through chromosomes 
    for( CHR in 1:22 ){
      
      # Create job name and log file name
      jobname <- paste0( job.id, ".", CHR )
      logfile <- paste0( maindir, "/logs/magma", CHR, ".log" )
      
      # Run
      if(do.local){
        message2( "Running job locally for chromosome: ", CHR )
        cmd <- paste( "/projects/0/prjs0817/repos/brett/e1_run_magma_hrc.sh",
                      CHR, maindir, population )
      }else{
        message2( "Submitting job to the cluster for chromosome: ", CHR )
        cmd <- paste( "sbatch",
                      "-J", jobname,
                      "-o", logfile,
                      "-e", logfile,
                      "/projects/0/prjs0817/repos/brett/e1_run_magma_hrc.sh",
                      CHR, maindir, population )
      }
      system(cmd)
    }
  }else if( ld_panel == "g1000" ){
    
    # 1000 Genomes
    logfile <- paste0( maindir, "/logs/magma.log" )
    if(do.local){
      message2("Running job locally")
      cmd <- paste( "/projects/0/prjs0817/repos/brett/e2_run_magma_g1000.sh",
                    maindir )
    }else{
      message2("Submitting job to the cluster")
      cmd <- paste( "sbatch",
                    "-J", job.id,
                    "-o", logfile,
                    "-e", logfile,
                    "/projects/0/prjs0817/repos/brett/e2_run_magma_g1000.sh",
                    maindir )
    }
    system(cmd)
  }else{
    stop("ld_panel must be either 'hrc' or 'g1000'")
  }
  
  # Wait until all jobs are finished before proceeding
  ensure_finished_jobs(job.id)
  
  # Check that all jobs successfully completed
  mag_chr_files <- paste0( "chr", 1:22, ".genes.out" )
  if( !file.exists(mag_chr_files) ) stop("Not all MAGMA output files exist")
}


#-------------------------------------------------------------------------------
#   run_magma_efficient
#-------------------------------------------------------------------------------

run_magma_pt <- function(maindir){
  cmd <- paste( "sbatch",
                "/projects/0/prjs0817/repos/brett/e3_run_magma_parallel_thin.sh",
                maindir )
  system(cmd)
}


#-------------------------------------------------------------------------------
#   collate_magma
#-------------------------------------------------------------------------------

collate_magma <- function(maindir){
  
  # Read in MAGMA association test files
  library(data.table)
  mag_dir <- file.path( maindir, "magma" )
  assoc_files <- list.files( path=mag_dir, pattern=".genes.out$", full.names=TRUE )
  assoc0 <- lapply( X=assoc_files, FUN=fread )
  assoc  <- do.call( rbind, assoc0 )
  assoc  <- assoc[ order( assoc$CHR, assoc$START ) , ]
  
  # Write out a collated MAGMA association test file
  assoc_outfile <- file.path( maindir, "magma.genes.out" )
  fwrite( x=assoc, file=assoc_outfile, sep="\t" )
  
  # Find and sort covariance files
  cov_files <- list.files( path=mag_dir, pattern=".genes.raw$", full.names=TRUE )
  chr       <- as.integer( sub( pattern=".*magma/chr([[:digit:]]+).genes.raw", 
                                replacement="\\1", x=cov_files ) )
  cov_files <- cov_files[ order(chr) ]
  
  # Write out a collated MAGMA covariance file
  cov_outfile <- file.path( maindir, "magma.genes.raw" )
  cmd1 <- paste( "head -n2", cov_files[1], ">", cov_outfile )
  cmd2 <- paste( "cat", paste( cov_files, collapse=" " ), "| grep -v '#' >>", cov_outfile )
  system(cmd1)
  system(cmd2)
}


#-------------------------------------------------------------------------------
#   magma_plots
#-------------------------------------------------------------------------------

magma_plots <- function( maindir, loci_dir ){
  
  # Read in peaks
  message2("Read in peaks")
  library(data.table)
  loci_file <- file.path( loci_dir, "loci_cs.tsv" )
  peaks <- fread(loci_file)
  pattern <- "^chr[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_(.*)$"
  peaks$snp <- sub( pattern=pattern, replacement="\\1", x=peaks$hit )
  # peaks <- peaks[ order( peaks$chr, peaks$bp ) , ]
  
  # Read in MAGMA files
  message2("Read in MAGMA files")
  mag_file <- file.path( maindir, "magma.genes.out" )
  mag  <- fread(mag_file)
  
  # Use z-score or P value as the plot's y-axis, depending on choice
  message2("Establish plot y-axis (z-score or -log10 P value)")
  z.or.p <- "z"
  if( z.or.p == "z" ){
    mag$Y <- abs(mag$ZSTAT)
    ylab <- "MAGMA Z-score"
  }else if( z.or.p == "p" ){
    mag$Y <- -log10(mag$P)
    ylab <- "-log10 P value"
  }
  
  # Determine P value thresholds
  zbonf <- p_to_z( 0.05 / NROW(mag) )
  znom  <- p_to_z( 0.05 )
  
  # Re-scale positions to Mbp
  message2("Re-scale positions to Mbp")
  mag$START <- mag$START/1e6
  mag$STOP  <- mag$STOP/1e6
  
  # Read in gene locations
  message2("Read in gene locations")
  gene_file <- "/projects/0/prjs0817/projects/pops/data/gene_locations.tsv"
  genes     <- fread(gene_file)
  names(genes) <- c( "GENE", "CHR", "START", "STOP", "TSS", "STRAND", "NAME" )
  
  # Add HGNC gene names to MAGMA
  message2("Add gene locations to MAGMA")
  mag$NAME <- genes$NAME[ match( mag$GENE, genes$GENE ) ]
  
  # Make an output directory to store plots
  mag_plot_dir <- file.path( maindir, "plots", "magma" )
  dir.create( path=mag_plot_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Loop through loci
  message2("Loop through loci, making a plot for each")
  for( i in seq_along( peaks$hit ) ){
    
    # Subset MAGMA results
    # message2("Starting locus: ", i, "/", NROW(peaks) )
    xmin <- peaks$lo[i]/1e6
    xmax <- peaks$hi[i]/1e6
    locus <- mag[ mag$CHR == peaks$chr[i] & 
                    mag$STOP  > xmin &
                    mag$START < xmax , ]
    
    # If there are no genes in the locus, expand the locus boundaries by +/- 2Mb
    if( NROW(locus) == 0 ){
      locus <- mag[ mag$CHR == peaks$chr[i] & 
                      mag$STOP  > xmin - 2 &
                      mag$START < xmax + 2 , ]
    }
    
    # Initialize the output file
    mag_plot_file <- file.path( mag_plot_dir, paste0( "region_", i, "_", 
                                                      peaks$snp[i], ".jpg" ) )
    jpeg( filename=mag_plot_file, width=480*4, height=270*4, res=75*4 )
    
    # Set up the plot
    ymax <- max(locus$Y) *1.08
    xlab <- paste0( "Chromosome ", peaks$chr[i], " (Mb)")
    par( mar=c( 4, 3.5, 0.5, 1.5 ) )
    plot( x=locus$START, y=locus$Y, xlim=c(xmin,xmax), ylim=c(0,ymax),
          xlab=xlab, ylab="", las=1, type="n" )
    title( ylab=ylab, line=2.3 )
    
    # Add horizontal lines
    abline( h=0,     lwd=2, col="grey70" )
    abline( h=zbonf, lwd=2, col="grey70", lty=2 )
    abline( h=znom,  lwd=2, col="grey70", lty=2 )
    
    # Add bars for each gene
    for( j in seq_len( NROW(locus) ) ){
      lines( x = c( locus$START[j], locus$STOP[j] ),
             y = c( locus$Y[j],     locus$Y[j] ),
             lwd=3, col="#70AD47" )
    }
    
    # Add gene names
    # Take special pains to include genes that don't lie entirely within the plot
    in_window <- locus$START > xmin & locus$STOP < xmax
    out_left  <- locus$START < xmin
    out_right <- locus$START > xmin & locus$STOP > xmax
    x_buffer <- ( xmax - xmin ) * 0.03
    xs <- ( locus$START + locus$STOP ) / 2
    ys <- locus$Y + ymax*0.05
    if( sum(in_window) > 0 ){
      text( x=xs[in_window], y=ys[in_window], labels=locus$NAME[in_window], 
            adj=c(0.5,0.5), col="#70AD47" )
    }
    if( sum(out_left) > 0 ){
      text( x=xmin - x_buffer, y=ys[out_left],  labels=locus$NAME[out_left],  
            adj=c(0,0.5), col="#70AD47" )
    }
    if( sum(out_right) > 0 ){
      text( x=xmax + x_buffer, y=ys[out_right], labels=locus$NAME[out_right], 
            adj=c(1,0.5), col="#70AD47" )
    }
    dev.off()
  }
}


#-------------------------------------------------------------------------------
#   run_pops
#-------------------------------------------------------------------------------

run_pops <- function(maindir){
  
  # Create job name and log file name
  jobname <- paste0( "p", sample( x=1:999, size=1 ) )
  logfile <- file.path( maindir, "logs/pops.log" )
  
  # Set up the command
  bash_script <- "/projects/0/prjs0817/repos/brett/f_run_pops.sh"
  cmd <- paste( "sbatch", 
                "-J", jobname,
                "-o", logfile,
                "-e", logfile,
                bash_script, 
                maindir )
  
  # Run
  system(cmd)
  ensure_finished_jobs(jobname)
}


#-------------------------------------------------------------------------------
#   pops_plots
#-------------------------------------------------------------------------------

pops_plots <- function( maindir, loci_dir ){
  
  # Read in peaks
  message2("Read in peaks")
  library(data.table)
  loci_file <- file.path( loci_dir, "loci_cs.tsv" )
  peaks <- fread(loci_file)
  pattern <- "^chr[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_(.*)$"
  peaks$snp <- sub( pattern=pattern, replacement="\\1", x=peaks$hit )
  # peaks <- peaks[ order( peaks$chr, peaks$bp ) , ]
  
  # Read in gene locations
  message2("Read in gene locations")
  gene_file <- "/projects/0/prjs0817/projects/pops/data/gene_locations.tsv"
  genes     <- fread(gene_file)
  
  # Read in POPS results, calculate P values
  message2("Read in POPS results, calculate P values")
  pops_file <- file.path( maindir, "pops.preds" )
  pops <- fread(pops_file)
  names(pops)[ names(pops) == "PoPS_Score" ] <- "pops"
  pops$p <- z_to_p( z=pops$pops )
  
  # Add gene names and locations to POPS results
  pops$gene  <- genes$NAME[  match( pops$ENSGID, genes$ENSGID ) ]
  pops$CHR   <- genes$CHR[   match( pops$ENSGID, genes$ENSGID ) ]
  pops$START <- genes$START[ match( pops$ENSGID, genes$ENSGID ) ]
  pops$STOP  <- genes$END[   match( pops$ENSGID, genes$ENSGID ) ]
  
  # Use z-score or P value as the plot's y-axis, depending on choice
  z.or.p <- "z"
  if( z.or.p == "z" ){
    pops$Y <- pops$pops
    ylab <- "PoP Score"
  }else if( z.or.p == "p" ){
    pops$Y <- -log10(pops$p)
    ylab <- "-log10 PoPS P value"
  }
  
  # Establish the maximum and minimum y-axis values
  # ymin <- min(pops$Y)
  # ymax <- max(pops$Y) * 1.05
  
  # Determine the "significance" threshold
  y_idx <- order( pops$Y, decreasing=TRUE )
  ysig1 <- pops$Y[ y_idx[ round( NROW(pops) * 0.01 ) ] ]
  ysig5 <- pops$Y[ y_idx[ round( NROW(pops) * 0.05 ) ] ]
  
  # Re-scale positions to Mbp
  pops$START <- pops$START/1e6
  pops$STOP  <- pops$STOP/1e6
  
  # Make an output directory to store plots
  pops_plot_dir <- file.path( maindir, "plots", "pops" )
  dir.create( path=pops_plot_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Loop through loci
  for( i in seq_along(peaks$hit) ){
    
    # Subset PoPS results
    # message2("Starting locus: ", i, "/", NROW(peaks) )
    xmin <- peaks$lo[i]/1e6
    xmax <- peaks$hi[i]/1e6
    idx  <- pops$CHR == peaks$chr[i] & 
      pops$STOP  > xmin &
      pops$START < xmax
    locus <- pops[ idx , ]
    
    # If there are no genes in the locus, expand the locus boundaries by +/- 2Mb
    if( NROW(locus) == 0 ){
      idx  <- pops$CHR == peaks$chr[i] & 
        pops$STOP  > xmin - 2 &
        pops$START < xmax + 2
      locus <- pops[ idx , ]
    }
    
    # Initialize the output file
    pops_plot_file <- file.path( pops_plot_dir, paste0( "region_", i, "_", 
                                                        peaks$snp[i], ".jpg" ) )
    jpeg( filename=pops_plot_file, width=480*4, height=270*4, res=75*4 )
    
    # Set up the plot
    ymin <- min( c( 0, locus$Y ) )
    ymax <- max( c( ysig1, locus$Y ) ) * 1.08
    xlab <- paste0( "Chromosome ", peaks$chr[i], " (Mb)")
    par( mar=c( 4, 3.5, 0.5, 1.5 ) )
    plot( x=locus$START, y=locus$Y, xlim=c(xmin,xmax), ylim=c(ymin,ymax),
          xlab=xlab, ylab="", las=1, type="n" )
    title( ylab=ylab, line=2.3 )
    
    # Add horizontal lines
    abline( h=0,    lwd=2, col="grey70" )
    abline( h=ysig1, lwd=2, col="grey70", lty=2 )
    abline( h=ysig5, lwd=2, col="grey70", lty=2 )
    
    # Add bars for each gene
    for( j in seq_len( NROW(locus) ) ){
      lines( x = c( locus$START[j], locus$STOP[j] ),
             y = c( locus$Y[j],     locus$Y[j] ),
             lwd=3, col="#70AD47" )
    }
    
    # Add gene names
    # Take special pains to include genes that don't lie entirely within the plot
    in_window <- locus$START > xmin & locus$STOP < xmax
    out_left  <- locus$START < xmin
    out_right <- locus$START > xmin & locus$STOP > xmax
    x_buffer <- ( xmax - xmin ) * 0.03
    xs <- ( locus$START + locus$STOP ) / 2
    ys <- locus$Y + ymax*0.05
    if( sum(in_window) > 0 ){
      text( x=xs[in_window], y=ys[in_window], labels=locus$gene[in_window], 
            adj=c(0.5,0.5), col="#70AD47" )
    }
    if( sum(out_left) > 0 ){
      text( x=xmin - x_buffer, y=ys[out_left],  labels=locus$gene[out_left],  
            adj=c(0,0.5), col="#70AD47" )
    }
    if( sum(out_right) > 0 ){
      text( x=xmax + x_buffer, y=ys[out_right], labels=locus$gene[out_right], 
            adj=c(1,0.5), col="#70AD47" )
    }
    dev.off()
  }
}


#-------------------------------------------------------------------------------
#   wrapper
#-------------------------------------------------------------------------------

brett <- function( maindir    = "/projects/0/prjs0817/projects/analyses/pd",
                   ld_panel   = "hrc",
                   population = "eur",
                   gw_file    = "/projects/0/prjs0817/projects/pops/analyses/pd/meta5_raw.tab.gz",
                   loci_dir   = NULL,
                   check_args = TRUE ){
  
  
  #-------------------------------------------------------------------------------
  #   Input descriptions
  #-------------------------------------------------------------------------------
  
  #   maindir:    Main directory in which to store results
  #   ld_panel:   Which LD reference panel should be used? Options are either: 
  #               "hrc" or "g1000".
  #   population: Which populations should be used? Options are either: 
  #               "eur" or "eas".
  #   gw_file:    GWAS file name
  #   loci_dir:   Directory containing a file showing boundaries for each locus
  #   check_args: Logical. Check whether arguments are valid?
  
  
  #-------------------------------------------------------------------------------
  #   Print inputs
  #-------------------------------------------------------------------------------
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/brett/z_brett.R")
  message2("Load libraries and sources")
  suppressPackageStartupMessages( library(data.table) )
  
  # Print inputs
  message_header("Print inputs")
  message2( "Main directory: ",     maindir )
  message2( "LD reference panel: ", ld_panel )
  message2( "Population: ",         population )
  message2( "GWAS file: ",          gw_file )
  message2( "Locus directory: ",    loci_dir )
  
  
  #-------------------------------------------------------------------------------
  #   Check arguments
  #-------------------------------------------------------------------------------
  
  if(check_args){
    message_header("Check arguments")
    message2("Checking arguments")
    check_arguments( ld_panel   = ld_panel,
                     population = population,
                     gw_file    = gw_file )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Assign output file names
  #-------------------------------------------------------------------------------
  
  # Make output and log directories
  logdir <- file.path( maindir, "logs" )
  dir.create( path=logdir, showWarnings=FALSE, recursive=TRUE )
  
  # Assign output file names
  clean_gw_file         <- file.path( maindir, "gwas_pvalues.tsv" )
  snp_loc_file          <- file.path( maindir, "snp_locations.tsv" )
  raw_snp_map_file      <- file.path( maindir, "snps_mapped_to_genes.genes.annot" )
  clean_snp_map_file    <- file.path( maindir, "snps_mapped_to_genes_filtered.genes.annot" )
  mag_ss_collated_file  <- file.path( maindir, "magma.genes.out" )
  mag_cov_collated_file <- file.path( maindir, "magma.genes.raw" )
  mag_plots_file        <- file.path( maindir, "magma_plots.pdf" )
  pops_coef_file        <- file.path( maindir, "pops.coefs" )
  pops_marg_file        <- file.path( maindir, "pops.marginals" )
  pops_pred_file        <- file.path( maindir, "pops.preds" )
  pops_plots_file       <- file.path( maindir, "pops_plots.pdf" )
  html_file             <- file.path( maindir, "report.html" )
  
  # Assign output file names that are dependent on the reference panel
  if( ld_panel == "hrc" ){
    mag_sumstats_files    <- file.path( maindir, "magma", paste0( "chr", 1:22, ".genes.out" ) )
    mag_covar_files       <- file.path( maindir, "magma", paste0( "chr", 1:22, ".genes.raw" ) )
  }else if( ld_panel == "g1000" ){
    mag_sumstats_files    <- mag_ss_collated_file
    mag_covar_files       <- mag_cov_collated_file
  }else{
    stop("ld_panel must be either 'hrc' or 'g1000'")
  }
  
  
  #-------------------------------------------------------------------------------
  #   Format the GWAS
  #   Write the SNP location file
  #   Write the peaks file
  #-------------------------------------------------------------------------------
  
  message_header("Format the GWAS, write the SNP location and peaks files")
  if( all( file.exists( clean_gw_file, 
                        snp_loc_file ) ) ){
    message2("Output files exist, skipping")
    
  }else{
    format_gwas_and_snp_loc_files( maindir    = maindir,
                                   ld_panel   = ld_panel,
                                   population = population,
                                   gw_file    = gw_file )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Map SNPs to genes
  #-------------------------------------------------------------------------------
  
  message_header("Map SNPs to genes")
  if( file.exists(raw_snp_map_file) ){
    message2("Output file exists, skipping")
  }else{
    map_snps_to_genes( maindir = maindir )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Remove genes without enough SNPs
  #-------------------------------------------------------------------------------
  
  message_header("Remove genes without enough SNPs")
  if( file.exists(clean_snp_map_file) ){
    message2("Output file exists, skipping")
  }else{
    message2("Removing genes without enough SNPs") 
    rm_genes_without_enough_snps( maindir = maindir )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Run MAGMA
  #-------------------------------------------------------------------------------
  
  message_header("Run MAGMA")
  if( all( file.exists( mag_sumstats_files, 
                        mag_covar_files ) ) ){
    message2("Output files exist, skipping")
  }else{
    run_magma( maindir    = maindir,
               ld_panel   = ld_panel,
               population = population )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Collate MAGMA results
  #-------------------------------------------------------------------------------
  
  message_header("Collate MAGMA results")
  if( all( file.exists( mag_ss_collated_file,
                        mag_cov_collated_file ) ) ){
    message2("Output file exists, skipping")
  }else{
    message2("Collating MAGMA results")
    collate_magma( maindir = maindir )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Make MAGMA plots
  #-------------------------------------------------------------------------------
  
  message_header("Make MAGMA plots")
  if( file.exists(mag_plots_file) ){
    message2("Output file exists, skipping")
  }else{
    message2("Making MAGMA plots")
    magma_plots( maindir  = maindir,
                 loci_dir = loci_dir )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Run POPS
  #-------------------------------------------------------------------------------
  
  message_header("Run POPS")
  if( all( file.exists( pops_coef_file,
                        pops_marg_file, 
                        pops_pred_file ) ) ){
    message2("Output files exist, skipping")
  }else{
    message2("Running POPS")
    run_pops( maindir = maindir )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Make POPS plots
  #-------------------------------------------------------------------------------
  
  message_header("Make POPS plots")
  if( file.exists(pops_plots_file) ){
    message2("Output file exists, skipping")
  }else{
    message2("Making POPS plots")
    pops_plots( maindir  = maindir,
                loci_dir = loci_dir )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Render an HTML report
  #-------------------------------------------------------------------------------
  
  message_header("Render an HTML report")
  if( file.exists(html_file) ){
    message2("Output file exists, skipping")
  }else{
    message2("Rendering an HTML report")
    library(rmarkdown)
    rmd_file <- file.path( maindir, "report.Rmd" )
    file.copy( from = "/projects/0/prjs0817/repos/brett/g_brett_template.Rmd",
               to   = rmd_file, overwrite=TRUE )
    args <- list( maindir    = maindir,
                  ld_panel   = ld_panel,
                  gw_file    = gw_file,
                  loci_dir   = loci_dir,
                  check_args = check_args )
    render( input       = rmd_file, 
            params      = args, 
            output_file = html_file )
  }
  message_header("Done")
}


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------



