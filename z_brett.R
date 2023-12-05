
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
# z_to_p:               Convert a z-score into a P value


# ensure_finished_jobs: Make sure all qsub jobs finish running before proceeding
ensure_finished_jobs <- function(identifier){
  
  message( date(), "   Make sure all jobs finish running before proceeding" )
  external.call <- paste0( "squeue | grep ", identifier, " | wc -l" )
  running.jobs  <- as.numeric( system( external.call, intern=TRUE ) ) 
  while( running.jobs > 0){
    Sys.sleep(10)
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

check_arguments <- function( ld.panel   = NULL, 
                             population = NULL,
                             gw.file    = NULL,
                             chr.bp.col = NULL,
                             chr.col    = NULL,
                             bp.col     = NULL,
                             a1.col     = NULL,
                             a2.col     = NULL,
                             p.col      = NULL,
                             eaf.col    = NULL,
                             n1.col     = NULL,
                             n0.col     = NULL,
                             n.col      = NULL,
                             n          = NULL,
                             z.or.p     = NULL ){
  
  # ld.panel must be either "hrc" or "g1000"
  if( ld.panel != "hrc" & ld.panel != "g1000" ) stop("ld.panel must be either 'hrc' or 'g1000'")
  
  # population must be either "eur" or "eas"
  if( population != "eur" & population != "eas" ) stop("population must be either 'eur' or 'eas'")
  
  # Does the GWAS file exist?
  if( !file.exists(gw.file) )  stop("GWAS file does not exist")
  
  # Either chr.bp.col or (chr.col + bp.col) must be specified
  if( is.null(chr.bp.col) ){
    if( is.null(chr.col) | is.null(bp.col) ){
      stop("chr.bp.col is not specified so both chr.col and bp.col must be specified")
    }
  }else{
    if( !is.null(chr.col) | !is.null(bp.col) ){
      stop("chr.bp.col is specified so chr.col and bp.col must not be specified")
    }
  }
  
  # One of the following must be specified:
  # n or n.col or (n1.col + n0.col)
  if( !is.null(n) ){
    if( !is.null(n.col) | !is.null(n1.col) | !is.null(n0.col) ){
      stop("n is specified so n.col, n1.col, and n0.col must not be specified")
    }
  }else if( !is.null(n.col) ){
    if( !is.null(n1.col) | !is.null(n0.col) ){
      stop("n.col is specified so n1.col and n0.col must not be specified")
    }
  }else{
    if( is.null(n1.col) | is.null(n0.col) ){
      stop("Neither n nor n.col are specified so both n1.col and n0.col must be specified")
    }
  }
  
  # z.or.p must be either "z" or "p"
  if( z.or.p != "z" & z.or.p != "p" ) stop("z.or.p must be either 'z' or 'p'")
  
  # Check that GWAS file column names exist
  col.names <- c( chr.bp.col, chr.col, bp.col, a1.col, a2.col, 
                  p.col, eaf.col, n1.col, n0.col, n.col )
  library(data.table)
  gw <- fread( file=gw.file, nrows=100 )
  bad.col.names <- setdiff( col.names, names(gw) )
  if( length(bad.col.names) > 0 ){
    stop("The following specified column names do not exist: ", 
         paste( bad.col.names, collapse=", " ) )
  }
  
  # Position must be a positive integer
  if( !is.null(bp.col) ){
    if( !is.integer( gw[[bp.col]] ) ) stop("Positions must be integers")
    if( any( gw[[bp.col]] ) < 1  )  stop("Positions must be > 0")
  }
  
  # Alleles must be characters
  if( !is.character( gw[[a1.col]] ) ) stop("Effect alleles must be characters")
  if( !is.character( gw[[a2.col]] ) ) stop("Non-effect alleles must be characters")
  
  # P value must be >= 0 and <= 1
  if( !is.numeric( gw[[p.col]] ) ) stop("P values must be numeric")
  if( any(   gw[[p.col]] ) < 0 )   stop("P values must be >= 0")
  if( any(   gw[[p.col]] ) > 1 )   stop("P values must be <= 1")
  
  # Effect allele frequency must be >= 0 and <= 1
  if( eaf.col %in% names(gw) ){
    if( !is.numeric( gw[[eaf.col]] ) ) stop("Effect allele frequencies must be numeric")
    if( any(   gw[[eaf.col]] ) < 0 )   stop("Effect allele frequencies must be >= 0")
    if( any(   gw[[eaf.col]] ) > 1 )   stop("Effect allele frequencies must be <= 1")
  }
  
  # Case counts must be positive integers
  if( !is.null(n1.col) ){
    if( !is.integer( gw[[n1.col]] ) ) stop("Case counts must be integers")
    if( any( gw[[n1.col]] < 1 ) )     stop("Case counts must be > 0")
  }
  
  # Control counts must be positive integers
  if( !is.null(n0.col) ){
    if( !is.integer( gw[[n0.col]] ) ) stop("Control counts must be integers")
    if( any( gw[[n0.col]] < 1 ) )     stop("Control counts must be > 0")
  }
  
  # Effective sample sizes must be positive numbers
  if( !is.null(n.col) ){
    if( !is.numeric( gw[[n.col]] ) ) stop("Effective sample sizes must be numbers")
    if( any( gw[[n.col]] <= 0 ) )    stop("Effective sample sizes must be > 0")
  }
  
  # Effective sample size must be a positive number
  if( !is.null(n) ){
    if( !is.numeric(n) ) stop("Effective sample size must be a number")
    if( n <= 0 )         stop("Effective sample size must be > 0")
  }
}


#-------------------------------------------------------------------------------
#   format_gwas_and_snp_loc_files
#-------------------------------------------------------------------------------

format_gwas_and_snp_loc_files <- function( maindir    = "/projects/0/prjs0817/projects/pops/analyses/pd",
                                           ld.panel   = "hrc",
                                           population = "eur",
                                           gw.file    = "/projects/0/prjs0817/projects/pops/analyses/pd/meta5_raw.tab.gz",
                                           chr.bp.col = "SNP",
                                           chr.col    = NULL,
                                           bp.col     = NULL,
                                           a1.col     = "A1",
                                           a2.col     = "A2",
                                           p.col      = "p",
                                           eaf.col    = "freq",
                                           n1.col     = "N_cases",
                                           n0.col     = "N_controls",
                                           n.col      = NULL,
                                           n          = NULL ){
  
  #-------------------------------------------------------------------------------
  #   Input descriptions
  #-------------------------------------------------------------------------------
  
  #   maindir:    Main directory in which to store results
  #   ld.panel:   Which LD reference panel should be used? Options are either: 
  #               "hrc" or "g1000".
  #   population: Which populations should be used? Options are either: 
  #               "eur" or "eas".
  #   gw.file:    GWAS file name
  #   chr.bp.col: Optional. The name of a GWAS column containing chromosome and bp 
  #               information separated by a punctuation character. Must be
  #               specified if chr.col and bp.col are not specified.
  #   chr.col:    Optional. The name of a GWAS column containing chromosome
  #               information. Must be specified if chr.bp.col is not specified.
  #   bp.col:     Optional. The name of a GWAS column containing position (bp)
  #               information. Must be specified if chr.bp.col is not specified.
  #   a1.col:     The name of a GWAS column containing the effect (alt) allele
  #   a2.col:     The name of a GWAS column containing the non-effect (ref) allele
  #   p.col:      The name of a GWAS column containing the P value
  #   eaf.col     Optional. The name of a GWAS column containing the frequency of
  #               the effect (alt) allele. If not specified, all palindromic SNPs
  #               will be removed. Otherwise palindromic SNPs with similar 
  #               frequencies in GWAS and HRC will be preserved.
  #   n1.col:     The name of a GWAS column containing the per-SNP number of 
  #               cases. Must be specified with n0.col, or must specify n.col 
  #               or n.
  #   n0.col:     The name of a GWAS column containing the per-SNP number of 
  #               controls. Must be specified with n1.col, or must specify n.col 
  #               or n.
  #   n.col:      The name of a GWAS column containing the per-SNP effective
  #               sample size. If not specified, must specify both n1.col and 
  #               n0.col, or n.
  #   n:          If per-SNP sample size information is not available, this
  #               specifies the study-wide effective sample size. If not 
  #               specified, must specify both n1.col and n0.col, or n.col.
  
  
  #-------------------------------------------------------------------------------
  #   Read in GWAS and HRC, format columns, subset to shared SNPs
  #-------------------------------------------------------------------------------
  
  # Load libraries and sources
  library(data.table)
  
  # Read in reference panel SNPs
  if( ld.panel == "hrc" ){                                  ### HRC
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
  }else if( ld.panel == "g1000" ){                          ### 1000G
    message2("Read in 1000 Genomes SNPs with EUR MAC >= 10")
    hrc <- fread("/projects/0/prjs0817/projects/pops/data/g1000_eur_snps_mac_ge_10.tsv")
  }else{
    stop("ld.panel must be either 'hrc' or 'g1000'")
  }
  
  # Read in GWAS
  message2("Read in GWAS")
  gw <- fread(gw.file)
  
  # Re-name GWAS columns
  message2("Re-name GWAS columns")
  names(gw)[ names(gw) == chr.bp.col ] <- "chr.bp"
  names(gw)[ names(gw) == chr.col    ] <- "chr"
  names(gw)[ names(gw) == bp.col     ] <- "bp"
  names(gw)[ names(gw) == a1.col     ] <- "a1"
  names(gw)[ names(gw) == a2.col     ] <- "a2"
  names(gw)[ names(gw) == p.col      ] <- "P"
  names(gw)[ names(gw) == eaf.col    ] <- "eaf"
  names(gw)[ names(gw) == n1.col     ] <- "n1"
  names(gw)[ names(gw) == n0.col     ] <- "n0"
  names(gw)[ names(gw) == n.col      ] <- "N"
  
  # Make columns for chromosome and position
  if( "chr.bp" %in% names(gw) ){
    message2("Make columns for chromosome and position")
    gw$chr <- as.integer( sub( pattern     = "^chr([[:alnum:]]+)[[:punct:]]([[:digit:]]+)$", 
                               replacement = "\\1", 
                               x=gw$chr.bp ) )
    gw$bp  <- as.integer( sub( pattern     = "^chr([[:alnum:]]+)[[:punct:]]([[:digit:]]+)$", 
                               replacement = "\\2", 
                               x=gw$chr.bp ) )
  }
  
  # Create a per-SNP effective sample size column
  #   If this column already exists, do nothing
  #   If an overall study N is provided, use it for all SNPs
  #   Otherwise, compute the effective N from the number of cases and controls
  if( "N" %in% names(gw) ){
    message2("An effective sample size column has been provided and will be used")
  }else if( !is.null(n) ){
    message2("A study-wide effective sample size has been provided and will be applied to all SNPs")
    gw$N <- n
  }else if( "n1" %in% names(gw) & "n0" %in% names(gw) ){
    message2("Columns for number of cases and controls have been provided for each SNP, computing the effective sample size")
    gw$N <- n_eff( gw$n1, gw$n0 )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Harmonize GWAS and reference panel alleles: without allele frequency information
  #-------------------------------------------------------------------------------
  
  # Subset GWAS and reference panel to shared SNPs based on chromosome and position
  message2("Subset GWAS and reference panel to shared SNPs based on chromosome and position")
  cpab_gw  <- paste( gw$chr,  gw$bp,  
                     ifelse( gw$a1   < gw$a2,   gw$a1,   gw$a2  ), 
                     ifelse( gw$a1    < gw$a2,   gw$a2,   gw$a1  ),  sep="_" )
  cpab_hrc <- paste( hrc$chr, hrc$bp, 
                     ifelse( hrc$alt < hrc$ref, hrc$alt, hrc$ref ), 
                     ifelse( hrc$alt < hrc$ref, hrc$ref, hrc$alt ), sep="_" )
  cpab_both <- intersect( cpab_hrc, cpab_gw )
  hrc2 <- hrc[ match( cpab_both, cpab_hrc ) , ]
  gw2  <- gw[  match( cpab_both, cpab_gw  ) , ]
  message2( "Of the ", NROW(gw), " GWAS SNPs, ", NROW(gw2), 
           " (", round( 100*NROW(gw2)/NROW(gw), 2 ), "%) were found in the reference panel" )
  message2( "Of the ", NROW(hrc), " reference panel SNPs, ", NROW(hrc2), 
           " (", round( 100*NROW(hrc2)/NROW(hrc), 2 ), "%) were found in the GWAS" )
  
  # Find palindromic SNPs and SNPs with alleles that are flipped in 
  # the reference panel v. GWAS ('discordant')
  message2("Find palindromic SNPs and SNPs with alleles that are flipped in ",
           "the reference panel v. GWAS ('discordant')")
  pal     <- ( hrc2$alt=="A" & hrc2$ref=="T" ) | 
             ( hrc2$alt=="T" & hrc2$ref=="A" ) | 
             ( hrc2$alt=="C" & hrc2$ref=="G" ) | 
             ( hrc2$alt=="G" & hrc2$ref=="C" )
  discord <- hrc2$alt != gw2$a1
  message2( sum(pal),     "/", NROW(hrc2), " (", round( 100 * sum(pal)     / NROW(hrc2), 2 ), "%) SNPs are palindromic" )
  message2( sum(discord), "/", NROW(hrc2), " (", round( 100 * sum(discord) / NROW(hrc2), 2 ), "%) SNPs have discordant alleles" )
  
  # For discordant non-palindromic SNPs: flip GWAS alleles
  message2("For discordant non-palindromic SNPs: flip GWAS alleles")
  disc_nonpal <- discord & !pal
  original_gwas_ref <- gw2$a2
  original_gwas_alt <- gw2$a1
  gw2$a2[disc_nonpal] <- original_gwas_alt[disc_nonpal]
  gw2$a1[disc_nonpal] <- original_gwas_ref[disc_nonpal]
  if( "eaf" %in% names(gw2) ){
    gw2$eaf[disc_nonpal] <- 1  - gw2$eaf[disc_nonpal]
  }
  message2( sum(disc_nonpal), "/", sum(discord), 
           " (", round( 100 * sum(disc_nonpal) / sum(discord), 2 ), 
           "%) discordant SNPs were non-palindromic, flipping alleles" )
  
  
  #-------------------------------------------------------------------------------
  #   Harmonize GWAS and HRC alleles: with allele frequency information
  #-------------------------------------------------------------------------------
  
  # If allele frequency data is available for the GWAS
  if( "eaf" %in% names(gw2) ){
    
    # For discordant palindromic SNPs with compatible AFs: flip GWAS alleles
    diff_af_disc    <- abs( ( 1 - gw2$eaf ) - hrc2$af ) > 0.2 | hrc2$af > 0.4
    disc_pal_compat <- discord & pal & !diff_af_disc
    gw2$a2[disc_pal_compat]  <- original_gwas_alt[disc_pal_compat]
    gw2$a1[disc_pal_compat]  <- original_gwas_ref[disc_pal_compat]
    gw2$eaf[disc_pal_compat] <- 1  - gw2$eaf[disc_pal_compat]
    message2( sum(disc_pal_compat), "/", sum(discord), 
             " (", round( 100 * sum(disc_pal_compat) / sum(discord), 2 ), 
             "%) discordant SNPs were palindromic but with AFs that were clearly ",
             "compatible with the reference panel, flipping alleles" )
    
    # Flag discordant palindromic SNPs with incompatible AFs for removal
    disc_pal_incompat <- discord & pal & diff_af_disc
    message2( sum(disc_pal_incompat), "/", sum(discord), 
              " (", round( 100 * sum(disc_pal_incompat) / sum(discord), 2 ), 
             "%) discordant SNPs were palindromic and had AFs that were not ",
             "clearly compatible with the reference panel, flagging for removal" )
    
    # For concordant palindromic SNPs with clearly incompatible AFs: flip GWAS alleles
    diff_af_conc <- abs( gw2$eaf - hrc2$af ) > 0.2
    common_af    <- hrc2$af > 0.4
    conc_pal_incompat <- !discord & pal & diff_af_conc & !common_af
    gw2$eaf[conc_pal_incompat] <- 1  - gw2$eaf[conc_pal_incompat]
    message2( sum(conc_pal_incompat), "/", sum( !discord & pal ), 
             " (", round( 100 * sum(conc_pal_incompat) / sum( !discord & pal ), 2 ), 
             "%) concordant palindromic SNPs had AFs that were clearly ",
             "different from the reference panel, flipped alleles" )
    
    # Report the number of concordant palindromic SNPs with clearly compatible AFs
    conc_pal_compat <- !discord & pal & !diff_af_conc & !common_af
    message2( sum(conc_pal_compat), "/", sum( !discord & pal ), 
             " (", round( 100 * sum(conc_pal_compat) / sum( !discord & pal ), 2 ), 
             "%) concordant palindromic SNPs had AFs that were clearly similar ",
             "to the reference panel, no action" )
    
    # Flag concordant palindromic SNPs with HRC MAF > 40% for removal
    conc_pal_ambig <- !discord & pal & common_af
    message2( sum(conc_pal_ambig), "/", sum( !discord & pal ), 
             " (", round( 100 * sum(conc_pal_ambig) / sum( !discord & pal ), 2 ), 
             "%) concordant palindromic SNPs had reference panel MAF > 40%, ",
             "flagging for removal" )
    
    # Remove flagged SNPs
    message2("Remove flagged SNPs")
    bad_snps <- disc_pal_incompat | conc_pal_ambig
    gw3  <- gw2[  !bad_snps , ]
    hrc3 <- hrc2[ !bad_snps , ]
    
  }else{
    
    # Remove palindromic SNPs
    message2("Remove palindromic SNPs")
    gw3  <- gw2[  !pal , ]
    hrc3 <- hrc2[ !pal , ]
  }
  
  
  #-------------------------------------------------------------------------------
  #   Wrap up harmonization
  #-------------------------------------------------------------------------------
  
  # Report the change in number of SNPs
  message2( "After harmonizing GWAS and reference panel SNPs, ", 
            NROW(gw3), "/", NROW(gw2), " (", round( 100 * NROW(gw3) / NROW(gw2), 2 ), 
            "%) remain" )
  
  # Check that CPRA is 100% identical now
  cpra_gw  <- paste( gw3$chr,  gw3$bp,  gw3$a2,   gw3$a1,   sep="_" )
  cpra_hrc <- paste( hrc3$chr, hrc3$bp, hrc3$ref, hrc3$alt, sep="_" )
  if( all( cpra_gw == cpra_hrc ) ){
    message2("CPRA is now identical for all GWAS and reference panel SNPs")
  }else{
    diff_cpra <- head( which( cpra_gw != cpra_hrc ) )
    stop( paste( "Error: not all CPRA are identical for GWAS and reference panel",
                 "SNPs. Here are (up to) the first 6:", diff_cpra, collapse=" " ) )
  }
  
  # Replace GWAS SNP names with reference panel names
  message2("Replace GWAS SNP names with reference panel names")
  gw3$SNP <- hrc3$snp
  
  
  #-------------------------------------------------------------------------------
  #   Find peaks
  #-------------------------------------------------------------------------------
  
  # Subset to significant SNPs
  message2("Subset to significant SNPs")
  sig_ss <- gw3[ gw3$P < 5e-8 , ]
  
  # Loop through loci and remove until none remain
  message2("Loop through loci and remove until none remain")
  peaks0 <- list()
  n_sig_snps <- NROW(sig_ss)
  while( n_sig_snps > 0 ){
    
    # Find top remaining SNP
    top_idx <- which( sig_ss$P == min(sig_ss$P) )[1]
    top_snp <- sig_ss$SNP[top_idx]
    top_chr <- sig_ss$chr[top_idx]
    top_bp  <- sig_ss$bp[top_idx]
    top_p   <- sig_ss$P[top_idx]
    # message2( "Analyzing locus ", length(peaks0)+1, ": chr", top_chr, ":", top_bp )
    
    # Subset to +/- 2Mb around it
    boundaries_right <- top_bp + seq( 0, 2e6, 2e5 )
    boundaries_left  <- top_bp - seq( 0, 2e6, 2e5 )
    top_sig_ss <- sig_ss[ sig_ss$chr == top_chr & 
                            sig_ss$bp < tail( boundaries_right, 1 ) &
                            sig_ss$bp > tail( boundaries_left, 1 ) , ]
    
    # Check a series of 200kb windows around the hit for significant SNPs
    any_sig_in_bin_right <- any_sig_in_bin_left <- list()
    for( i in seq_len( length(boundaries_right)  -  1 ) ){
      any_sig_in_bin_right[[i]] <- any( top_sig_ss$chr == top_chr & 
                                          top_sig_ss$bp >= boundaries_right[i] & 
                                          top_sig_ss$bp < boundaries_right[i+1] )
      any_sig_in_bin_left[[i]]  <- any( top_sig_ss$chr == top_chr & 
                                          top_sig_ss$bp <= boundaries_left[i] & 
                                          top_sig_ss$bp > boundaries_left[i+1] )
    }
    any_sig_in_bin_right <- unlist(any_sig_in_bin_right)
    any_sig_in_bin_left  <- unlist(any_sig_in_bin_left)
    
    # Find the last window to the right and left that have any P < 5e-8
    bin_idx_rightmost <- tail( which(any_sig_in_bin_right), 1 )
    bin_idx_leftmost  <- tail( which(any_sig_in_bin_left),  1 )
    
    # Take the most distant SNPs with P < 5e-8 as the boundaries
    bin_rightmost <- top_sig_ss[ top_sig_ss$chr == top_chr &
                                   top_sig_ss$bp >= boundaries_right[ bin_idx_rightmost ] &
                                   top_sig_ss$bp <  boundaries_right[ bin_idx_rightmost + 1 ] , ]
    bin_leftmost  <- top_sig_ss[ top_sig_ss$chr == top_chr &
                                   top_sig_ss$bp <= boundaries_left[ bin_idx_leftmost ] &
                                   top_sig_ss$bp >  boundaries_left[ bin_idx_leftmost + 1 ] , ]
    bp_rightmost <- tail( bin_rightmost$bp, 1 )
    bp_leftmost  <- head( bin_leftmost$bp,  1 )
    
    # Record boundaries
    peaks0[[top_snp]] <- data.frame( snp=top_snp, chr=top_chr, bp=top_bp, 
                                     lo=bp_leftmost, hi=bp_rightmost, p=top_p )
    
    # Remove this region from the GWAS sumstats
    in_this_locus <- sig_ss$chr == top_chr & 
                     sig_ss$bp <= bp_rightmost &
                     sig_ss$bp >= bp_leftmost
    sig_ss <- sig_ss[ !in_this_locus , ]
    n_sig_snps <- NROW(sig_ss)
  }
  peaks <- as.data.table( do.call( rbind, peaks0 ) )
  
  
  #-------------------------------------------------------------------------------
  #   Format and write outputs
  #-------------------------------------------------------------------------------
  
  # Dump GWAS P values
  # Column names/order: SNP, P, N
  message2("Dump GWAS P values")
  gw_outfile <- file.path( maindir, "gwas_pvalues.tsv" )
  gw_out <- gw3[ , c( "SNP", "P", "N" ) ]
  fwrite( x=gw_out, file=gw_outfile, sep="\t" )
  
  # Dump a file of SNP locations
  # No header, but column order: SNP, CHR, BP
  message2("Dump a file of SNP locations")
  snp_loc_outfile <- file.path( maindir, "snp_locations.tsv" )
  snp_loc_out <- gw3[ , c( "SNP", "chr", "bp" ) ]
  fwrite( x=snp_loc_out, file=snp_loc_outfile, sep="\t", col.names=FALSE )
  
  # Dump peaks
  message2("Dump peaks")
  peaks_outfile <- file.path( maindir, "peaks.tsv" )
  fwrite( x=peaks, file=peaks_outfile, sep="\t" )
  
  
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
  cmd <- paste( "/projects/0/prjs0817/software2/magma/magma",
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

run_magma <- function( maindir, ld.panel, population ){
  
  # Create a job identifier
  job.id <- paste0( "m", sample( x=1:999, size=1 ) )
  
  # If using HRC, run each chromosome separately
  # If using 1000 Genomes, run all at once
  if( ld.panel == "hrc" ){
    
    # HRC: loop through chromosomes 
    for( CHR in 1:22 ){
      
      # Create job name and log file name
      jobname <- paste0( job.id, ".", CHR )
      logfile <- paste0( maindir, "/logs/magma", CHR, ".log" )
      
      # Run
      message2( "Submitting job to the cluster for chromosome: ", CHR )
      cmd <- paste( "sbatch",
                    "-J", jobname,
                    "-o", logfile,
                    "-e", logfile,
                    "/projects/0/prjs0817/repos/brett/e1_run_magma_hrc.sh",
                    CHR, maindir, population )
      system(cmd)
    }
  }else if( ld.panel == "g1000" ){
    
    # 1000 Genomes
    message2("Submitting job to the cluster")
    logfile <- paste0( maindir, "/logs/magma.log" )
    cmd <- paste( "sbatch",
                  "-J", job.id,
                  "-o", logfile,
                  "-e", logfile,
                  "/projects/0/prjs0817/repos/brett/e2_run_magma_g1000.sh",
                  maindir )
    system(cmd)
  }else{
    stop("ld.panel must be either 'hrc' or 'g1000'")
  }
  
  # Wait until all jobs are finished before proceeding
  ensure_finished_jobs(job.id)
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

magma_plots <- function( maindir, z.or.p="z" ){
  
  # Read in peaks
  message2("Read in peaks")
  library(data.table)
  peaks_file <- file.path( maindir, "peaks.tsv" )
  peaks <- fread(peaks_file)
  # peaks <- peaks[ order( peaks$chr, peaks$bp ) , ]
  
  # Read in MAGMA files
  message2("Read in MAGMA files")
  mag_file <- file.path( maindir, "magma.genes.out" )
  mag  <- fread(mag_file)
  
  # Use z-score or P value as the plot's y-axis, depending on choice
  message2("Establish plot y-axis (z-score or -log10 P value)")
  if( z.or.p == "z" ){
    mag$Y <- abs(mag$ZSTAT)
    ylab <- "Z-score"
  }else if( z.or.p == "p" ){
    mag$Y <- -log10(mag$P)
    ylab <- "-log10 P value"
  }else{
    stop("z.or.p must be either 'z' or 'p'")
  }
  
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
  for( i in seq_len( NROW(peaks) ) ){
    
    # Subset MAGMA results
    # message2("Starting locus: ", i, "/", NROW(peaks) )
    buffer_Mbp <- 0.2
    xmin <- peaks$lo[i]/1e6 - buffer_Mbp
    xmax <- peaks$hi[i]/1e6 + buffer_Mbp
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
    xlab <- paste0( "Chr", peaks$chr[i], " position (Mbp)")
    par( mar=c( 4, 4, 0.5, 0.5 ) )
    plot( x=locus$START, y=locus$Y, xlim=c(xmin,xmax), ylim=c(0,ymax),
          xlab=xlab, ylab=ylab, las=1, type="n" )
    
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
  
  # Run
  bash_script <- "/projects/0/prjs0817/repos/brett/f_run_pops.sh"
  cmd <- paste( "sbatch", 
                "-J", jobname,
                "-o", logfile,
                "-e", logfile,
                bash_script, 
                maindir )
  system(cmd)
  ensure_finished_jobs(jobname)
}


#-------------------------------------------------------------------------------
#   pops_plots
#-------------------------------------------------------------------------------

pops_plots <- function( maindir, z.or.p="z" ){
  
  # Read in peaks
  message2("Read in peaks")
  library(data.table)
  peaks_file <- file.path( maindir, "peaks.tsv" )
  peaks <- fread(peaks_file)
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
  if( z.or.p == "z" ){
    pops$Y <- pops$pops
    ylab <- "PoP Score"
  }else if( z.or.p == "p" ){
    pops$Y <- -log10(pops$p)
    ylab <- "-log10 PoPS P value"
  }else{
    stop("z.or.p must be either 'z' or 'p'")
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
  for( i in seq_len( NROW(peaks) ) ){
    
    # Subset PoPS results
    # message2("Starting locus: ", i, "/", NROW(peaks) )
    buffer_Mbp <- 0.2
    xmin <- peaks$lo[i]/1e6 - buffer_Mbp
    xmax <- peaks$hi[i]/1e6 + buffer_Mbp
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
    xlab <- paste0( "Chr", peaks$chr[i], " position (Mbp)")
    par( mar=c( 4, 4, 0.5, 0.5 ) )
    plot( x=locus$START, y=locus$Y, xlim=c(xmin,xmax), ylim=c(ymin,ymax),
          xlab=xlab, ylab=ylab, las=1, type="n" )
    
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
                   ld.panel   = "hrc",
                   population = "eur",
                   gw.file    = "/projects/0/prjs0817/projects/pops/analyses/pd/meta5_raw.tab.gz",
                   chr.bp.col = "SNP",
                   chr.col    = NULL,
                   bp.col     = NULL,
                   a1.col     = "A1",
                   a2.col     = "A2",
                   p.col      = "p",
                   eaf.col    = "freq",
                   n1.col     = "N_cases",
                   n0.col     = "N_controls",
                   n.col      = NULL,
                   n          = NULL,
                   z.or.p     = "z",
                   check.args = TRUE ){
  
  
  #-------------------------------------------------------------------------------
  #   TODOs
  #-------------------------------------------------------------------------------
  
  # format_gwas_and_snp_loc_files
  #    1. Remove duplicated CPRA in the GWAS file
  
  
  #-------------------------------------------------------------------------------
  #   Input descriptions
  #-------------------------------------------------------------------------------
  
  #   maindir:    Main directory in which to store results
  #   ld.panel:   Which LD reference panel should be used? Options are either: 
  #               "hrc" or "g1000".
  #   population: Which populations should be used? Options are either: 
  #               "eur" or "eas".
  #   gw.file:    GWAS file name
  #   chr.bp.col: Optional. The name of a GWAS column containing chromosome and bp 
  #               information separated by a punctuation character. Must be
  #               specified if chr.col and bp.col are not specified.
  #   chr.col:    Optional. The name of a GWAS column containing chromosome
  #               information. Must be specified if chr.bp.col is not specified.
  #   bp.col:     Optional. The name of a GWAS column containing position (bp)
  #               information. Must be specified if chr.bp.col is not specified.
  #   a1.col:     The name of a GWAS column containing the effect (alt) allele
  #   a2.col:     The name of a GWAS column containing the non-effect (ref) allele
  #   p.col:      The name of a GWAS column containing the P value
  #   eaf.col     Optional. The name of a GWAS column containing the frequency of
  #               the effect (alt) allele. If not specified, all palindromic SNPs
  #               will be removed. Otherwise palindromic SNPs with similar 
  #               frequencies in GWAS and HRC will be preserved.
  #   n1.col:     The name of a GWAS column containing the per-SNP number of 
  #               cases. Must be specified with n0.col, or must specify n.col 
  #               or n.
  #   n0.col:     The name of a GWAS column containing the per-SNP number of 
  #               controls. Must be specified with n1.col, or must specify n.col 
  #               or n.
  #   n.col:      The name of a GWAS column containing the per-SNP effective
  #               sample size. If not specified, must specify both n1.col and 
  #               n0.col, or n.
  #   n:          If per-SNP sample size information is not available, this
  #               specifies the study-wide effective sample size. If not 
  #               specified, must specify both n1.col and n0.col, or n.col.
  #   z.or.p:     Should plots use the POPS (and MAGMA) z-scores or P values.
  #               Must be either "z" or "p".
  #   check.args: Logical. Check whether arguments are valid?
  
  
  #-------------------------------------------------------------------------------
  #   Print inputs
  #-------------------------------------------------------------------------------
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/brett/z_brett.R")
  message2("Load libraries and sources")
  suppressPackageStartupMessages( library(data.table) )
  
  # Print inputs
  message_header("Print inputs")
  message2( "Main directory: ", maindir )
  message2( "LD reference panel: ", ld.panel )
  message2( "Population: ", population )
  message2( "GWAS file: ", gw.file )
  message2( "Chromosome/position column name: ", chr.bp.col)
  message2( "Chromosome column name: ", chr.col )
  message2( "Position column name: ", bp.col )
  message2( "Effect allele column name: ", a1.col )
  message2( "Non-effect allele column name: ", a2.col )
  message2( "P value column name: ", p.col )
  message2( "Effect allele frequency column name: ", eaf.col )
  message2( "Number of cases column name: ", n1.col )
  message2( "Number of controls column name: ", n0.col )
  message2( "Effective sample size column name: ", n.col )
  message2( "Effective sample size: ", n )
  message2( "Use z-score or P value?: ", z.or.p )
  
  
  #-------------------------------------------------------------------------------
  #   Check arguments
  #-------------------------------------------------------------------------------
  
  if(check.args){
    message_header("Check arguments")
    message2("Checking arguments")
    check_arguments( ld.panel   = ld.panel,
                     population = population,
                     gw.file    = gw.file,
                     chr.bp.col = chr.bp.col,
                     chr.col    = chr.col,
                     bp.col     = bp.col,
                     a1.col     = a1.col,
                     a2.col     = a2.col,
                     p.col      = p.col,
                     eaf.col    = eaf.col,
                     n1.col     = n1.col,
                     n0.col     = n0.col,
                     n.col      = n.col,
                     n          = n,
                     z.or.p     = z.or.p )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Assign output file names
  #-------------------------------------------------------------------------------
  
  # Make output and log directories
  logdir <- file.path( maindir, "logs" )
  dir.create( path=logdir, showWarnings=FALSE, recursive=TRUE )
  
  # Assign output file names
  clean_gw.file         <- file.path( maindir, "gwas_pvalues.tsv" )
  snp_loc_file          <- file.path( maindir, "snp_locations.tsv" )
  peaks_file            <- file.path( maindir, "peaks.tsv" )
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
  if( ld.panel == "hrc" ){
    mag_sumstats_files    <- file.path( maindir, "magma", paste0( "chr", 1:22, ".genes.out" ) )
    mag_covar_files       <- file.path( maindir, "magma", paste0( "chr", 1:22, ".genes.raw" ) )
  }else if( ld.panel == "g1000" ){
    mag_sumstats_files    <- mag_ss_collated_file
    mag_covar_files       <- mag_cov_collated_file
  }else{
    stop("ld.panel must be either 'hrc' or 'g1000'")
  }
  
  
  #-------------------------------------------------------------------------------
  #   Format the GWAS
  #   Write the SNP location file
  #   Write the peaks file
  #-------------------------------------------------------------------------------
  
  message_header("Format the GWAS, write the SNP location and peaks files")
  if( all( file.exists( clean_gw.file, 
                        snp_loc_file,
                        peaks_file ) ) ){
    message2("Output files exist, skipping")
    
  }else{
    format_gwas_and_snp_loc_files( maindir    = maindir,
                                   ld.panel   = ld.panel,
                                   population = population,
                                   gw.file    = gw.file,
                                   chr.bp.col = chr.bp.col,
                                   chr.col    = chr.col,
                                   bp.col     = bp.col,
                                   a1.col     = a1.col,
                                   a2.col     = a2.col,
                                   p.col      = p.col,
                                   eaf.col    = eaf.col,
                                   n1.col     = n1.col,
                                   n0.col     = n0.col,
                                   n.col      = n.col,
                                   n          = n )
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
               ld.panel   = ld.panel,
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
    magma_plots( maindir = maindir, 
                 z.or.p  = z.or.p )
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
    pops_plots( maindir = maindir, 
                z.or.p  = z.or.p )
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
                  ld.panel   = ld.panel,
                  gw.file    = gw.file,
                  chr.bp.col = chr.bp.col,
                  chr.col    = chr.col,
                  bp.col     = bp.col,
                  a1.col     = a1.col,
                  a2.col     = a2.col,
                  p.col      = p.col,
                  eaf.col    = eaf.col,
                  n1.col     = n1.col,
                  n0.col     = n0.col,
                  n.col      = n.col,
                  n          = n,
                  z.or.p     = z.or.p,
                  check.args = check.args )
    render( input       = rmd_file, 
            params      = args, 
            output_file = html_file )
  }
  message_header("Done")
}


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------



