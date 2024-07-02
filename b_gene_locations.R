
#   b_gene_locations.R


#--------------------------------------------------------------------------------
#   Load POPS genes
#--------------------------------------------------------------------------------

# Read in POPS genes
library(data.table)
pops_file <- "/projects/0/prjs0817/projects/pops/data/features_split/features.1"
pops_cmd <- paste( "awk '{print $1}'", pops_file )
pops <- fread( cmd=pops_cmd )


#--------------------------------------------------------------------------------
#   Load gencode genes, create columns
#--------------------------------------------------------------------------------

# Read in gencode v44 file
gc_file <- "/projects/0/prjs0817/projects/pops/data/gencode.v44.grch37.gff3.gz"
gc_cmd  <- paste( 'zcat', gc_file, '| grep -v "#"' )
gc      <- fread( cmd=gc_cmd )
names(gc) <- c( "CHR", "src", "type", "START", "END", "score", "STRAND", "phase", "attr" )

# Subset to genes
gc2 <- gc[ gc$type == "gene" , ]

# Extract gene ID, type, and HGNC name
gc2$ENSGID <- sub( pattern     = "^ID=(.*)\\.[[:digit:]]*;gene_id=.*",
                   replacement = "\\1",
                   x           = gc2$attr)
gc2$TYPE <- sub( pattern     = "^.*gene_type=([[:alpha:]_]+);gene.*",
                 replacement = "\\1",
                 x           = gc2$attr)
gc2$NAME <- sub( pattern     = ".*gene_name=(.*);level=.*",
                 replacement = "\\1",
                 x           = gc2$attr)

# Make CHR an integer
gc2$CHR <- sub( pattern     = "^chr([[:alnum:]]+)$",
                replacement = "\\1",
                x           = gc2$CHR )

# Make a TSS column (START if on + strand, END if on - strand)
gc2$TSS <- ifelse( gc2$STRAND == "+",
                   gc2$START,
                   gc2$END )


#--------------------------------------------------------------------------------
#   Subset to rows and columns, write
#--------------------------------------------------------------------------------

# Subset
gtypes <- c( "protein_coding", "lncRNA", "miRNA", "snRNA", "snoRNA",
             "transcribed_unprocessed_pseudogene", 
             "transcribed_processed_pseudogene",
             "transcribed_unitary_pseudogene" )
gcols <- c( "ENSGID", "CHR", "START", "END", "TSS", "STRAND", "NAME", "TYPE" )
gc3 <- gc2[ gc2$TYPE %in% gtypes , ..gcols ]
gc4 <- gc3[ match( pops$ENSGID, gc2$ENSGID ) , ..gcols ]
gc4 <- gc4[ !is.na(gc4$ENSGID) , ]

# Write
outfile1 <- "/projects/0/prjs0817/projects/pops/data/gene_locations.tsv"
outfile2 <- "/projects/0/prjs0817/projects/pops/data/gene_locations_inc_nc.tsv"
fwrite( x=gc4, file=outfile1, sep="\t" )
fwrite( x=gc3, file=outfile2, sep="\t" )


#--------------------------------------------------------------------------------
#   Done
#--------------------------------------------------------------------------------



