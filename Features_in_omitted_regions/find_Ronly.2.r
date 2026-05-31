#!/usr/bin/env Rscript

# Usage: 
# Rscript find_Ronly.1.r --genomeQ='GCF_045343885' --genomeR='GCF_030291915'
# map pattern to subject

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--genomeQ", type="character", help="genomeQ, query (Required)", metavar="genomeID"),
	make_option("--genomeR", type="character", help="genomeR, ref (Required)", metavar="genomeID"),
	make_option("--gff_dir", type="character", default=".", help="source data dir for gff.gz", metavar="dirpath"),
	make_option("--gff_ext", type="character", default=".gff.gz", help="gff filename extension", metavar="gff_ext")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")
suppressMessages(library(ape))
packageVersion("ape")
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## function ########################



######################## main ########################
# get names
genomeQ <- opt$genomeQ
genomeR <- opt$genomeR
gff_dir <- opt$gff_dir

## check and prepare source data
# gff
if (!file.exists(paste0(gsub("/$", "", gff_dir),'/', genomeR, opt$gff_ext))) {
	cat('genomeID gff NOT exists:', genomeR, '\n')
	quit(save = "no", status = 1)    
} else {
	# read input
	gff <- read.gff(paste0(gsub("/$", "", gff_dir),'/', genomeR, opt$gff_ext), na.strings = c(".", "?"), GFF3 = TRUE)
	gff <- gff[gff$type!='region',] # remove first row for whole genome
}

aligned_q_ori <- readRDS(paste0(genomeQ,'_',genomeR,'_aligned_q.rds'))
aligned_r_ori <- readRDS(paste0(genomeQ,'_',genomeR,'_aligned_r.rds'))

# discard any reads with any matches less than the cutoff (0.99) and greater than some minimum threshold (0.9). These hits create ambiguity.
reads2rm <- unique(c(aligned_q_ori[aligned_q_ori$PID<0.99,'Pattern'], aligned_r_ori[aligned_r_ori$PID<0.99,'Pattern']))
aligned_q_ori <- aligned_q_ori[!aligned_q_ori$Pattern %in% reads2rm,]
aligned_r_ori <- aligned_r_ori[!aligned_r_ori$Pattern %in% reads2rm,]
readIDs_map2q <- aligned_q_ori$Pattern
readIDs_map2r <- aligned_r_ori$Pattern
df_Ronly <- aligned_r_ori[aligned_r_ori$Pattern %in% setdiff(unique(readIDs_map2r), unique(readIDs_map2q)), ]
rm(aligned_q_ori, aligned_r_ori)
gc()

# normalize read coverage: weighted coverage
# weight = 1/num of hits for that read
hits_per_read <- table(df_Ronly$Pattern) # num of hits per read
df_Ronly$weight <- 1 / hits_per_read[as.character(df_Ronly$Pattern)] # assign weight
# add coverage
ir_Ronly <- IRanges(start = df_Ronly$SubjectStart, end = df_Ronly$SubjectEnd)
cov <- coverage(ir_Ronly, weight = df_Ronly$weight)
cov_vec <- as.numeric(cov)
print(summary(cov_vec)) # read mapped to same regions?
ir_Ronly_reduced <- reduce(ir_Ronly)
mode_mapped_read_len <- as.numeric(names(sort(-table(round(df_Ronly$AlignmentLength))))[1])

###### compute values for plotting frac CDS vs. NC ######
gff_CDS <- gff[(gff$type=='CDS' & !grepl('pseudo=true',gff$attributes)),] # non-pseudo CDS
ir_wholegenome_CDS <- IRanges(start = gff_CDS$start, end = gff_CDS$end)
ir_wholegenome_CDS_reduced <- reduce(ir_wholegenome_CDS)

# annotate Ronly all gff
ir_gff <- IRanges(start = gff$start, end = gff$end)
# find overlaps
overlap_hits <- findOverlaps(ir_gff, ir_Ronly_reduced) # findOverlaps(query, subject)
# subset GFF to only contains Ronly regions
gff_Ronly <- gff[unique(queryHits(overlap_hits)), ]


###### compute values for plotting COG enrichment ######
# annotate Ronly CDS gff
# find overlaps
overlap_hits <- findOverlaps(ir_wholegenome_CDS, ir_Ronly_reduced) # findOverlaps(query, subject)
# subset GFF to only contains Ronly regions
gff_Ronly_CDS <- gff_CDS[unique(queryHits(overlap_hits)), ]

# extract list of non-pseudogene CDS to get COG
if (nrow(gff_Ronly_CDS)>0) {
	gff_Ronly_CDS$protein_id <- gsub(".*protein_id=([^;]+).*", "\\1", gff_Ronly_CDS$attributes)
	gff_Ronly_CDS$product <- gsub('.*product=([^;]+).*', '\\1', gff_Ronly_CDS$attributes)
	gff_Ronly_CDS$genomeQ <- genomeQ
	gff_Ronly_CDS$genomeR <- genomeR
	save(df_Ronly, cov_vec, gff_Ronly, gff_Ronly_CDS, file=paste0(genomeQ,'_',genomeR,'_Ronly.RData'), compress='xz')
} else {
	cat('!!!No CDS in Ronly regions.\n')
	save(df_Ronly, cov_vec, gff_Ronly, gff_Ronly_CDS, file=paste0(genomeQ,'_',genomeR,'_Ronly.RData'), compress='xz')
}