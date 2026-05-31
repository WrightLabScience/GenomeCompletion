#!/usr/bin/env Rscript

# Usage: 
# Rscript count_feature_bases.3.r --prefix='GCF_045343885'
# Count CDS, ncRNA, pseudoCDS, non-annotated bases
# split pseudo CDS by types

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--prefix", type="character", help="prefix for input and output (Required)"),
	make_option("--Ronly", type="character", default=NA, help="Rdata with Ronly results"),
	make_option("--gff_ext", type="character", default=".gff.gz", help="gff filename extension")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(IRanges))
suppressMessages(library(ape))
options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## function ########################
# for each position in ir_target, count how many types cover it
# use coverage() to get per-base type membership efficiently
get_coverage_vec <- function(ir, ir_target) {
	target_len <- max(end(ir_target))
	if (length(ir) == 0) return(rep(0, sum(width(ir_target))))
	as.integer(coverage(ir, width = target_len)[ir_target] > 0)
}

######################## main ########################
# read input
prefix <- opt$prefix
if (!is.na(opt$Ronly)) {
	load(opt$Ronly) # df_Ronly, cov_vec, gff_Ronly, gff_Ronly_CDS
	if (nrow(df_Ronly)>0) {
		gff_all <- gff_Ronly
		ir_target <- reduce(IRanges(start = df_Ronly$SubjectStart, end = df_Ronly$SubjectEnd))
	} else {
		cat('Ronly is empty!!!\n')
		quit(save = "no", status = 0)
	}
} else { # genomeR
	gff_all <- read.gff(paste0(prefix, opt$gff_ext), na.strings = c(".", "?"), GFF3 = TRUE)
	ir_target <- IRanges(start=1, end=gff_all[gff_all$type=='region','end'])
	writeLines(as.character(unique(gff_all$type)), paste0(prefix, "_gff_types.txt")) # report possible types
}

# output
output_title <- paste(
	'### Title',			# 0
	'prefix',				# 1
	'total_bases',			# 2
	'total_CDS_bases',		# 3
	'total_ncRNA_bases',	# 4
	'total_other_bases',	# 5
	'pseudoFS_bases',		# 6
	'pseudoIS_bases',		# 7
	'pseudoTrunc_bases',	# 8
	sep='\t')

# gff for ranges
if (nrow(gff_all)>0) {
	gff_CDS <- gff_all[(gff_all$type=='CDS' & !grepl('pseudo=true',gff_all$attributes)),] # non-pseudo CDS
	gff_ncRNA <- gff_all[(grepl('RNA',gff_all$type) | grepl('ribo',gff_all$type)),] # ncRNA
	gff_pseudo <- gff_all[(gff_all$type=='CDS' & grepl('pseudo=true',gff_all$attributes, ignore.case=TRUE)),] # pseudo CDS
	# split pseudo by types
	gff_pseudo$Note <- URLdecode(gsub('.*Note=([^;]+).*','\\1',gff_pseudo$attributes))
	gff_pseudo$frameshifted <- grepl("frameshift", gff_pseudo$Note, ignore.case = TRUE)
	gff_pseudo$internal_stop <- grepl("internal stop", gff_pseudo$Note, ignore.case = TRUE)
	gff_pseudo$incomplete <- grepl("incomplete", gff_pseudo$Note, ignore.case = TRUE) # missing N-terminus and/or C-terminus
	# convert to IRanges
	ir_CDS <- if (nrow(gff_CDS) > 0) reduce(IRanges(start = gff_CDS$start, end = gff_CDS$end)) else IRanges()
	ir_CDS <- IRanges::intersect(ir_target, ir_CDS)
	ir_ncRNA <- if (nrow(gff_ncRNA) > 0) reduce(IRanges(start = gff_ncRNA$start, end = gff_ncRNA$end)) else IRanges()
	ir_ncRNA <- IRanges::intersect(ir_target, ir_ncRNA)
	ir_pseudo <- if (nrow(gff_pseudo) > 0) reduce(IRanges(start = gff_pseudo$start, end = gff_pseudo$end)) else IRanges()
	ir_pseudo <- IRanges::intersect(ir_target, ir_pseudo)
	ir_other <- IRanges::setdiff(ir_target, reduce(c(ir_CDS, ir_pseudo, ir_ncRNA))) # unannotated
	ir_pseudoFS <- if (any(gff_pseudo$frameshifted)) reduce(IRanges(start = gff_pseudo$start[gff_pseudo$frameshifted], end = gff_pseudo$end[gff_pseudo$frameshifted])) else IRanges()
	ir_pseudoFS <- IRanges::intersect(ir_target, ir_pseudoFS)
	ir_pseudoIS <- if (any(gff_pseudo$internal_stop)) reduce(IRanges(start = gff_pseudo$start[gff_pseudo$internal_stop], end = gff_pseudo$end[gff_pseudo$internal_stop])) else IRanges()
	ir_pseudoIS <- IRanges::intersect(ir_target, ir_pseudoIS)
	ir_pseudoTrunc <- if (any(gff_pseudo$incomplete)) reduce(IRanges(start = gff_pseudo$start[gff_pseudo$incomplete], end = gff_pseudo$end[gff_pseudo$incomplete])) else IRanges()
	ir_pseudoTrunc <- IRanges::intersect(ir_target, ir_pseudoTrunc)
	# for each feature type, get intersected IRanges with ir_target
	ir_types <- list(
		CDS             = ir_CDS,
		ncRNA           = ir_ncRNA,
		pseudoFS        = ir_pseudoFS,
		pseudoIS        = ir_pseudoIS,
		pseudoTrunc     = ir_pseudoTrunc,
		other           = ir_other
	)
	# build matrix (columns = types, rows = positions)
	mat <- do.call(cbind, lapply(ir_types, get_coverage_vec, ir_target=ir_target))
	# normalize by rowSums (positions covered by multiple types get fractional credit)
	rs <- rowSums(mat)
	rs[rs == 0] <- 1  # avoid division by zero for uncovered positions (shouldn't happen)
	mat_norm <- mat / rs
	# column sums = normalized base counts per type
	result <- colSums(mat_norm)
	
	output_value <- paste(
		ifelse(is.na(opt$Ronly), '### Result_genomeR', '### Result_Ronly'),	# 0
		prefix,					# 1
		nrow(mat_norm),			# 2
		result['CDS'],			# 3
		result['ncRNA'],		# 4
		result['other'],		# 5
		result['pseudoFS'],		# 6
		result['pseudoIS'],		# 7
		result['pseudoTrunc'],	# 8
		sep='\t')
	cat('\n', output_title, '\n', output_value, '\n\n\n', sep='')
} else {
	cat('GFF is empty!!!\n')
	output_value <- paste(
		ifelse(is.na(opt$Ronly), '### Result_genomeR', '### Result_Ronly'),	# 0
		prefix,					# 1
		sum(width(ir_target)),	# 2
		0,						# 3
		0,						# 4
		sum(width(ir_target)),	# 5
		0,						# 6
		0,						# 7
		0,						# 8
		sep='\t')
	cat('\n', output_title, '\n', output_value, '\n\n\n', sep='')
}