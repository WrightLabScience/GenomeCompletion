#!/usr/bin/env Rscript

# Usage: 
# Rscript unmap_annot.1.r --genomeQ='GCF_045343885' --genomeR='GCF_030291915' --assembly_summary=assembly_summary.rds --annot
# map pattern to subject

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--genomeQ", type="character", help="genomeQ, query (Required)"),
	make_option("--genomeR", type="character", help="genomeR, ref (Required)"),
	make_option("--ir_gaps_query", type="character", help="ir_gaps_query from nucmer, genomeQ_genomeR.rds (Required)"),
	make_option("--min_gap_size", type="integer", default=150, help="minimum gap size"),
	make_option("--gff_dir", type="character", default=".", help="source data dir for gff.gz"),
	make_option("--gff_ext", type="character", default=".gff.gz", help="gff filename extension"),
	make_option("--assembly_summary", type="character", default=NA, help="assembly_summary.rds file with ftp_path for genomes")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")
suppressMessages(library(SynExtend))
packageVersion("SynExtend")
suppressMessages(library(ape))
packageVersion("ape")
options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## function ########################
# download source data
fetch_data <- function(seqID, out_ext, FTP_dir, FTP_ext, attempt_n=5) {
	out_file <- paste0(seqID, out_ext)
	FTPadd <- paste0(FTP_dir, gsub('.*\\/(.+)\\/$', '\\1', FTP_dir), FTP_ext)
	cat("\nDownloading source data from: ", FTPadd,'\n')
	# try downloading the file until succeed
	attempt <- 1
	while (!file.exists(out_file) && attempt <= attempt_n) {
		attempt <- attempt + 1
		Sys.sleep(30)
		try(download.file(FTPadd, out_file))
	}
	if (!file.exists(out_file)) {
		cat('\n!!!!!!!!!!!!!!!', seqID, 'failed to download source data from FTP.\n')
	} else {
		cat('# Done downloading source data file: ', out_file, '\n')
	}
}

######################## main ########################
# parameters
min_gap_size <- opt$min_gap_size
genomeQ <- opt$genomeQ
genomeR <- opt$genomeR
assembly_summary <- readRDS(opt$assembly_summary)

# gaps between merged query blocks
ir_gaps_query <- readRDS(opt$ir_gaps_query)
cat('\n# Number of gaps before size filtering:', length(ir_gaps_query), '\n')
ir_unmapped <- ir_gaps_query[width(ir_gaps_query)>=min_gap_size]
cat('\n# Number of gaps after size filtering:', length(ir_unmapped), '\n')

### overlay with gff feature ###
# gff for subject
if (file.exists(paste0(gsub("/$", "", opt$gff_dir),'/', genomeR, opt$gff_ext))) {
	cat('genomeR gff exists:', genomeR, '\n')
} else {
	FTP_dir <- assembly_summary[genomeR, 'ftp_path']
	fetch_data(genomeR, out_ext=opt$gff_ext, FTP_dir, FTP_ext='_genomic.gff.gz', attempt_n=5)
}
# read gff
gff_Q <- read.gff(paste0(gsub("/$", "", opt$gff_dir),'/', genomeR, opt$gff_ext), na.strings = c(".", "?"), GFF3 = TRUE)
genome_regions <- gff_Q[gff_Q$type=='region',c('seqid','start','end')]

gff_Q <- gff_Q[gff_Q$type!='region',] # remove first row for whole genome
if (length(unique(gff_Q$seqid))>1) {
	cat('STOP!!!!!\n')
	quit(save = "no", status = 0, runLast = FALSE)
}
ir_gff <- IRanges(start = gff_Q$start, end = gff_Q$end)

for (type_h in unique(gff_Q$type)) {
	cat('# Type: ', type_h, '\n', sep='')
}

# find genes overlap with unmapped regions
overlap_hits <- findOverlaps(ir_gff, ir_unmapped) # findOverlaps(query, subject)
# subset GFF to only contains unmapped genes
unmapped_gff <- gff_Q[unique(queryHits(overlap_hits)), ]
if (nrow(unmapped_gff)>0) {
	saveRDS(unmapped_gff, file = paste0(genomeQ,'_',genomeR,"_unmapped_gff.rds"))
	# extract list of non-pseudogene CDS to get COG
	unmap_CDS_gff <- unmapped_gff[(unmapped_gff$type=='CDS' & !grepl('pseudo=true',unmapped_gff$attributes)),] # non-pseudo CDS
	if (sum(gff_Q$type=='CDS')<1) {
		cat('!!!!!!!!!!!!!!! NO CDS in GFF.\n')
	} else if (nrow(unmap_CDS_gff)!=0) {
		unmap_CDS_gff$protein_id <- gsub(".*protein_id=([^;]+).*", "\\1", unmap_CDS_gff$attributes)
		unmap_CDS_gff$product <- gsub('.*;product=([^;]+);.*', '\\1', unmap_CDS_gff$attributes)
		unmap_CDS_gff$genomeQ <- genomeQ
		unmap_CDS_gff$genomeR <- genomeR
		write.table(x = unmap_CDS_gff[,c('protein_id','product','genomeQ','genomeR')],file = paste0(genomeQ,'_',genomeR,"_protein_id.tsv"), quote = FALSE, row.names = FALSE, append = FALSE, sep = "\t", col.names = FALSE)
	} else {
		cat('!!!No CDS in unmapped regions.\n')
	}
} else {
	cat('!!!No gene in unmapped regions.\n')
}


###### calculate values for plotting Fig2C ######
gff_CDS <- gff_Q[(gff_Q$type=='CDS' & !grepl('pseudo=true',gff_Q$attributes)),] # non-pseudo CDS
ir_wholegenome_CDS <- IRanges(start = gff_CDS$start, end = gff_CDS$end)
ir_wholegenome_CDS <- reduce(ir_wholegenome_CDS)
# a = total number of coding bases
a = sum(width(ir_wholegenome_CDS))
# b = total number of coding bases unmapped
b = sum(width(intersect(ir_unmapped, ir_wholegenome_CDS)))
# X = b/a = fraction coding ummapped
X = b/a

# c = total number of non-coding bases
c = genome_regions[1,'end'] - a # whole genome size - total number of coding bases
# d = total number of non-coding bases unmapped
d = sum(width(ir_unmapped)) - b # unmapped bases - unmapped CDS bases
# Y =  d/c = fraction non-coding unmmaped
Y = d/c


# report summary
output_title <- paste(
	'### Title',					# 0
	'genomeQ',						# 1
	'genomeR',						# 2
	'total_CDS_bases_a',
	'unmapped_CDS_bases_b',
	'fract_CDS_unmapped_X',
	'total_NC_bases_c',
	'unmapped_NC_bases_d',
	'frac_NC_unmapped_Y',
	'genomeQ_size',
	'unmapped_bases',				# 5
	'total_gene_count',				# 10
	'total_pseudogene_count',		# 11
	'total_CDS_count',				# 12
	'total_HP_count',				# 13
	'total_rRNA_count',				# 14
	'total_tRNA_count',				# 15
	'total_tmRNA_count',			# 16
	'total_SRP_RNA_count',			# 17
	'total_RNase_P_RNA_count',		# 18
	'total_riboswitch_count',		# 19
	'unmapped_gene_count',			# 21
	'unmapped_pseudogene_count',	# 22
	'unmapped_CDS_count',			# 23
	'unmapped_HP_count',			# 24
	'unmapped_rRNA_count',			# 25
	'unmapped_tRNA_count',			# 26
	'unmapped_tmRNA_count',			# 27
	'unmapped_SRP_RNA_count',		# 28
	'unmapped_RNase_P_RNA_count',	# 29
	'unmapped_riboswitch_count',	# 30
	sep='\t')
output_value <- paste(
	'### Result_unmap',									# 0
	genomeQ,								# 1
	genomeR,								# 2
	a,
	b,
	X,
	c,
	d,
	Y,
	genome_regions[1,'end'],
	sum(width(ir_unmapped)),					# 5
	sum(gff_Q$type=='gene'),					# 10
	sum(gff_Q$type=='pseudogene'),			# 11
	sum(gff_Q$type=='CDS'),					# 12
	sum(grepl("hypothetical protein", gff_Q$attributes, ignore.case = TRUE), na.rm = TRUE), # 13
	sum(gff_Q$type=='rRNA'),					# 14
	sum(gff_Q$type=='tRNA'),					# 15
	sum(gff_Q$type=='tmRNA'),					# 16
	sum(gff_Q$type=='SRP_RNA'),				# 17
	sum(gff_Q$type=='RNase_P_RNA'),			# 18
	sum(gff_Q$type=='riboswitch'),			# 19
	sum(unmapped_gff$type=='gene'),					# 21
	sum(unmapped_gff$type=='pseudogene'),			# 22
	sum(unmapped_gff$type=='CDS'),					# 23
	sum(grepl("hypothetical protein", unmapped_gff$attributes, ignore.case = TRUE), na.rm = TRUE), # 24
	sum(unmapped_gff$type=='rRNA'),					# 25
	sum(unmapped_gff$type=='tRNA'),					# 26
	sum(unmapped_gff$type=='tmRNA'),				# 27
	sum(unmapped_gff$type=='SRP_RNA'),				# 28
	sum(unmapped_gff$type=='RNase_P_RNA'),			# 29
	sum(unmapped_gff$type=='riboswitch'),			# 30
	sep='\t')
cat('\n', output_title, '\n', output_value, '\n\n\n', sep='')