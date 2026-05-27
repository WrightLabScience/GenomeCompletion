#!/usr/bin/env Rscript

# Usage: 
# Rscript count_CDS_bases.1.r --genomeID='GCF_045343885' --faa_dir=./ --assembly_summary=assembly_summary_1ctg_s50.rds
# get protein fasta from protein id table

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--genomeID", type="character", help="genomeID (Required)"),
	make_option("--gff_dir", type="character", default=".", help="source data dir for gff.gz", metavar="dirpath"),
    make_option("--gff_ext", type="character", default=".gff.gz", help="gff filename extension", metavar="gff_ext"),
	make_option("--assembly_summary", type="character", default="assembly_summary.rds", help="assembly_summary.rds file with ftp_path for genomes", metavar="assembly_summary")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(IRanges))
suppressMessages(library(ape))
options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## function ########################
# download source data
fetch_data <- function(genomeID, out_ext, FTP_dir, FTP_ext, attempt_n=5) {
	out_file <- paste0(genomeID, out_ext)
	FTPadd <- paste0(FTP_dir,'/', gsub('.*\\/(.*)$', '\\1', FTP_dir), FTP_ext)
	cat("\nDownloading source data from: ", FTPadd,'\n')
	# try downloading the file until succeed
	attempt <- 1
	while (!file.exists(out_file) && attempt <= attempt_n) {
		attempt <- attempt + 1
		Sys.sleep(1)
		try(download.file(FTPadd, out_file))
	}
	if (!file.exists(out_file)) {
		cat('\n!!!!!!!!!!!!!!!', genomeID, 'failed to download source data from FTP.\n')
	} else {
		cat('# Done downloading source data file: ', out_file, '\n')
	}
}

######################## main ########################
# read input
genomeID <- opt$genomeID


## check and prepare source data: gff
if (file.exists(paste0(gsub("/$", "", opt$gff_dir),'/', genomeID, opt$gff_ext))) {
	cat('# genomeID gff exists:', genomeID, '\n')
} else {
	assembly_summary <- readRDS(opt$assembly_summary)
	FTP_dir <- assembly_summary[genomeID, 'ftp_path']
	fetch_data(genomeID, out_ext=opt$gff_ext, FTP_dir, FTP_ext='_genomic.gff.gz', attempt_n=5)
}

# list of non-pseudo CDS proteinIDs
gff_path <- paste0(genomeID, opt$gff_ext)
gff_h <- read.gff(gff_path, na.strings = c(".", "?"), GFF3 = TRUE)
gff_CDS <- gff_h[(gff_h$type=='CDS' & !grepl('pseudo=true',gff_h$attributes)),] # non-pseudo CDS
gff_CDS$protein_id <- gsub(".*protein_id=([^;]+).*", "\\1", gff_CDS$attributes)
cat('# genomeID:', genomeID, '; total_CDS_count:', length(gff_CDS$protein_id), '; unique_CDS_count:', length(unique(gff_CDS$protein_id)), '\n')
writeLines(gff_CDS$protein_id, paste0(genomeID, '_proteinIDs.txt'))

# calculate genic(CDS) vs intergenic regions
genome_regions <- gff_h[gff_h$type=='region',c('seqid','start','end')]
total_genome_bases <- sum(genome_regions$end)
total_CDS_bases <- 0
total_NC_bases <- 0
for (seqid_h in genome_regions$seqid) {
    gff_CDS_h <- gff_CDS[gff_CDS$seqid==seqid_h,]
    ir <- IRanges(start = gff_CDS_h$start, end = gff_CDS_h$end)
    merged_ir <- reduce(ir)
    CDS_bases <- sum(width(merged_ir))
    total_CDS_bases <- total_CDS_bases+CDS_bases
    contig_size <- as.numeric(genome_regions[genome_regions$seqid==seqid_h,'end'])
    # uncovered_ir <- gaps(merged_ir, start=1, end=contig_size)
    # print(sum(width(uncovered_ir)))
    NC_bases <- contig_size - CDS_bases
    total_NC_bases <- total_NC_bases+NC_bases
}
cat('# genomeID (genome_size):', genomeID, total_genome_bases,'; total_CDS_bases (count/frac):', total_CDS_bases, total_CDS_bases/total_genome_bases, '; total_NC_bases (count/frac):', total_NC_bases, total_NC_bases/total_genome_bases, '\n')