#!/usr/bin/env Rscript

# Usage: 
# Rscript get_COGs.3.r --protein_id_tsv='GCF_045343885_GCF_030291915_protein_id.tsv' --assembly_summary=assembly_summary_1ctg_s50.rds
# get protein fasta from protein id table

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--protein_id_tsv", type="character", help="protein_id_tsv (Required)", metavar="genomeID"),
	make_option("--faa_dir", type="character", default=".", help="source data dir for faa.gz", metavar="dirpath"),
    make_option("--faa_ext", type="character", default=".faa.gz", help="faa filename extension", metavar="faa_ext"),
	make_option("--db", type="character", default="~/scrdir/cdd_11242025/Cog", help="db path for rpsblast", metavar="db"),
	make_option("--cog_tab", type="character", default="cog24col.tsv", help="COG table", metavar="COG"),
	make_option("--assembly_summary", type="character", help="assembly_summary.rds file with ftp_path for genomes", metavar="assembly_summary")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
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
# read input
unmapped_CDS_df <- read.table(opt$protein_id_tsv, sep="\t", quote="", header=FALSE, fill=FALSE, col.names=c('protein_id','product','genomeQ','genomeR'))
cat('# unmapped_CDS_df dim:', dim(unmapped_CDS_df), '\n')

genomeID <- unique(unmapped_CDS_df$genomeQ)
assembly_summary <- readRDS(opt$assembly_summary)
prefix <- sub('_protein_id.tsv','',opt$protein_id_tsv)


## check and prepare source data
if (file.exists(paste0(gsub("/$", "", opt$faa_dir),'/', genomeID, opt$faa_ext))) {
	cat('genomeID faa exists:', genomeID, '\n')
} else {
	FTP_dir <- assembly_summary[genomeID, 'ftp_path']
	fetch_data(genomeID, out_ext=opt$faa_ext, FTP_dir, FTP_ext='_protein.faa.gz', attempt_n=5)
}

# export unmapped faa seqs
faa <- readAAStringSet(paste0(gsub("/$", "", opt$faa_dir),'/', genomeID, opt$faa_ext))
faa_names <- sub(" .*", "", names(faa)) # protein_id
names(faa) <- sub("\\s\\[[^]]+\\]$", "", sub(".*MULTISPECIES:\\s", "", names(faa))) # product
unmapped_faa <- faa[faa_names %in% unmapped_CDS_df$protein_id]
writeXStringSet(unmapped_faa, paste0(prefix,'.faa.gz'), compress=TRUE)

names(unmapped_faa) <- faa_names[faa_names %in% unmapped_CDS_df$protein_id]
writeXStringSet(unmapped_faa, 'unmapped.faa')

# run rpsblast and process result in bash
RPSBLAST <- paste0('rpsblast',
                   ' -query unmapped.faa',
				   ' -db ', opt$db,
				   ' -out results_rps.tsv',
				   ' -evalue 1e-5',
				   ' -outfmt "6 qseqid evalue bitscore stitle"')
cat('### Running BLAST ...\n')
start_time <- Sys.time()
system(command = RPSBLAST, timeout = 999999999)
end_time <- Sys.time()
elapsed_sec <- as.numeric(end_time - start_time, units = "secs")
cat(sprintf("### BLAST time: %.2f min (%.2f sec)\n", elapsed_sec/60, elapsed_sec))
system(command = "sort -k1,1 -k2,2g -k3,3gr results_rps.tsv | awk '!seen[$1]++' > results_rps_f.tsv")

# read blast hits into R
cog24.fun <- read.table(opt$cog_tab, sep="\t", quote="", comment.char='', header=TRUE)

lines <- readLines('results_rps_f.tsv', warn = FALSE)
if (length(lines) == 0 || all(grepl("^\\s*$", lines))) {
    message("## BLAST file empty or contains no data.")
	cat('\n### Result_COGs\t', prefix, '\t', nrow(unmapped_CDS_df), '\t', paste(rep(0,26), collapse='\t'),'\t',nrow(unmapped_CDS_df),'\n\n\n', sep='')
} else {
    rps <- read.table('results_rps_f.tsv', sep="\t", quote="", header=FALSE, col.names=c('protein_id','evalue','bitscore','cog_info'))
	rps2 <- rps %>%
		extract(cog_info,
				into = c("COG","gene","description","categories"),
				regex = "^(COG[0-9]+), ([^,]+), ([^\\[]+) \\[([^]]+)\\].*",
				remove=TRUE)
	for (categoryID_h in rownames(cog24.fun)) {
		description_h <- cog24.fun[categoryID_h, 'description']
		rps2[,categoryID_h] <- grepl(description_h, rps2$categories, ignore.case = TRUE)
	}
	unmapped_rps_df <- merge(unmapped_CDS_df, rps2, by='protein_id', all.x=TRUE)
	unmapped_rps_df[,c(11:36)][is.na(unmapped_rps_df[,c(11:36)])] <- 0
	df_h <- unmapped_rps_df[,c(11:36)]
	df_h$noCOG <- ifelse(rowSums(df_h)<1, 1, 0) # for no blast hits
	df_h <- df_h/rowSums(df_h) # normalize
	unmapped_rps_df <- cbind(unmapped_rps_df[,c(1:10)], df_h)
	saveRDS(unmapped_rps_df, paste0(prefix, '_unmapped_COG_df.rds'))
	# report summary
	summary_value <- colSums(df_h)
	output_title <- paste(names(summary_value), collapse='\t')
	output_value <- paste(summary_value, collapse='\t')
	cat('\n### Title\tgenome_pair\tunmapped_CDS_count\t', output_title, '\n### Result_COGs\t', prefix, '\t', nrow(unmapped_rps_df), '\t', output_value, '\n\n\n', sep='')
}