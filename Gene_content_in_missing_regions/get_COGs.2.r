#!/usr/bin/env Rscript

# Usage: 
# Rscript get_COGs.2.r --proteinIDs='GCF_030291915_proteinIDs.tsv' --assembly_summary=assembly_summary_1ctg_s50.rds
# get protein fasta from protein id table

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--proteinIDs", type="character", help="genomeID_proteinIDs.txt (Required)"),
	make_option("--faa_dir", type="character", default=".", help="source data dir for faa.gz", metavar="dirpath"),
    make_option("--faa_ext", type="character", default=".faa.gz", help="faa filename extension", metavar="faa_ext"),
	make_option("--db", type="character", default="~/scrdir/cdd_11242025/Cog", help="db path for rpsblast", metavar="db"),
	make_option("--cog_tab", type="character", default="cog24col.tsv", help="COG table", metavar="COG"),
	make_option("--assembly_summary", type="character", default="assembly_summary.rds", help="assembly_summary.rds file with ftp_path for genomes", metavar="assembly_summary")
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
proteinIDs <- scan(opt$proteinIDs, character(), quote = "")
genomeID <- sub('_proteinIDs.txt','',opt$proteinIDs)

## check and prepare source data
if (file.exists(paste0(gsub("/$", "", opt$faa_dir),'/', genomeID, opt$faa_ext))) {
	cat('# genomeID faa exists:', genomeID, '\n')
} else {
	assembly_summary <- readRDS(opt$assembly_summary)
	FTP_dir <- assembly_summary[genomeID, 'ftp_path']
	fetch_data(genomeID, out_ext=opt$faa_ext, FTP_dir, FTP_ext='_protein.faa.gz', attempt_n=5)
}

# faa
faa <- readAAStringSet(paste0(gsub("/$", "", opt$faa_dir),'/', genomeID, opt$faa_ext))
proteinID_df <- data.frame(proteinID = sub(" .*", "", names(faa)), product = sub("\\s\\[[^]]+\\]$", "", sub(".*MULTISPECIES:\\s", "", names(faa))), stringsAsFactors = FALSE)
proteinID_df$product <- mapply(function(id, product) {sub(paste0("^", id, "\\s+"), "", product)}, proteinID_df$proteinID, proteinID_df$product)
# export selected faa seqs
sel_faa <- faa[proteinID_df$proteinID %in% proteinIDs]
names(sel_faa) <- sub(" .*", "", names(sel_faa))
writeXStringSet(sel_faa, 'sel.faa')
# merge product table
sel_df <- merge(data.frame(proteinID = proteinIDs), proteinID_df, by = "proteinID", all.x = TRUE)
#cat(length(proteinIDs), length(faa), dim(sel_df), '\n')

# run rpsblast and process result in bash
RPSBLAST <- paste0('rpsblast',
                   ' -query sel.faa',
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
rps <- read.table('results_rps_f.tsv', sep="\t", quote="", header=FALSE)
colnames(rps) <- c('proteinID','evalue','bitscore','cog_info')
rps2 <- rps %>%
  extract(cog_info,
          into = c("COG","gene","description","categories"),
          regex = "^(COG[0-9]+), ([^,]+), ([^\\[]+) \\[([^]]+)\\].*",
          remove=TRUE)
#cat('#rps2 dim: ', dim(rps2), '\n')
for (categoryID_h in rownames(cog24.fun)) {
    description_h <- cog24.fun[categoryID_h, 'description']
    rps2[,categoryID_h] <- grepl(description_h, rps2$categories, ignore.case = TRUE)
}
#cat('#rps2 dim: ', dim(rps2), '\n')
sel_rps_df <- merge(sel_df, rps2, by='proteinID', all.x=TRUE)
#cat(dim(rps2), dim(sel_rps_df),'\n')
sel_rps_df[,c(9:34)][is.na(sel_rps_df[,c(9:34)])] <- 0
df_h <- sel_rps_df[,c(9:34)]
df_h$noCOG <- ifelse(rowSums(df_h)<1, 1, 0) # for no blast hits
df_h <- df_h/rowSums(df_h) # normalize
sel_rps_df <- cbind(sel_rps_df[,c(1:8)], df_h)
saveRDS(sel_rps_df, paste0(genomeID, '_sel_COG_df.rds'))

# report summary
summary_value <- colSums(df_h)
output_title <- paste(names(summary_value), collapse='\t')
output_value <- paste(summary_value, collapse='\t')
cat('\n### Title\tgenomeID\ttotal_CDS_count\t', output_title, '\n### Result\t', genomeID, '\t', length(proteinIDs), '\t', output_value, '\n\n\n', sep='')