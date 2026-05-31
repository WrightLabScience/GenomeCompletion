# Usage: 
# Rscript fetch_source_data.3.r assembly_summary.rds GCF_000026945 ftp_path gbff fna gff assembly_stats

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
assembly_summary <- ARGS[1]
seqID <- ARGS[2]
ftp_colname <- ARGS[3]
file_types <- ARGS[4:length(ARGS)]

print(ARGS)

# load libraries

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## function ########################
fetch_data <- function(seqID, out_ext, FTP_dir, FTP_ext, attempt_n=5) {
	out_file <- paste(seqID, out_ext, sep='')
	FTPadd <- paste(FTP_dir, gsub('.*\\/(.+)\\/$', '\\1', FTP_dir), FTP_ext,sep='')
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
# download source data
assembly_summary <- readRDS(assembly_summary)

cat(seqID, '\n')
FTP_dir <- assembly_summary[seqID, ftp_colname]
if (FTP_dir!='' | FTP_dir!='na') {
	if (is.element('gbff', file_types)) {
		fetch_data(seqID, out_ext='.gbff.gz', FTP_dir, FTP_ext='_genomic.gbff.gz', attempt_n=5)
	}
	if (is.element('fna', file_types)) {
		fetch_data(seqID, out_ext='.fna.gz', FTP_dir, FTP_ext='_genomic.fna.gz', attempt_n=5)
	}
	if (is.element('gff', file_types)) {
		fetch_data(seqID, out_ext='.gff.gz', FTP_dir, FTP_ext='_genomic.gff.gz', attempt_n=5)
	}
	if (is.element('faa', file_types)) {
		fetch_data(seqID, out_ext='.faa.gz', FTP_dir, FTP_ext='_protein.faa.gz', attempt_n=5)
	}
    if (is.element('assembly_stats', file_types)) {
		fetch_data(seqID, out_ext='_assembly_stats.txt', FTP_dir, FTP_ext='_assembly_stats.txt', attempt_n=5)
	}
}