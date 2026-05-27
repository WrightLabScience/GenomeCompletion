#!/usr/bin/env Rscript

# Usage: 
# Rscript get_SRA_info.1.r --input=assembly_summary.rds --output=job_01.tsv

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--input", type="character", default=NA, help="input", metavar="input"),
	make_option("--start_i", type="integer", default=1, help="begin row_i", metavar="start_i"),
	make_option("--end_i", type="integer", default=100, help="end row_i", metavar="end_i"),
	make_option("--output", type="character", default="out.tsv", help="output", metavar="output")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))


######################## function ########################
# sra_search <- entrez_search(db = "sra", term = biosample_h)
safe_entrez_search <- function(db, term, max_attempts = 5) {
	for (attempt_i in seq_len(max_attempts)) {
		sra_search <- tryCatch(
			{
				entrez_search(db = db, term = term)
			},
			#if an error occurs, tell me the error
			error=function(e) {
				message("Error occurred: ", conditionMessage(e))
				message("Retry search attempt ", attempt_i)
				return(NULL)
			},
			#if a warning occurs, tell me the warning
			warning=function(w) {
				message('Warning occurred: ', conditionMessage(w))
			}
		)
		if (!is.null(sra_search)) {
			return(sra_search)
		}
		Sys.sleep(2^attempt_i)
		message("Retry attempt ", attempt_i)
	}
	return(NULL)
}

safe_entrez_fetch <- function(db, id, rettype, retmode, max_attempts = 5) {
	for (attempt_i in seq_len(max_attempts)) {
		runinfo <- tryCatch(
			{
				entrez_fetch(db = db, id = id, rettype=rettype, retmode=retmode)
			},
			#if an error occurs, tell me the error
			error=function(e) {
				message("Error occurred: ", conditionMessage(e))
				message("Retry search attempt ", attempt_i)
				return(NULL)
			},
			#if a warning occurs, tell me the warning
			warning=function(w) {
				message('Warning occurred: ', conditionMessage(w))
			}
		)
		if (!is.null(runinfo)) {
			return(runinfo)
		}
		Sys.sleep(2^attempt_i)
		message("Retry attempt ", attempt_i)
	}
	return(NULL)
}

######################## main ########################
library(rentrez)
# API key
set_entrez_key("55a0d4f7439395d7c765ecefe3deb00d1c09")

df_input <- readRDS(opt$input)

system.time({
	df_list <- list()
	done_j <- 1
	for (row_i in seq(from = opt$start_i, to = opt$end_i)) {
		genomeID_h <- df_input[row_i,'genomeID']
		biosample_h <- df_input[row_i,'biosample']
		sra_search <- safe_entrez_search(db = "sra", term = biosample_h)
		if (!is.null(sra_search) && sra_search$count > 0) {
			runinfo <- safe_entrez_fetch(db = "sra", id = paste(sra_search$ids, collapse = ","), rettype = "runinfo", retmode = "text")
			if (!is.null(runinfo)) {
				runinfo_table <- read.csv(text = runinfo, stringsAsFactors = FALSE)
				if (nrow(runinfo_table) == 0) {
					message(row_i, '\t', genomeID_h, ' FAILED!!! empty runinfo\n')
					next
				}
				runinfo_table$genomeID <- genomeID_h
				runinfo_table$biosample <- biosample_h
				df_list[[done_j]] <- runinfo_table[,c('genomeID', 'biosample', 'Run', 'Platform', 'LibraryLayout', 'Model')]
				done_j=done_j+1
				message(row_i, '\t', genomeID_h, ' done.\n')
			} else {
				message(row_i, '\t', genomeID_h, ' FAILED!!! entrez_fetch\n')
				next
			}
		} else {
			message(row_i, '\t', genomeID_h, ' FAILED!!! entrez_search\n')
			next
		}
		Sys.sleep(0.12)
		if (row_i %% 500 == 0) {
			tmp_df <- do.call(rbind, df_list)
			write.table(tmp_df, file=paste0(opt$output,".checkpoint"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
		}
	}
	df_output <- do.call(rbind, df_list)
})

write.table(df_output, file = opt$output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)