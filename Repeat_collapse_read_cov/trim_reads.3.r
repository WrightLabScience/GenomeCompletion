#!/usr/bin/env Rscript

# Usage: 
# Rscript trim_reads.3.r --infile_1=fastp_1.fastq.gz --infile_2=fastp_2.fastq.gz

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--infile_1", type="character", default=NA, help="infile_1", metavar="infile_1"),
	make_option("--infile_2", type="character", default=NA, help="infile_2", metavar="infile_2"),
	make_option("--maxAverageError", type="double", default=0.01, help="maxAverageError", metavar="maxAverageError"),
	make_option("--minWidth", type="integer", default=NA, help="minWidth", metavar="minWidth"),
	make_option("--output", type="character", default='', help="output prefix", metavar="output"),
	make_option("--compress", action = "store_true", default=FALSE, help="compress output as .gz"),
	make_option("--fasta", action = "store_true", default=FALSE, help="write fasta")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))


######################## main ########################
# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

# set seed
set.seed(123)

Read_1 <- readQualityScaledDNAStringSet(filepath = opt$infile_1)
if (is.na(opt$minWidth)) {
	minWidth <- max(width(Read_1)) %/% 2
} else {
	minWidth <- opt$minWidth
}

# trim with trimDNA
Read_1 <- TrimDNA(myDNAStringSet = Read_1, leftPatterns = "", rightPatterns = "", type = "sequences", quality = Read_1@quality, maxAverageError = opt$maxAverageError, minWidth = minWidth)

if (is.na(opt$infile_2)) { # SE
	message('# Single end reads.\n')
	if (length(Read_1) == 0) {
		print("No reads remaining after quality filtering.")
		q(save = "no")
	}
	# filter with average quality >= Q20
	quals1a <- mean(as(Read_1@quality, "IntegerList"))
	w1 <- quals1a >= 20
	Read_1 <- Read_1[w1]
	if (length(Read_1) == 0) {
		print("No reads remaining after quality filtering.")
		q(save = "no")
	}
	if (opt$fasta) {
		if (opt$compress) {
			writeXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA.fasta.gz'), compress = TRUE)
		} else {
			writeXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA.fasta'), compress = FALSE)
		}
	} else {
		if (opt$compress) {
			writeQualityScaledXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA.fastq.gz'), compress = TRUE)
		} else {
			writeQualityScaledXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA.fastq'), compress = FALSE)
		}
	}
} else { # PE
	message('# Paired end reads.\n')
	Read_2 <- readQualityScaledDNAStringSet(filepath = opt$infile_2)
	# trim with trimDNA
	Read_2 <- TrimDNA(myDNAStringSet = Read_2, leftPatterns = "", rightPatterns = "", type = "sequences", quality = Read_2@quality, maxAverageError = opt$maxAverageError, minWidth = minWidth)
	# remove unpaired
	partner1 <- do.call(rbind, strsplit(x = names(Read_1), split = " ", fixed = TRUE))
	partner1 <- partner1[, 1L]
	partner2 <- do.call(rbind, strsplit(x = names(Read_2), split = " ", fixed = TRUE))
	partner2 <- partner2[, 1L]
	Read_1 <- Read_1[partner1 %in% partner2]
	Read_2 <- Read_2[partner2 %in% partner1]
	if (length(Read_1) == 0 | length(Read_2) == 0) {
		print("No reads remaining after quality filtering.")
		q(save = "no")
	}
	rm(list = c("partner1", "partner2"))
	# filter with average quality >= Q20
	quals1a <- mean(as(Read_1@quality, "IntegerList"))
	quals2a <- mean(as(Read_2@quality, "IntegerList"))
	w1 <- quals1a >= 20 & quals2a >= 20
	Read_1 <- Read_1[w1]
	Read_2 <- Read_2[w1]
	if (length(Read_1) == 0 | length(Read_2) == 0) {
		print("No reads remaining after quality filtering.")
		q(save = "no")
	}
	if (opt$fasta) {
		if (opt$compress) {
			writeXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA_1.fasta.gz'), compress = TRUE)
			writeXStringSet(x = Read_2, filepath = paste0(opt$output, 'trimDNA_2.fasta.gz'), compress = TRUE)
		} else {
			writeXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA_1.fasta'), compress = FALSE)
			writeXStringSet(x = Read_2, filepath = paste0(opt$output, 'trimDNA_2.fasta'), compress = FALSE)
		}
	} else {
		if (opt$compress) {
			writeQualityScaledXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA_1.fastq.gz'), compress = TRUE)
			writeQualityScaledXStringSet(x = Read_2, filepath = paste0(opt$output, 'trimDNA_2.fastq.gz'), compress = TRUE)
		} else {
			writeQualityScaledXStringSet(x = Read_1, filepath = paste0(opt$output, 'trimDNA_1.fastq'), compress = FALSE)
			writeQualityScaledXStringSet(x = Read_2, filepath = paste0(opt$output, 'trimDNA_2.fastq'), compress = FALSE)
		}
	}
}