#!/usr/bin/env Rscript

# Usage: 
# Rscript TR_total_len.1.r --genome='GCF_000834965.fna.gz'
# find exact TR and get distribution of element length
# ref: https://support.bioconductor.org/p/68996/#69062 by Hervé Pagès

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--genome", type="character", help="genome file(Required)"),
	make_option("--min_length", type="integer", default=24, help="min.length"),
	make_option("--min_period", type="integer", default=2, help="min.period"),
	make_option("--max_period", type="integer", default=100, help="max.period")
); 


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(Biostrings))
packageVersion("Biostrings")
options(timeout=999999999) # default is 60 sec, will fail if the genome is too big


######################## function ########################
## Find all *exact* tandem repeats with a period equal to or a divisor of
## 'period.multiple'.
.findTandemRepeats0 <- function(subject, period.multiple=2,
										 include.period1=FALSE, min.length=opt$min_length)
{
	if (!isSingleNumber(period.multiple))
		stop("'period.multiple' must be a single integer")
	if (!is.integer(period.multiple))
		period.multiple <- as.integer(period.multiple)
	if (period.multiple < 2L)
		stop("'period.multiple' must be >= 2")
	if (!isTRUEorFALSE(include.period1))
		stop("'include.period1' must be TRUE or FALSE")
	if (!isSingleNumber(min.length))
		stop("'min.length' must be a single integer")
	if (!is.integer(min.length))
		min.length <- as.integer(min.length)
	if (min.length < 12L)
		stop("'min.length' must be >= 12")
	s1 <- subseq(subject, start=1L+period.multiple)
	s2 <- subseq(subject, end=-1L-period.multiple)
	
	# fix: use IRanges() directly instead of coercing to NormalIRanges
	match_vec <- as.raw(s1) == as.raw(s2)
	ir <- IRanges(start = which(diff(c(FALSE, match_vec, FALSE)) == 1),
				  end   = which(diff(c(FALSE, match_vec, FALSE)) == -1) - 1)
	
	ir <- ir[width(ir) >= period.multiple]
	end(ir) <- end(ir) + period.multiple
	ir <- ir[width(ir) >= min.length]
	repeats <- Views(subject, ir)
	
	af <- alphabetFrequency(repeats, baseOnly=TRUE)
	ok <- af[ , "other"] == 0L  # has no IUPAC ambiguities
	if (!include.period1)
		ok <- ok & rowSums(af[ , DNA_BASES] != 0L) >= 2L
	repeats[which(ok)]
}

## Find all *exact* tandem repeats of period <= 12.
findTandemRepeats <- function(subject, include.period1=FALSE, min.period=7, max.period=12)
{
	trs_list <- lapply(min.period:max.period,
		function(period.multiple)
			ranges(.findTandemRepeats0(subject,
									   period.multiple,
									   include.period1))
	)
	trs <- sort(unique(do.call("c", trs_list)))
	Views(subject, trs)
}


######################## main ########################
# read genomic DNA fna
genome <- readDNAStringSet(opt$genome)

system.time({
	v_TR_lengths <- c()
	for (i in seq_along(genome)) {
		cat("Processing contig", i, "of", length(genome), "\n")
		TRs <- findTandemRepeats(genome[[i]], include.period1=TRUE, min.period=opt$min_period, max.period=opt$max_period)
		v_TR_lengths <- c(v_TR_lengths, width(TRs))
	}
})
print(summary(v_TR_lengths))
print(length(v_TR_lengths))
saveRDS(v_TR_lengths, paste0(sub("(\\.[^.]+)+$", "", basename(opt$genome)), '.rds'))