#!/usr/bin/env Rscript

# Usage: 
# Rscript IdentifyMissingReads.7.r --genome_file='GCF_037113485.fna.gz' --short_read_file='GCF_037113485_trimDNA.fasta.gz'
# comparing K-mer counts in reads and assembly

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--genome_file", type="character", help="genome_file (Required)"),
	make_option("--short_read_file", type="character", help="short_read_file (Required)"),
	make_option("--MAX_COVERAGE", default=Inf, type="double", help="subsample reads to a maximum coverage (>> 0)"),
	make_option("--FACTR", default=5, type="double", help="factor determining bounds of fitted region (>= 1)"),
	make_option("--MIN_mapped_frac", default=0.5, type="double", help="minimam fraction of reads mapped to genome"),
	make_option("--save_plot", action = "store_true", default=FALSE, help="save plot")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

######################## function ########################



######################## main ########################
# provide files (possibly trimmed, otherwise missing may include adapter)
short_read_file <- scan(text = opt$short_read_file, sep = ",", what = character(), quiet = TRUE)
genomeID <- sub("^(.+)\\.fna.*", "\\1", basename(opt$genome_file))
genome <- readDNAStringSet(opt$genome_file, format="fasta")

MAX_COVERAGE <- opt$MAX_COVERAGE
FACTR <- opt$FACTR
MIN_mapped_frac <- opt$MIN_mapped_frac
K <- as.integer(ceiling(log(sum(width(genome)), 4)))

# compute read k-mer counts
reads <- readDNAStringSet(short_read_file, format="fasta")
if (sum(width(reads))/sum(width(genome)) > MAX_COVERAGE) {
	reads <- sample(reads)
	reads <- reads[seq_len(which.max(cumsum(as.numeric(width(reads))) >= sum(width(genome))*MAX_COVERAGE))]
}
counts_reads <- oligonucleotideFrequency(reads, K, simplify.as="collapse")
reads <- reverseComplement(reads)
counts_reads <- counts_reads + oligonucleotideFrequency(reads, K, simplify.as="collapse")

# estimate coverage
N <- min(1e6L, length(reads))
index <- IndexSeqs(c(genome, reverseComplement(genome)), K=K, step=K, processors=NULL)
hits <- SearchIndex(sample(reads, N), index, perPatternLimit=1, processors=NULL)
mapped_fraction <- nrow(hits)/N
cat('# mapped_fraction:', mapped_fraction, '\n')
if (mapped_fraction < MIN_mapped_frac) {
	cat("## Insufficient mapped reads:", genomeID, mapped_fraction, '\n')
	quit(save = "no", status = 0)
}
total_length <- sum(width(reads))
coverage <- mapped_fraction*total_length/sum(width(genome))
cat('# coverage:', coverage, '\n')

# find the region to fit
t <- tabulate(counts_reads)
d <- diff(log(t))
w <- which(head(d, coverage/2) < 0)
w <- w[which.max(t[w] <= max(t[w])/FACTR)]:w[which.max(t[w] <= max(FACTR, FACTR*min(t[w])))]
w <- w[t[w] > 0]
y <- t[w]

# initialize optimality function
SSE <- function(params) {
	mult <- params[1L]
	size <- params[2L]
	mean <- params[3L]
	Y <- mult*dnbinom(w, size=size, mu=mean)
	sum((log(y) - log(Y))^2) # sum of squared error
}

# grid search for initial parameters
ini <- expand.grid(mult=max(y)*2^(-5:25),
	size=2^seq(-10, 10, 0.25),
	mean=2^seq(-10, 10, 0.25))
val <- numeric(nrow(ini))
for (i in seq_len(nrow(ini)))
	val[i] <- SSE(unlist(ini[i,]))

# fit truncated NB distribution
best <- Inf
i <- c(head(order(val), 100),
	tapply(seq_along(val), ini$mult, function(x) x[which.min(val[x])]))
for (j in i) {
	temp <- try(optim(unlist(ini[j,]),
			SSE,
			lower=2^c(-5, -10, -10),
			method="L-BFGS-B"),
		silent=TRUE)
	if (is(temp, "try-error")) # fall back to Nelder-Mead
		temp <- optim(unlist(ini[j,]), SSE)
	if (temp$value < best) {
		o <- temp
		best <- temp$value
	}
}
cat('# fitted SSE:', o$value, '\n')
cat('# relative to bounds:', genomeID, t(as.matrix(sapply(ini, range)))/o$par, '\n')

# determine the threshold of significance
threshold <- o$par[1L]*pnbinom(seq_along(t), size=o$par[2L], mu=o$par[3L], lower.tail=FALSE)
threshold <- which.max(threshold < 1L)
cat('# threshold:', threshold, '\n')

# plot observed counts and fit
if (opt$save_plot) {
	pdf(file = paste0(genomeID, '_cov', as.character(MAX_COVERAGE), '.pdf'))
	plot(t,
		ylim=c(1, max(t)),
		xlab=paste0(K, "-mer count"),
		ylab="Occurrences in reads",
		log="xy",
		type="l",
		main=paste0('IdentifyMissingReads_v7 (', genomeID, ')'))
	lines(w, t[w], col="red")
	lines(w[1L]:threshold,
		o$par[1L]*dnbinom(w[1L]:threshold, size=o$par[2L], mu=o$par[3L]),
		col="green")
	abline(v=threshold, lty=2)
	dev.off()
}

# compute genome k-mer counts
counts_genome <- oligonucleotideFrequency(genome, K, simplify.as="collapse") +
	oligonucleotideFrequency(reverseComplement(genome), K, simplify.as="collapse")

# identify missing k-mers
missing <- sum(counts_reads >= threshold &
	counts_genome == 0L)
cat('# bases missing:', missing, '\n')
cat('# fraction missing:', missing/sum(width(genome)), '\n')

output_title <- paste('### Title_Kmer', 'genomeID', 'MAX_COVERAGE', 'K', 'mapped_fraction', 'coverage', 'w', 'par_mult', 'par_size', 'par_mean', 'fitted_SSE', 'threshold', 'genome_size', 'bases_missing', 'frac_missing', sep='\t')
output_value <- paste('### Result_Kmer', genomeID, MAX_COVERAGE, K, mapped_fraction, coverage, 
                      ifelse(length(w) == 0, NA, w),
                      o$par[1L], o$par[2L], o$par[3L], o$value,
                      ifelse(is.na(threshold), NA, threshold),
                      sum(width(genome)),
                      ifelse(is.na(missing), NA, missing),
                      ifelse(is.na(missing/sum(width(genome))), NA, missing/sum(width(genome))),
                      sep='\t')
cat('\n', output_title, '\n', output_value, '\n\n', sep='')