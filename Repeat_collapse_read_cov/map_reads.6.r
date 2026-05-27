#!/usr/bin/env Rscript

# Usage: 
# Rscript map_reads.6.r --reads='reads.fastq' --genomeQ='GCF_000834965.fna.gz' --K=8 --proc=8
# map reads to genome assembly

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--reads", type="character", help="reads file, pattern, query (Required)", metavar="reads"),
	make_option("--genomeQ", type="character", help="genomeQ file, subject, query (Required)", metavar="genomeQ"),
	make_option("--K", type="integer", default=8, help="K-mer length", metavar="K"),
	make_option("--minScore", type="integer", default=NA, help="minScore for SearchIndex", metavar="minScore"),
	make_option("--proc", type="integer", default=1, help="the number of processors to use", metavar="proc"),
	make_option("--cutoff_frac_mapped", type="double", default=0.95, help="cutoff of fraction of genomeQ mapped"),
	make_option("--cutoff_frac", type="double", default=0.5, help="set fraction of mode as left censoring cutoff"),
	make_option("--alpha", type="double", default=0.05, help="cutoff for significant p-val"),
	make_option("--cutoff_len", type="integer", default=50, help="cutoff for high coverage region length"),
	make_option("--cutoff_fold", type="integer", default=2, help="cutoff for high coverage region fold")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")
suppressMessages(library(dplyr))
options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## function ########################
SearchIndex_AlignPairs <- function(pattern, subject, index, minScore, proc, chunk_size=25000, rc=FALSE) {
	read_lens <- width(pattern)
	n <- length(pattern)
	chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
	
	aligned_list <- lapply(chunks, function(idx) {
		message('### SearchIndex chunk ', idx[1], '-', idx[length(idx)], ' ...')
		chunk_pattern <- if (rc) reverseComplement(pattern[idx]) else pattern[idx]
		hits <- SearchIndex(
			pattern=chunk_pattern, invertedIndex=index, minScore=minScore,
			perPatternLimit=0, perSubjectLimit=0, # keeping everything now and save best-scored later
			maskRepeats=FALSE, maskLCRs=FALSE, proc=proc
		)
		if (nrow(hits) == 0) return(NULL)
		
		message('### AlignPairs chunk ...')
		aligned <- AlignPairs(
			pattern=chunk_pattern, subject=subject, pairs=hits, proc=proc,
			perfectMatch=5, misMatch=-15, gapOpening=-20, gapExtension=-10
		)
		
		# fix Pattern index to global index
		aligned$Pattern <- idx[aligned$Pattern]
		
		aligned$PID <- aligned$Matches / read_lens[aligned$Pattern]
		aligned$gap_length <- abs(aligned$SubjectEnd - aligned$SubjectStart) - abs(aligned$PatternEnd - aligned$PatternStart)
		aligned <- aligned[((aligned$PID > 0.9) & (aligned$gap_length < 10)), ]
		gc()
		return(aligned[, c("Pattern","PatternStart","PatternEnd","Subject","SubjectStart","SubjectEnd","Matches","Mismatches","AlignmentLength","Score","PID")])
	})
	
	do.call(rbind, Filter(Negate(is.null), aligned_list))
}

map_reads <- function(pattern, subject, K, minScore, proc) {
	# IndexSeqs
	message('### IndexSeqs ...')
	index <- IndexSeqs(subject=subject, K=K, step=K, maskRepeats=FALSE, maskLCRs=FALSE, maskNumerous=FALSE)
	# SearchIndex: both strands
	message('### SearchIndex_AlignPairs same ...')
	aligned_same <- SearchIndex_AlignPairs(pattern=pattern, subject=subject, index=index, minScore=minScore, proc=proc, rc=FALSE)
	saveRDS(aligned_same, file='tmp_same')
	rm(aligned_same); gc()

	message('### SearchIndex_AlignPairs rc ...')
	aligned_rc <- SearchIndex_AlignPairs(pattern=pattern, subject=subject, index=index, minScore=minScore, proc=proc, rc=TRUE)

	message('### rbind merge aligned ...')
	aligned_same <- readRDS('tmp_same')
	aligned_merge <- rbind(aligned_same, aligned_rc)
	rm(aligned_same, aligned_rc); gc()
	for (contig_i in seq_along(subject)) {
		idx <- aligned_merge$Subject == contig_i # row index for contig_i
		cat('contig_i:', contig_i, '\n')
		if (sum(idx)>0) {
			ir <- IRanges(start = pmin(aligned_merge$SubjectStart[idx], aligned_merge$SubjectEnd[idx]), end = pmax(aligned_merge$SubjectStart[idx], aligned_merge$SubjectEnd[idx]))
			ir <- reduce(ir)  # merge overlaps
			frac_mapped <- sum(width(ir)) / width(subject)[contig_i]
		} else {
			frac_mapped <- 0
		}
		cat('# frac_mapped:', names(subject)[contig_i], '\t', width(subject)[contig_i], '\t', frac_mapped, '\n')
	}
	return(list(aligned_merge = aligned_merge, frac_mapped = frac_mapped))
}

# MLE optimize parameters for NB for left truncated read coverage
fit_truncated_nb <- function(x, cutoff, init_mu=NULL, init_size=NULL) {
	# Keep only truncated data
	x_trunc <- x[x >= cutoff]
	n <- length(x_trunc)

	if (n == 0) stop("No data above cutoff.")

	# Initial values
	if (is.null(init_mu)) init_mu <- mean(x_trunc)
	if (is.null(init_size)) {
		v <- var(x_trunc)
		init_size <- ifelse(v > init_mu, init_mu^2 / (v - init_mu), 10)
	}
	
	# Negative log-likelihood
	nll <- function(par) {
		mu <- exp(par[1])
		size <- exp(par[2])

		if (!is.finite(mu) || !is.finite(size) || mu <= 0 || size <= 0) {
			return(Inf)  # return large penalty instead of NaN
		}
		
		# log PMF for observed values
		ll_obs <- sum(dnbinom(x_trunc, size=size, mu=mu, log=TRUE))
		
		# log normalization constant
		log_tail <- pnbinom(cutoff - 1, size=size, mu=mu, lower.tail=FALSE, log.p=TRUE)
		
		if (!is.finite(log_tail)) return(Inf)
		
		ll <- ll_obs - n * log_tail
		
		return(-ll)
	}
	
	fit <- optim(
		par = log(c(init_mu, init_size)),
		fn = nll,
		method = "L-BFGS-B"
	)
	
	list(
		mu = exp(fit$par[1]),
		size = exp(fit$par[2]),
		logLik = -fit$value,
		convergence = fit$convergence
	)
}


######################## main ########################
pattern_q_reads <- tryCatch(
	readDNAStringSet(opt$reads, format = "fasta"),
	error = function(e) {
		message("FASTA read failed, trying FASTQ...")
		tryCatch(
			readDNAStringSet(opt$reads, format = "fastq"),
			error = function(e2) {
				cat("# Cannot read reads file:", opt$reads, "\n")
				quit(save = "no", status = 1)
			}
		)
	}
)

subject_q_genome <- tryCatch(
	readDNAStringSet(opt$genomeQ),
	error = function(e) {
		cat("# Cannot read genomeQ file:", opt$genomeQ, "\n")
		quit(save = "no", status = 1)
	}
)
cat('# Genome length:', width(subject_q_genome), '\n')

output_prefix <- paste0(sub("(\\.[^.]+)+$", "", basename(opt$genomeQ)))

# map reads to genome
message('### map_reads q ...')
result_q <- map_reads(pattern=pattern_q_reads, subject=subject_q_genome, K=opt$K, minScore=opt$minScore, proc=opt$proc)
aligned_merge_q    <- result_q$aligned_merge
frac_mapped_q      <- result_q$frac_mapped
if (frac_mapped_q < opt$cutoff_frac_mapped) {
	cat("# frac_mapped_q too low:", opt$genomeQ, frac_mapped_q, "\n")
	quit(save = "no", status = 1)
}
rm(subject_q_genome, pattern_q_reads); gc()
saveRDS(aligned_merge_q, file = paste0(output_prefix, "_aligned_q.rds"))

##### calculate read coverage #####
# filtering: keep best hits and PID>=0.99
cat('# dim before filter:', dim(aligned_merge_q), '\n')
aligned_q_best <- aligned_merge_q %>%
	group_by(Pattern) %>%
	filter(Score == max(Score)) %>%
	ungroup()
aligned_q_best <- aligned_q_best[aligned_q_best$PID>=0.99,]
cat('# dim after filter:', dim(aligned_q_best), '\n')
if (nrow(aligned_q_best) < 1) {
	cat('# No reads mapped after filtering:', output_prefix, '\n')
	quit(save = "no", status = 0)  # status=0 = clean exit, status=1 = error exit
}
rm(aligned_merge_q); gc()

# normalize read coverage: weighted coverage
# weight = 1/num of hits for that read
hits_per_read <- table(aligned_q_best$Pattern) # num of hits per read
aligned_q_best$weight <- 1 / hits_per_read[as.character(aligned_q_best$Pattern)] # assign weight
# add coverage
ir <- IRanges(start = aligned_q_best$SubjectStart, end = aligned_q_best$SubjectEnd)
cov <- coverage(ir, weight = aligned_q_best$weight)
rm(aligned_q_best)
cov_vec <- as.numeric(cov)

######## use negative binomial distribution to calculate p-value ########
# calculate cutoff for left censoring
cov_smooth <- stats::filter(cov_vec, rep(1/5, 5), sides = 2) # smoothing w/ window size=5bp
mode_cov <- as.numeric(names(sort(-table(round(na.omit(cov_smooth)))))[1])
rm(cov_smooth)
threshold_trunc    <- max(1, round(mode_cov * opt$cutoff_frac))
threshold_foldcov  <- mode_cov * opt$cutoff_fold

# fit the model
fit <- fit_truncated_nb(x=round(cov_vec), cutoff=threshold_trunc, init_mu=NULL, init_size=NULL)

# compute p-values
pvals <- pnbinom(cov_vec - 1, size = fit$size, mu = fit$mu, lower.tail = FALSE)
padj <- p.adjust(pvals, method = "BH")
cat('# padj:\n')
print(summary(padj))

# high cov regions
width_genomeSeq <- length(cov_vec)

ir_high_cov <- reduce(IRanges(start = which(padj < opt$alpha), width = 1))
ir_high_cov <- ir_high_cov[width(ir_high_cov) >= opt$cutoff_len]
frac_high_cov = sum(width(ir_high_cov)) / width_genomeSeq
print(ir_high_cov)

ir_ltfoldcov <- reduce(IRanges(which(cov_vec >= threshold_foldcov)))
print(ir_ltfoldcov)
if (length(ir_high_cov) > 0 & length(ir_ltfoldcov) > 0) {
	ir_fold_high_cov <- IRanges::intersect(ir_high_cov, ir_ltfoldcov)
} else {
	ir_fold_high_cov <- IRanges()
}
frac_fold_high_cov <- sum(width(ir_fold_high_cov)) / width_genomeSeq
save(cov_vec, ir_high_cov, ir_fold_high_cov, padj, file=paste0(output_prefix,'_highcov.RData'), compress='xz')

# output
output_title <- paste('# Title:', 'genomeID', 'frac_high_cov', 'n_high_cov', 'frac_fold_high_cov', 'n_fold_high_cov', 'genome_length', 'frac_mapped', sep='\t')
output_value <- paste('# Results:', output_prefix, frac_high_cov, length(ir_high_cov), frac_fold_high_cov, length(ir_fold_high_cov), width_genomeSeq, frac_mapped_q, sep='\t')
cat('\n', output_title, '\n', output_value, '\n\n\n', sep='')