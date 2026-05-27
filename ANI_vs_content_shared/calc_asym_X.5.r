#!/usr/bin/env Rscript

# Usage: 
# Rscript calc_asym_X.4.r --genus='genus' --genomeQ='GCF_000834965.fna.gz' --genomeR='GCF_000833165.fna.gz' --coords=coords.out.visual --fragLen=3000
# Calculate degree of asymmety and Kuiper's test X for plotting

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--genomeQ_ID", type="character", help="genomeQ ID, query (Required)"),
	make_option("--genomeR_ID", type="character", help="genomeR ID, ref (Required)"),
	make_option("--coords", type="character", help="nucmer.coords (Required)"),
	make_option("--genus", type="character", default='NA', help="genus"),
	make_option("--min_gap_size", type="integer", default=50, help="minimum gap size")
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
### replicates ###
replicates <- function(v_gaps) {
	N <- 1e4L # number of replicates
	pvals <- numeric(N)
	for (i in seq_len(N)) {
		v_gaps <- sample(c(F, T), 100, prob=c(0.99, 0.01), replace=TRUE)
		X <- sum(abs(range(cumsum(v_gaps)/sum(v_gaps) - seq_along(v_gaps)/length(v_gaps))))*sqrt(sum(v_gaps))
		pvals[i] <- (2*(4*X^2 - 1) - 8*X/(3*sqrt(sum(v_gaps)))*(4*X^2 - 3))*exp(-2*X^2)
	}
	cat('# pvals range:', range(pvals, na.rm=TRUE), '\n')
	cat('# mean pvals < 0.05:', mean(pvals < 0.05, na.rm=TRUE), '\n') # want this to be 0.05
}

######################## main ########################
# parameters
min_gap_size <- opt$min_gap_size

# inputs
genomeQ_ID <- opt$genomeQ_ID
genomeR_ID <- opt$genomeR_ID
coords <- read.table(opt$coords, skip = 4, header = FALSE, col.names = c("start_R","end_R","start_Q","end_Q","fraglen_R","fraglen_Q","perc_id","LENR","LENQ","COVR","COVQ","genomeR","genomeQ"))
lenQ <- coords$LENQ[1]
lenR <- coords$LENR[1]

# (1) How much of Genome Q is covered — merge overlapping query intervals
ir_query <- IRanges(start = pmin(coords$start_Q, coords$end_Q), end = pmax(coords$start_Q, coords$end_Q))
ir_query <- reduce(ir_query)
cov_query <- sum(width(ir_query))
perc_cov_query <- round(cov_query * 100 / lenQ, 4)
cat('\n# Number of continous matched ranges:', length(ir_query), '\n')

# (2) Count contiguous indel blocks — gaps between merged query blocks
ir_gaps_query <- gaps(ir_query, start = 1, end = lenQ)
cat('\n# Number of gaps before size filtering:', length(ir_gaps_query), '\n')
saveRDS(ir_gaps_query, paste0(genomeQ_ID, '_', genomeR_ID, '.rds'))

ir_gaps_query_f <- ir_gaps_query[width(ir_gaps_query)>=min_gap_size]
cat('\n# Number of gaps after size filtering:', length(ir_gaps_query_f), '\n')

cat('\n# Blocks of genomeQ missing in genomeR:\n')
df <- as.data.frame(ir_gaps_query)
df$keep <- df$width>=min_gap_size
print(df)

### save v_gaps and num of block events
if (length(ir_gaps_query_f)>0) { # if there are gaps
	num_indels <- length(ir_gaps_query_f)
	cov <- coverage(ir_gaps_query_f, width = lenQ)
	v_gaps <- as.logical(cov)
} else { # no gaps
	num_indels <- 0
	v_gaps <- logical(lenQ)
}
cat('\n# Length of genomeQ:', length(v_gaps), '\n')
cat('\n# Fraction of genomeQ missing in genomeR:', mean(v_gaps), '\n')

# plot synteny map
pdf( paste0(genomeQ_ID, '_', genomeR_ID, ".pdf", sep="") )
plot(NULL, xlim = c(0,lenQ), ylim = c(1,lenR), xlab = paste('qry:', genomeQ_ID), ylab = paste('ref:', genomeR_ID), main = paste0('nucmer synteny map\n(perc_covQ=',perc_cov_query, '; n_block_missing=', num_indels, ')'))
segments(x0 = coords$start_Q, y0 = coords$start_R, x1 = coords$end_Q, y1 = coords$end_R, col = ifelse(coords$start_Q<coords$end_Q, 'black','red'))
if (length(ir_gaps_query)>0) { # before filtering
	for(i in seq(length(ir_gaps_query))) {
		rect(start(ir_gaps_query)[i], 0, end(ir_gaps_query)[i], lenR, col=rgb(1,0,0,0.2), border=NA)
	}
}
if (length(ir_gaps_query_f)>0) { # after filtering
	for(i in seq(length(ir_gaps_query_f))) {
		rect(start(ir_gaps_query_f)[i], 0, end(ir_gaps_query_f)[i], lenR, col=rgb(0,0,1,0.2), border=NA)
	}
}
dev.off()

################ Kuiper's test ################
# X = test statistic is V*sqrt(N), where V is sum of absolute deviations and N is the number of events:
X_Kuiper <- sum(abs(range(cumsum(v_gaps)/sum(v_gaps) - seq_along(v_gaps)/length(v_gaps))))*sqrt(sum(v_gaps))
# for x > 1.2, the p-value (P) can reasonably be approximated as:
# NOTE: need to confirm X > 1.2 (otherwise define P as 1), and P is bounded 0 to 1
Pval <- ifelse(X_Kuiper > 1.2, (2*(4*X_Kuiper^2 - 1) - 8*X_Kuiper/(3*sqrt(sum(v_gaps)))*(4*X_Kuiper^2 - 3))*exp(-2*X_Kuiper^2), 1)
logPval <- ifelse(X_Kuiper > 1.2, log(2*(4*X_Kuiper^2 - 1) - 8*X_Kuiper/(3*sqrt(sum(v_gaps)))*(4*X_Kuiper^2 - 3)) + (-2*X_Kuiper^2), 0)
logPval <- ifelse(is.na(logPval), -Inf, logPval)
replicates(v_gaps)

# output
output_title <- paste('# Title:', 'genus', 'genomeQ', 'genomeR', 'covQ', 'lenQ', 'Frac_missing', 'X_Kuiper', 'Pval', 'logPval', 'num_indels', sep='\t')
output_value <- paste('# Results:', opt$genus, genomeQ_ID, genomeR_ID, cov_query, lenQ, mean(v_gaps), X_Kuiper, Pval, logPval, num_indels, sep='\t')
cat('\n', output_title, '\n', output_value, '\n\n\n', sep='')