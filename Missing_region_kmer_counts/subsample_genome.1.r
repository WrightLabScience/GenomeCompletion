# Usage: 
# Rscript subsample_genome.1.r --in_fasta=sequence.fasta --out_fasta=out.fasta --fraction=0.9 --max_fragments=100

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--in_fasta", type="character", default="sequence.fasta", help="sequence.fasta", metavar="sequence.fasta"),
	make_option("--out_fasta", type="character", default="out.fasta", help="out.fasta", metavar="out.fasta"),
	make_option("--fraction", type="double", default=0.9, help="fraction retain", metavar="fraction"),
	make_option("--chunk_size", type="integer", default=1000, help="chunk size for sampling", metavar="chunk_size")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################

genome <- readDNAStringSet(opt$in_fasta)
full_seq <- unlist(genome) # flatten into one contiguous sequence
genome_len <- length(full_seq)

# split into chunks
starts <- seq(1, genome_len, by=opt$chunk_size)
ends <- starts + opt$chunk_size - 1
ends[length(ends)] <- genome_len
chunks <- IRanges(start=starts, end=ends)

# Number of chunks to keep
keep_n <- round(opt$fraction * length(chunks))

# Randomly pick chunks
keep_idx <- sort(sample(seq_along(chunks), keep_n))
keep_ranges <- chunks[keep_idx]
keep_ranges <- reduce(keep_ranges)

# Extract kept contigs
contigs <- DNAStringSet(NULL)
for (out_frag_i in seq_along(keep_ranges)) {
	contigs <- c(contigs, DNAStringSet(subseq(full_seq, start=start(keep_ranges)[out_frag_i], end=end(keep_ranges)[out_frag_i])))
}
names(contigs) <- paste0("contig", seq_along(contigs))
writeXStringSet(contigs, opt$out_fasta)

out_text <- paste0('n_contigs = ', length(contigs), ', frac_subsampled = ', opt$fraction, ', len_ori = ', genome_len, ', len_subsampled = ', sum(width(contigs)), '\n')
cat(out_text)

if (abs(sum(width(contigs))/genome_len-opt$fraction) > 0.05) {
	cat('##### !!!!! Fraction Error !!!!!\n', sum(width(contigs))/genome_len, opt$fraction, '\n')
}