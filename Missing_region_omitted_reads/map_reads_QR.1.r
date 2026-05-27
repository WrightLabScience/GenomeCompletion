#!/usr/bin/env Rscript

# Usage: 
# Rscript map_reads_QR.1.r --reads='reads.fastq' --genomeQ='GCF_000834965.fna.gz' --genomeR='GCF_000834965.fna.gz' --K=8 --proc=8
# map reads to genome assembly

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--reads", type="character", help="reads file, pattern, source (Required)", metavar="reads"),
	make_option("--genomeQ", type="character", help="genomeQ file, subject, source (Required)", metavar="genomeQ"),
	make_option("--genomeR", type="character", help="genomeR file, subject, partner (Required)", metavar="genomeR"),
	make_option("--K", type="integer", default=8, help="K-mer length", metavar="K"),
	make_option("--minScore", type="integer", default=NA, help="minScore for SearchIndex", metavar="minScore"),
	make_option("--cutoff_frac_mapped", type="double", default=0.95, help="cutoff of fraction of genomeQ mapped"),
	make_option("--proc", type="integer", default=1, help="the number of processors to use", metavar="proc")
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
rotate_circular_genome <- function(genome, rotate_bp=301) {
    cat('# rotate_bp:', rotate_bp, '\n')
    DNAStringSet(c(genome[[1]], subseq(genome[[1]], 1, rotate_bp))) # append a short fragment of start region to the end
}

SearchIndex_AlignPairs <- function(pattern, subject, index, minScore, proc, chunk_size=25000, rc=FALSE) {
	read_lens <- width(pattern)
	n <- length(pattern)
	chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
	
	aligned_list <- lapply(chunks, function(idx) {
		message('### SearchIndex chunk ', idx[1], '-', idx[length(idx)], ' ...')
		chunk_pattern <- if (rc) reverseComplement(pattern[idx]) else pattern[idx]
		hits <- SearchIndex(
			pattern=chunk_pattern, invertedIndex=index, minScore=minScore,
			perPatternLimit=0, perSubjectLimit=0,
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
		message('# frac_mapped:', names(subject)[contig_i], '\t', width(subject)[contig_i], '\t', frac_mapped)
	}
	return(list(aligned_merge = aligned_merge, frac_mapped = frac_mapped))
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
max_read_len <- max(width(pattern_q_reads))
subject_q_genome <- tryCatch(
	readDNAStringSet(opt$genomeQ),
	error = function(e) {
		cat("# Cannot read genomeQ file:", opt$genomeQ, "\n")
		quit(save = "no", status = 1)
	}
)
cat('# GenomeQ length:', width(subject_q_genome), '\n')

subject_r_genome <- tryCatch(
	readDNAStringSet(opt$genomeR),
	error = function(e) {
		cat("# Cannot read genomeR file:", opt$genomeR, "\n")
		quit(save = "no", status = 1)
	}
)
cat('# GenomeR length:', width(subject_r_genome), '\n')
output_prefix <- paste0(sub("(\\.[^.]+)+$", "", basename(opt$genomeQ)), '_', sub("(\\.[^.]+)+$", "", basename(opt$genomeR)))

# rotate/append genome
subject_q_genome <- rotate_circular_genome(subject_q_genome, rotate_bp=max_read_len)
subject_r_genome <- rotate_circular_genome(subject_r_genome, rotate_bp=max_read_len)

# map reads to both genomes
message('### map_reads q ...')
result_q <- map_reads(pattern=pattern_q_reads, subject=subject_q_genome, K=opt$K, minScore=opt$minScore, proc=opt$proc)
aligned_merge_q    <- result_q$aligned_merge
frac_mapped_q      <- result_q$frac_mapped
if (frac_mapped_q < opt$cutoff_frac_mapped) {
	cat("# frac_mapped_q too low:", opt$genomeQ, frac_mapped_q, "\n")
	quit(save = "no", status = 1)
}
rm(subject_q_genome)
saveRDS(aligned_merge_q, file = paste0(output_prefix, "_aligned_q.rds"))

message('### map_reads r ...')
result_r <- map_reads(pattern=pattern_q_reads, subject=subject_r_genome, K=opt$K, minScore=opt$minScore, proc=opt$proc)
aligned_merge_r    <- result_r$aligned_merge
frac_mapped_r      <- result_r$frac_mapped
readIDs_q <- seq_along(pattern_q_reads)
rm(subject_r_genome, pattern_q_reads)
saveRDS(aligned_merge_r, file = paste0(output_prefix, "_aligned_r.rds"))

# discard any reads with any matches less than the cutoff (0.99) and greater than some minimum threshold (0.9). These hits create ambiguity.
reads2rm <- unique(c(aligned_merge_q[aligned_merge_q$PID<0.99,'Pattern'], aligned_merge_r[aligned_merge_r$PID<0.99,'Pattern']))
aligned_merge_q <- aligned_merge_q[!aligned_merge_q$Pattern %in% reads2rm,]
aligned_merge_r <- aligned_merge_r[!aligned_merge_r$Pattern %in% reads2rm,]
readIDs_map2q <- aligned_merge_q$Pattern
readIDs_map2r <- aligned_merge_r$Pattern

# report reads map to R only
message('# list of q reads map to r only:')
print(aligned_merge_r[aligned_merge_r$Pattern %in% setdiff(unique(readIDs_map2r), unique(readIDs_map2q)), ])
rm(aligned_merge_q, aligned_merge_r)
gc()

# record the difference in number of hits
#readIDs_q <- seq_along(pattern_q_reads)
readIDs_q <- readIDs_q[!readIDs_q %in% reads2rm]
n_q <- integer(length(readIDs_q))
n_r <- integer(length(readIDs_q))
names(n_q) <- names(n_r) <- readIDs_q
hits_q <- table(readIDs_map2q)
hits_r <- table(readIDs_map2r)
n_q[names(hits_q)] <- hits_q
n_r[names(hits_r)] <- hits_r
hit_table <- table(n_q, n_r)
hit_table <- as.data.frame.matrix(hit_table)
message('# record the difference in number of hits:')
print(hit_table)
saveRDS(hit_table, file = paste0(output_prefix,"_hit_table.rds"))

# report number of reads (Venn diagram)
readIDs_map2q <- unique(readIDs_map2q)
readIDs_map2r <- unique(readIDs_map2r)
QnR <- length(intersect(readIDs_map2q, readIDs_map2r)) # [1:,1:]
Qonly <- length(setdiff(readIDs_map2q, readIDs_map2r)) # [1:,0]
Ronly <- length(setdiff(readIDs_map2r, readIDs_map2q)) # [0,1:]
non <- length(setdiff(readIDs_q, union(readIDs_map2q, readIDs_map2r))) # [0,0]
discarded <- length(reads2rm)
frac_Ronly <- Ronly / (QnR + Qonly + Ronly)
# message('# reads map to q & r:', length(intersect(readIDs_map2q, readIDs_map2r)))
# message('# reads map to q but not r:', length(setdiff(readIDs_map2q, readIDs_map2r)))
# message('# reads map to r but not q (!):', length(setdiff(readIDs_map2r, readIDs_map2q)))
# message('# reads map to neither [1,1]:', length(setdiff(readIDs_q, union(readIDs_map2q, readIDs_map2r))))
# message('# discarded reads with ambiguity:', length(reads2rm))

# output
output_title <- paste('# Title:', 'genomeQ', 'genomeR', 'frac_mapped_Q', 'frac_mapped_R', 'QnR', 'Qonly', 'Ronly', 'non','discarded', 'frac_Ronly', sep='\t')
output_value <- paste('# Results:', sub("(\\.[^.]+)+$", "", basename(opt$genomeQ)), sub("(\\.[^.]+)+$", "", basename(opt$genomeR)), frac_mapped_q, frac_mapped_r, QnR, Qonly, Ronly, non, discarded, frac_Ronly, sep='\t')
cat('\n', output_title, '\n', output_value, '\n\n\n', sep='')