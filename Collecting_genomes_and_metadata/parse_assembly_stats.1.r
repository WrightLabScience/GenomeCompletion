#!/usr/bin/Rscript

# parse_assembly_stats.1.r

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read *_assembly_stats.txt and output tsv

# v1 2025/02/18

# Usage:
# Rscript parse_assembly_stats.1.r /Users/shc167/Documents/project/proj03/proj03.07/ _assembly_stats.txt /Users/shc167/Documents/project/proj03/proj03.07/assembly_stats_out.tsv


# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
in_dir <- ARGS[1]
in_filename_ext <- ARGS[2] # _assembly_stats.txt
out_file <- ARGS[3]

print(ARGS)


######################## function ########################



######################## main ########################
# get names
filename_pattern <- paste0('*', in_filename_ext)
filenames <- list.files(path=in_dir, pattern = filename_pattern)

colnames <- c(
	"Result_Header",
	"Assembly_name",
	"Organism_name",
	"Isolate",
	"Taxid",
	"BioSample",
	"BioProject",
	"Submitter",
	"Date",
	"Assembly_type",
	"Release_type",
	"Assembly_level",
	"Genome_representation",
	"Assembly_method",
	"Genome_coverage",
	"Sequencing_technology",
	"Excluded_from_RefSeq",
	"GenBank_assembly_accession",
	"RefSeq_assembly_accession",
	"RefSeq_assembly_and_GenBank_assemblies_identical",
	"GC_content")
cat(colnames, sep='\t')
cat('\n')
write(colnames, file = out_file, ncolumns = length(colnames), append = FALSE, sep = "\t")

# initiate values
Assembly_name <- NA
Organism_name <- NA
Isolate <- NA
Taxid <- NA
BioSample <- NA
BioProject <- NA
Submitter <- NA
Date <- NA
Assembly_type <- NA
Release_type <- NA
Assembly_level <- NA
Genome_representation <- NA
Assembly_method <- NA
Genome_coverage <- NA
Sequencing_technology <- NA
Excluded_from_RefSeq <- NA
GenBank_assembly_accession <- NA
RefSeq_assembly_accession <- NA
RefSeq_assembly_and_GenBank_assemblies_identical <- NA
GC_content <- NA

# opening result file
n_files <- length(filenames)
for (i in seq_along(filenames)) {
	filename_h <- filenames[i]
	lines <- readLines(filename_h)
	lines <- iconv(lines, "UTF-8", "UTF-8", sub = "byte")
	lines <- trimws(lines)
	for (line in lines) {
		if (startsWith(line, "# Assembly name:")) {
			Assembly_name <- gsub("# Assembly name:\\s+", "", line)
		} else if (startsWith(line, "# Organism name:")) {
			Organism_name <- gsub("# Organism name:\\s+", "", line)
		} else if (startsWith(line, "# Isolate:")) {
			Isolate <- gsub("# Isolate:\\s+", "", line)
		} else if (startsWith(line, "# Taxid:")) {
			Taxid <- gsub("# Taxid:\\s+", "", line)
		} else if (startsWith(line, "# BioSample:")) {
			BioSample <- gsub("# BioSample:\\s+", "", line)
		} else if (startsWith(line, "# BioProject:")) {
			BioProject <- gsub("# BioProject:\\s+", "", line)
		} else if (startsWith(line, "# Submitter:")) {
			Submitter <- gsub("# Submitter:\\s+", "", line)
		} else if (startsWith(line, "# Date:")) {
			Date <- gsub("# Date:\\s+", "", line)
		} else if (startsWith(line, "# Assembly type:")) {
			Assembly_type <- gsub("# Assembly type:\\s+", "", line)
		} else if (startsWith(line, "# Release type:")) {
			Release_type <- gsub("# Release type:\\s+", "", line)
		} else if (startsWith(line, "# Assembly level:")) {
			Assembly_level <- gsub("# Assembly level:\\s+", "", line)
		} else if (startsWith(line, "# Genome representation:")) {
			Genome_representation <- gsub("# Genome representation:\\s+", "", line)
		} else if (startsWith(line, "# Assembly method:")) {
			Assembly_method <- gsub("# Assembly method:\\s+", "", line)
		} else if (startsWith(line, "# Genome coverage:")) {
			Genome_coverage <- gsub("# Genome coverage:\\s+", "", line)
		} else if (startsWith(line, "# Sequencing technology:")) {
			Sequencing_technology <- gsub("# Sequencing technology:\\s+", "", line)
		} else if (startsWith(line, "# Excluded from RefSeq:")) {
			Excluded_from_RefSeq <- gsub("# Excluded from RefSeq:\\s+", "", line)
		} else if (startsWith(line, "# GenBank assembly accession:")) {
			GenBank_assembly_accession <- gsub("# GenBank assembly accession:\\s+", "", line)
		} else if (startsWith(line, "# RefSeq assembly accession:")) {
			RefSeq_assembly_accession <- gsub("# RefSeq assembly accession:\\s+", "", line)
		} else if (startsWith(line, "# RefSeq assembly and GenBank assemblies identical:")) {
			RefSeq_assembly_and_GenBank_assemblies_identical <- gsub("# RefSeq assembly and GenBank assemblies identical:\\s+", "", line)
		} else if (startsWith(line, "all\tall\tall\tall\tgc-perc")) {
			GC_content <- gsub("all\tall\tall\tall\tgc-perc\t", "", line)
		}
	}
	output <- c("Result_output", Assembly_name, Organism_name, Isolate, Taxid, BioSample, BioProject, Submitter, Date, Assembly_type, Release_type, Assembly_level, Genome_representation, Assembly_method, Genome_coverage, Sequencing_technology, Excluded_from_RefSeq, GenBank_assembly_accession, RefSeq_assembly_accession, RefSeq_assembly_and_GenBank_assemblies_identical, GC_content)
	# write results
	write(output, file = out_file, ncolumns = length(output), append = TRUE, sep = "\t")
	# progress report every 100 files
	if (i %% 100 == 0 || i == n_files) {
		cat(sprintf("[%d/%d] %.1f%% - last: %s\n", i, n_files, 100 * i / n_files, filename_h))
	}
	# reset values for next file
	Assembly_name <- NA; Organism_name <- NA; Isolate <- NA
	Taxid <- NA; BioSample <- NA; BioProject <- NA
	Submitter <- NA; Date <- NA; Assembly_type <- NA
	Release_type <- NA; Assembly_level <- NA; Genome_representation <- NA
	Assembly_method <- NA; Genome_coverage <- NA; Sequencing_technology <- NA
	Excluded_from_RefSeq <- NA; GenBank_assembly_accession <- NA
	RefSeq_assembly_accession <- NA
	RefSeq_assembly_and_GenBank_assemblies_identical <- NA; GC_content <- NA
}

cat("\nn_file_in = ", n_files, '\n', sep='')