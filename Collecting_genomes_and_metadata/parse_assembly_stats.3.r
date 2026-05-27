#!/usr/bin/env Rscript

# Usage: 
# Rscript parse_assembly_stats.3.r --inputXML=esummary.xml

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--inputXML", type="character", default="esummary.xml", help="inputXML", metavar="inputXML"),
    make_option("--output", type="character", default="out.txt", help="output", metavar="output")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

library(xml2)
library(XML)

######################## function ########################


######################## main ########################
esummary_h <- read_xml(opt$inputXML)
AssemblyAccession <- xml_text(xml_find_first(esummary_h, ".//AssemblyAccession"))
Taxid <- xml_text(xml_find_first(esummary_h, ".//Taxid"))
SpeciesTaxid <- xml_text(xml_find_first(esummary_h, ".//SpeciesTaxid"))
SpeciesName <- xml_text(xml_find_first(esummary_h, ".//SpeciesName"))
AssemblyStatus <- xml_text(xml_find_first(esummary_h, ".//AssemblyStatus"))
BioprojectAccn <- xml_text(xml_find_first(esummary_h, ".//BioprojectAccn"))
BioSampleAccn <- xml_text(xml_find_first(esummary_h, ".//BioSampleAccn"))
Coverage <- xml_text(xml_find_first(esummary_h, ".//Coverage"))

output_text <- c(
    AssemblyAccession,
    Taxid,
    SpeciesTaxid,
    SpeciesName,
    AssemblyStatus,
    BioprojectAccn,
    BioSampleAccn,
    Coverage
    )

write(output_text, file = opt$output, ncolumns = length(output_text), append = FALSE, sep = "\t")