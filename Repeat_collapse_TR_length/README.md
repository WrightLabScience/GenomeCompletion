# Investigating repeat collapse by comparing exact tandem repeat total lengths distributions

Related results: Figure 4D

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_Repeat_collapse_TR_length.txt](Analysis_Repeat_collapse_TR_length.txt)

## Dataset:
- `assembly_summary_balanced.rds` - assembly summary for `fetch_source_data.3.r`, genomes balanced by species taxid and sequencing technology
- `TR_total_len.2.zip` - zipped `TR_total_len.2/`, results from `TR_total_len.1.r`, tandem repeat total lengths of genomes

## Scripts:
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `TR_total_len.1.r` - get exact tandem repeat total lengths in genomes

## Result figures:
- Figure 4D - CDFs of the "total length" (not unit length) distribution of identical repeats by sequencing technologies
    - [CDF_TR_lens.3.pdf](CDF_TR_lens.3.pdf) - full range
    - [CDF_TR_lens.4.pdf](CDF_TR_lens.4.pdf) - zoomed-in