# Investigating repeat collapse by read coverage

Related results: Figure 4ABC

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_repeat_collapse_read_cov.txt](Analysis_repeat_collapse_read_cov.txt)

## Dataset:
- `assembly_summary_complete_1ctg_hs.rds` - assembly summary for `get_SRA_info.2.r`
- `SRA_info.1.tsv` - results from `get_SRA_info.2.r`
- `SRA_info.2.tsv` - `SRA_info.1` filtered by paired and illumina
- `frac_high_cov.1.tsv` - results from `map_reads.6.r`
- `frac_high_cov.1_merge_f95_balanced.rds` - `frac_high_cov.1` filtered by fraction of genome mapped, and balanced by species taxid and sequence technology

## Scripts:
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `get_SRA_info.2.r` - retrieve SRA info (run IDs) from NCBI
- `trim_reads.3.r` - trim reads by quality score and pairing
- `map_reads.6.r` - map reads to the source genome, get regions with statistically high coverage

## Result figures:
- Figure 4A - example of high cov region in GCF_001677195.1
    - [Cov_spike.GCF_001677195.1.pdf](Cov_spike.GCF_001677195.1.pdf)
- Figure 4B - example of fitting a negative binomial distribution to left-truncated read coverage histogram in GCF_011405515.1
    - [NB_MLE.GCF_011405515.pdf](NB_MLE.GCF_011405515.pdf)
- Figure 4C - CDFs of fraction high coverage by sequencing technology
    - [CDF_frac_high_cov.5.pdf](CDF_frac_high_cov.5.pdf) - full range
    - [CDF_frac_high_cov.6.pdf](CDF_frac_high_cov.6.pdf) - zoomed-in