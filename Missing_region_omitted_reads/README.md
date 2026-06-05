# Investigating putative missing region by mapping reads to source and high ANI partner genomes

Related results: Figure 5A

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_missing_region_omitted_reads.txt](Analysis_missing_region_omitted_reads.txt)

## Dataset:
- `df_fastANI_maxANI_f999.rds` - list of closely related genome pairs (ANI >= 99.9%), merged with sequencing technology and species taxid information
- `df_fastANI_maxANI_f999_mapQR.rds` - results of `map_reads_QR.1.r` merged with `df_fastANI_maxANI_f999`
- `cat_fastANI.2.tsv` - results from `fastANI`, for finding closely related partner genomes

## Scripts:
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `map_reads_QR.1.r` - map reads to genome assembly, both source genome (genomeQ) and closely related partner genome (genomeR)

## Result figures:
- Figure 5A - CDFs of fraction reads mapped to partner genome only (ommitted reads)
    - [CDF_frac_Ronly.3.pdf](CDF_frac_Ronly.3.pdf) - full range
    - [CDF_frac_Ronly.5.pdf](CDF_frac_Ronly.5.pdf) - zoomed-in