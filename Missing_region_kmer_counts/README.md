# Investigating putative missing region by comparing k-mer counts between reads and assembly

Related results: Figure 5CD, S1

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_missing_region_kmer_counts.txt](Analysis_missing_region_kmer_counts.txt)

## Dataset:
- `assembly_summary_complete_1ctg_hs.rds` - assembly summary used for fetching source data
- `IdentifyMissingReads.7.covInf_merge_mf95.rds` - results of kmer counts comparison from `IdentifyMissingReads.7.r` for real genomes
- `IdentifyMissingReads_test.5.rds` - results of kmer counts comparison from `IdentifyMissingReads.7.r` for genome subsampling test set
- `CheckM_taxon.list` - checkm taxon_list fetched on August 21 2025 for running CheckM (v1.2.3). Listing available taxonomic-specific marker sets.
- `subsample_CheckM.5.rds` - results from CheckM for genome subsampling test set

## Scripts:
- `IdentifyMissingReads.7.r` - comparing K-mer counts in reads and assembly
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `subsample_genome.1.r` - randomly subsample genome regions to retain a given fraction

## Result figures:
- Figure 5C - Distribution of K-mer counts with fitted
    - [GCF_017356805_covInf.pdf](GCF_017356805_covInf.pdf) - example plot from genome GCF_017356805.1
- Figure 5D - CDFs of estimated fraction missing by sequencing technology
    - [CDF_Kmer.9.pdf](CDF_Kmer.9.pdf) - full range, colored by sequencing technology
- Figure S1 - Genome completeness estimation from our method and CheckM
    - [sub_missing.6.pdf](sub_missing.6.pdf) - all chunk sizes (1st row), chunk size = 1000 nt (2nd row), chunk size = 500 nt (3rd row), our K-mer based method (left column), CheckM completeness (right column)