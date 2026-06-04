# Investigating gene features in putative missing regions with omitted reads

Related results: Figure 5B

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_features_in_omitted_regions.txt](Analysis_features_in_omitted_regions.txt)


## Dataset:
- `assembly_summary_complete_1ctg_hs.rds` - assembly summary used for fetching source data
- `feature_bases_merge.4.rds` - result data from count feature bases

## Scripts:
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `count_feature_bases.3.r` - count feature bases based on input annotation `*.gff` files
- `find_Ronly.2.r` - find regions in the closely related partner genome mapped by omitted reads from the source genome

## Result figures:
- Figure 5B - Odds of feature bases
    - [feature_frac_logodds.7.pdf](feature_frac_logodds.7.pdf) - with 99% confidence interval error bars
    - [feature_frac_logodds.8.pdf](feature_frac_logodds.8.pdf) - p < 0.01 significance shown as "*"