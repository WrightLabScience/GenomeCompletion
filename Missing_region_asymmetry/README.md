# Genome content asymmetry between closely related genome pair (ANI ≥ 99.9%)

Related results: Figure 3ABCD

Notes for how to used the provided dataset to generate the results: [Analysis_missing_region_asymmetry.txt](Analysis_missing_region_asymmetry.txt)

## Dataset:
- `fastANI_fl3000_maxANI_f999.rds` - Pairwise ANI results for complete genomes, highest average ANI pair per genome, filtered by average ANI > 99.9%, merged with genome content shared results sorted by genome sizes per pair

## Result figures:
- Figure 3A - genome region shared skewness, sorted by genome sizes
    - [skew_bySeqTech_nucmer.2.pdf](skew_bySeqTech_nucmer.2.pdf) - full range, colored by sequencing technology
    - [skew_noColor_nucmer.3.pdf](skew_noColor_nucmer.3.pdf) - full range, no color
    - [skew_noColor_nucmer.4.pdf](skew_noColor_nucmer.4.pdf) - zoomed-in
- Figure 3B - CDF of skew by sequencing technology, NOT sorted by genome sizes
    - [CDF_skew_bySeqTech_nucmer.5.pdf](CDF_skew_bySeqTech_nucmer.5.pdf) - full range, colored by sequencing technology
- Figure 3C - genome region shared vs. number of contiguous indel blocks
    - [AF_indelblocks_bySeqTech_nucmer.2.pdf](AF_indelblocks_bySeqTech_nucmer.2.pdf) -  full range, colored by sequencing technology
    - [AF_indelblocks_noColor_nucmer.3.pdf](AF_indelblocks_noColor_nucmer.3.pdf) - full range, no color
    - [AF_indelblocks_noColor_nucmer.4.pdf](AF_indelblocks_noColor_nucmer.4.pdf) - zoomed-in
- Figure 3D - CDF of number of contiguous indel blocks
    - [CDF_indelblocks_bySeqTech_nucmer.3.pdf](CDF_indelblocks_bySeqTech_nucmer.3.pdf) -  full range
    - [CDF_indelblocks_bySeqTech_nucmer.4.pdf](CDF_indelblocks_bySeqTech_nucmer.4.pdf) -  zoomed-in