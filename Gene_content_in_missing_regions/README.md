# Genomic content in regions not shared by high ANI genome pair

Related results: Figure 2CD

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_gene_content_in_missing_regions.txt](Analysis_gene_content_in_missing_regions.txt)

## Dataset:
- `cog-24.fun.tsv` - COG category IDs, group names, colors, and descriptions downloaded from [NCBI COG database](https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab) on November 24 2025
    - `cog24col.tsv` - `cog-24.fun.tsv` with row and column names
- `whole_genome_COG.1.tsv` - result data from `get_COGs.2.r` for whole genome COGs (background)
- `unmap_annot_CDSs.2.tsv` - result data from `unmap_annot.2.r` for regions not shared
- `unmap_annot_COGs.2.tsv` - result data from `get_COGs.3.r` for regions not shared

## Scripts:
- `count_CDS_bases.1.r` - count bases in CDS based on input annotation `*.gff` files
- `unmap_annot.2.r` - get gene annotations for regions not shared
- `get_COGs.2.r` - count bases in regions annotated as COGs for whole genome
- `get_COGs.3.r` - count bases in regions annotated as COGs for regions not shared

## Result figures:
- Figure 2C - CDS vs. non-coding
    - [CDS_NC_frac.4.pdf](CDS_NC_frac.4.pdf) - full range
    - [CDS_NC_frac.5.pdf](CDS_NC_frac.5.pdf) - zoomed-in
- Figure 2D - COG categories enrichment
    - [COG_frac_logodds.7.pdf](COG_frac_logodds.7.pdf)