# Comparing genome sequence similarity (ANI) vs. genome content shared

Related results: Figure 2AB

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_ANI_vs_content_shared.txt](Analysis_ANI_vs_content_shared.txt)


## Dataset:
- `refseq_bacteria_assembly_summary_04112026_complete.rds` - assembly summary of complete genomes
- `refseq_bacteria_assembly_summary_04112026_draft.rds` - assembly summary of draft genomes
- `fastANI_fl3000_maxANI.rds` - pairwise ANI results for complete genomes, highest average ANI pair per genome
- `fastANI_fl3000_maxANI_f999.rds` - `fastANI_fl3000_maxANI` filtered by average ANI > 99.9%
- `fastANI_draftonly_fl3000_maxANI.rds` - pairwise ANI results for draft genomes, highest average ANI pair per genome

## Scripts:
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `calc_asym_X.5.r` - parse nucmer's output as genome content shared result

## Result figures:
- Figure 2A - ANI vs. content shared (complete genomes)
    - [ANI_AF_fastANI3000_nucmer.1.pdf](ANI_AF_fastANI3000_nucmer.1.pdf) - full range
    - [ANI_AF_fastANI3000_nucmer.3.pdf](ANI_AF_fastANI3000_nucmer.3.pdf) - zoomed-in
    - [ANI_AF_fastANI3000_nucmer.5.pdf](ANI_AF_fastANI3000_nucmer.5.pdf) - zoomed-in

- Figure 2B - ANI vs. content shared (draft genomes)
    - [ANI_AF_draftonly_fastANI3000_nucmer.1.pdf](ANI_AF_draftonly_fastANI3000_nucmer.1.pdf) - full range
    - [ANI_AF_draftonly_fastANI3000_nucmer.3.pdf](ANI_AF_draftonly_fastANI3000_nucmer.3.pdf) - zoomed-in