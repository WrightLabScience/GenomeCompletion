# Collecting genome assemblies and corresponding metadata

Related results: Figure 1B

Notes for how to used the provided dataset and scripts to generate the results: [Analysis_genome_assembly_dataset.txt](Analysis_genome_assembly_dataset.txt)


## Dataset:
- `refseq_bacteria_assembly_summary_04112026.txt.zip` - NCBI RefSeq bacteria assembly summary (downloaded on April 11 2026)
- `CheckM_report_prokaryotes_04112026.txt.zip` - CheckM completeness/contamination report (downloaded on April 11 2026)
- `SeqTech_class.1.rds` - result of sequencing technology classification
- `refseq_bacteria_assembly_summary_04112026.rds` - parsed assembly summary with metadata appended

## Scripts:
- `fetch_source_data.3.r` - download source data (i.e., genome assemblies `*.fna`, metadata `*_assembly_stats.txt`, annotation `*.gff`) from NCBI
- `parse_assembly_stats.1.r` - parse metadata from `*_assembly_stats.txt`
- `parse_assembly_stats.3.r` - parse metadata from `esummary.xml`

## Result figures:
- Figure 1B left panel - CDFs of CheckM completeness
    - [CDF_checkm_comp.1.pdf](CDF_checkm_comp.1.pdf) - full range
    - [CDF_checkm_comp.2.pdf](CDF_checkm_comp.2.pdf) - zoomed-in
    - [CDF_checkm_comp.3.pdf](CDF_checkm_comp.3.pdf) - zoomed-in

- Figure 1B right panel: CDFs of CheckM contamination
    - [CDF_checkm_cont.1.pdf](CDF_checkm_cont.1.pdf) - full range
    - [CDF_checkm_cont.2.pdf](CDF_checkm_cont.2.pdf) - zoomed-in
    - [CDF_checkm_cont.3.pdf](CDF_checkm_cont.3.pdf) - zoomed-in