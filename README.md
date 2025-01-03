# ATACseq Analysis Notes

### Rotation Project: Chromatin Accessibility Changes in Leukemia NPM1 Mutant Cells

This repository documents the analysis pipeline for an ATAC-seq experiment investigating chromatin accessibility changes across the genome in leukemia NPM1 mutant cells. The study compares three control groups (DMSO) and three experimental groups (dTAG).

### Overview

The pipeline includes the following key steps:
1. **Quality Control (QC):** Ensuring data quality and integrity for downstream analysis.
2. **Peak Calling:** Identifying regions of open chromatin.
3. **Differential Accessibility Analysis:** Determining genomic regions with significant changes in accessibility between groups.
4. **Annotation:** Mapping identified peaks to genes and genomic features.

### References
- **QC Guide:** [ATACseqQC Bioconductor Vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html)  
- **Pipeline Workflow:** [ATAC-seq Workshop by Sean Davis](https://seandavi.github.io/AtacSeqWorkshop/articles/Workflow.html)

### Repository Contents
- **scripts:** Code for data preprocessing, alignment, and analysis.
- **results:** Key findings, including peak files and visualizations.

Feel free to explore the repository and adapt the pipeline for your ATAC-seq projects!
