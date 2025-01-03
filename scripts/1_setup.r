# Install Bioconductor Package Manager if Not Already Installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install and Load Required Bioconductor Packages for ATAC-seq Analysis
required_bioc_packages <- c(
    "Rsubread", "Rsamtools", "GenomicRanges", 
    "GenomicFeatures", "GenomicAlignments", 
    "rtracklayer", "ComplexHeatmap", "stringi", "matrixStats", "reticulate"
)


########## Troubleshooting Notes for Package Installation ##########
# 1. Delete any lock files from the R library folder:
#    Example: /path/to/R/library/00LOCK-<package_name>
# 2. Issues installing `stringi`? Use Conda for installation on Linux systems:
#    `conda install -c conda-forge r-stringi`
#    Reference: https://stackoverflow.com/questions/67260518/error-installing-the-stringi-package-on-r-4-on-linux-ubuntu
# 3. `BiocManager` is the preferred tool for managing Bioconductor packages. It replaces older methods such as `biocLite`.
