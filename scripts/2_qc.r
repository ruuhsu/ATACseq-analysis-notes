### Load required packages from 1_setup.r ###
library(Rsamtools)
library(ggplot2)
library(magrittr)
library(GenomicAlignments)

# Index creation for the reference genome
genome <- "~/ATAC_Workshop/ATAC_Data/ATAC_Reference/hg19_Genome.fa"
indexForSubread <- gsub("\\.fa$", "", genome)
# buildindex(indexForSubread, genome, indexSplit = FALSE)

# Example BAM files
read1 <- "~/ruuhsu/bam/IMS_DMSO_ATAC_1_S19_R1.PE2SE.nodup.repair.sorted.bam"

# Mapping rate calculation
library(Rsubread)
bam_files <- list(
    IMS_dTAG13_S22 = "~/ruuhsu/bam/IMS_dTAG13_ATAC_1_S22_R1.PE2SE.nodup.repair.sorted.bam",
    IMS_dTAG13_S23 = "~/ruuhsu/bam/IMS_dTAG13_ATAC_2_S23_R1.PE2SE.nodup.repair.sorted.bam",
    IMS_dTAG13_S24 = "~/ruuhsu/bam/IMS_dTAG13_ATAC_3_S24_R1.PE2SE.nodup.repair.sorted.bam"
)

mapped_rates <- lapply(bam_files, propmapped)
print(mapped_rates)

# Distribution of mapped reads
for (file in names(bam_files)) {
    idxstatsBam(bam_files[[file]]) %>% 
        ggplot(aes(seqnames, mapped, fill = seqnames)) + 
        geom_bar(stat = "identity") + 
        coord_flip() + 
        ggtitle(paste("Mapped Reads Distribution:", file)) +
        theme_minimal()
}

# Reading mapped reads
atac_reads <- lapply(bam_files, function(file) {
    readGAlignmentPairs(
        file, 
        param = ScanBamParam(
            mapqFilter = 1, 
            flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), 
            what = c("qname", "mapq", "isize")
        )
    )
})

# Coverage calculations
coverage_results <- lapply(atac_reads, coverage)
print(names(coverage_results))
