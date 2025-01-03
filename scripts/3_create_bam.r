#The coverage method for genomic ranges calculates, for each base, the number of overlapping features. In the case of a BAM file from ATAC-Seq converted to a GAlignmentPairs object, the coverage gives us an idea of the extent to which reads pile up to form peaks.
cvg = coverage(IMS_dTAG13_S22)
class(cvg)
names(cvg)
cvg$chr7

# length(atacReads)
IMS_dTAG13_S22_atacReads

#Retrieving insert sizes
IMS_dTAG13_S22_atacReads_1 <- GenomicAlignments::first(IMS_dTAG13_S22_atacReads)
IMS_dTAG13_S22_insertSizes <- abs(elementMetadata(IMS_dTAG13_S22_atacReads_1)$isize)
head(IMS_dTAG13_S22_insertSizes)


#Plotting insert sizes
library(magrittr)
library(dplyr)
library(ggplot2)
IMS_dTAG13_S22_fragLenPlot <- table(IMS_dTAG13_S22_insertSizes) %>% data.frame %>% rename(IMS_dTAG13_S22_insertSizes = IMS_dTAG13_S22_insertSizes, 
                                                                                      Count = Freq) %>% mutate(IMS_dTAG13_S22_insertSizes = as.numeric(as.vector(IMS_dTAG13_S22_insertSizes)), 
                                                                                                               Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = IMS_dTAG13_S22_insertSizes, y = Count)) + 
  geom_line() 

IMS_dTAG13_S22_fragLenPlot + theme_bw()
IMS_dTAG13_S22_fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()

#open, mono- and di-nucleosome profiles
#what's the cut-off for open, mono- and di- nucleosome???
IMS_dTAG13_S22_fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 
                                                                                                            437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw()


IMS_dTAG13_S22_fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 
                                                                                          247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()


#retriving tss

#1.tss on reference genome
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19_TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), fix = "start", 1)
hg19_TSSs

#nucleosome free
#BiocManager::install("soGGi")
library(soGGi)

IMS_dTAG13_S22_nucFree <- regionPlot(bamFile = IMS_dTAG13_S22, testRanges = hg19_TSSs, style = "point", 
                                   format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100, 
                                   forceFragment = 50) #takes time to run

IMS_dTAG13_S22_nucFree


# Mononucleosome
IMS_dTAG13_S22_monoNuc <- regionPlot(bamFile = IMS_dTAG13_S22, testRanges = hg19_TSSs, style = "point", 
                                   format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240, 
                                   forceFragment = 80)

# Dinucleosome
IMS_dTAG13_S22_diNuc <- regionPlot(bamFile = IMS_dTAG13_S22, testRanges = hg19_TSSs, style = "point", 
                                 format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437, 
                                 forceFragment = 160)

IMS_dTAG13_S22_nucFree_gL <- IMS_dTAG13_S22_nucFree 
IMS_dTAG13_S22_monoNuc_gL <- IMS_dTAG13_S22_monoNuc 
IMS_dTAG13_S22_diNuc_gL <- IMS_dTAG13_S22_diNuc
save(IMS_dTAG13_S22_monoNuc_gL,IMS_dTAG13_S22_nucFree_gL,IMS_dTAG13_S22_diNuc_gL,file='/storage/goodell/projects/ruuhsu/bam/IMS_dTAG13_S22_gL_soGGiResults.RData')


#Plotting ATAC-seq signal of TSSs (Plotting open, mono- and di-nucleosome signal profiles)
library(soGGi)
load(file = '/storage/goodell/projects/ruuhsu/bam/IMS_dTAG13_S22_gL_soGGiResults.RData')
plotRegion(IMS_dTAG13_S22_nucFree_gL, outliers = 0.01)
plotRegion(IMS_dTAG13_S22_diNuc_gL, outliers = 0.01)
plotRegion(IMS_dTAG13_S22_monoNuc_gL, outliers = 0.01)


#create BAM for insertSizes
IMS_dTAG13_S22_atacReads_Open <- IMS_dTAG13_S22_atacReads[IMS_dTAG13_S22_insertSizes<100, ]
IMS_dTAG13_S22_atacReads_Open
IMS_dTAG13_S22_atacReads_MonoNuc <- IMS_dTAG13_S22_atacReads[IMS_dTAG13_S22_insertSizes > 180 & IMS_dTAG13_S22_insertSizes < 240, ]
IMS_dTAG13_S22_atacReads_MonoNuc
IMS_dTAG13_S22_atacReads_diNuc <- IMS_dTAG13_S22_atacReads[IMS_dTAG13_S22_insertSizes > 315 & IMS_dTAG13_S22_insertSizes < 437, ]
IMS_dTAG13_S22_atacReads_diNuc

#Creating BAM files split by insert sizes.
IMS_dTAG13_S22_openRegionBam <- gsub("\\.bam", "_openRegions\\.bam", IMS_DMSO_S19)
IMS_dTAG13_S22_monoNucBam <- gsub("\\.bam", "_monoNuc\\.bam", IMS_DMSO_S19)
IMS_dTAG13_S22_diNucBam <- gsub("\\.bam", "_diNuc\\.bam", IMS_DMSO_S19)

library(rtracklayer)
export(IMS_dTAG13_S22_atacReads_Open, IMS_dTAG13_S22_openRegionBam, format = "bam")
export(IMS_dTAG13_S22_atacReads_MonoNuc, IIMS_dTAG13_S22_monoNucBam, format ="bam")
export(IMS_dTAG13_S22_atacReads_diNuc, IMS_dTAG13_S22_diNucBam,format = "bam")

#Creating an open region bigWig.
IMS_dTAG13_S22_openRegionBigWig <- gsub("\\.bam", "_openRegions\\.bw", IMS_dTAG13_S22)
IMS_dTAG13_S22_openRegionRPMBigWig <- gsub("\\.bam", "_openRegionsRPM\\.bw", IMS_dTAG13_S22)
IMS_dTAG13_S22_atacFragments_Open <- granges(IMS_dTAG13_S22_atacReads_Open)
export.bw(coverage(IMS_dTAG13_S22_atacFragments_Open), IMS_dTAG13_S22_openRegionBigWig)
