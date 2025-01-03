srun -n 4 -t 3:00:00 --pty -p mhgcp --mem=2000 macs3 callpeak \
    -t ~/ruuhsu/bam/IMS_dTAG13_ATAC_3_S24_R1.PE2SE.nodup.repair.sorted.bam \
    -f BAMPE --outdir ~/ruuhsu/bam/results \
    --name IMS_dTAG13_ATAC_3_S24 -g hs

# BAM to BED Conversion
# macs3 randsample -i file.bam -f BAMPE -p 100 -o output.bed
