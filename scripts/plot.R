library(RColorBrewer)
# library(magrittr)
# library(tidyverse)
# setwd(paste0('/share/users/yujw/20190711_sub_species_identifing/output/admixture/'))
library(readr)
k_5 <- as.matrix(read_table2(paste0(snakemake@input[["qfile"]]), 
                             col_names = FALSE))
sample_info<- read_table2(snakemake@input[["sampleinfo"]], col_names = F)
used_sample<- t(read_delim(paste0(snakemake@input[["samplename"]]), col_names = F, delim = ','))
rownames(sample_info)<-sample_info$X1
sample_info<- sample_info[as.vector(used_sample),]
rownames(k_5)<- used_sample
sample_info$X2<-factor(sample_info$X2)
pdf(file=snakemake@output[["plots"]], width=10, height=5)
barplot(t(k_5[order(sample_info$X2),]), border = NA, space = 0, col = brewer.pal(5, 'Set2'),las=2)
dev.off()