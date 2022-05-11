## Gene expression quantification
## Author:Zhao Peng (pengzhao@nwafu.edu.cn)
## Northwest A&F University
## This script briefly summarizes the codes in the analysis process, for reference only.

## Gene expression quantification
kallisto index -i ./IWGSC_V1.1_HC_LC ./IWGSC_v1.1_HC_LC_transcripts.fasta
kallisto quant -i ./IWGSC_V1.1_HC_LC -t 6 -o ./result/sample/ sample_1.fq.gz sample_2.fq.gz
## Summarized the expression levels
library("tximport")
samples <- read.table("sample_list.txt",sep="\t",header = T)
tx2gene <- read.table("gene_trans_HC_LC.txt",sep="\t",header = T)
dir <- "./"
files <- file.path(dir, "kallisto", samples$Sample, "abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene,countsFromAbundance = c("lengthScaledTPM"))
write.table(txi.kallisto$abundance,file = "TPM_HC_LC.txt",sep = "\t",quote = F)
