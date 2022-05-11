## Population genetic analysis
## Author:Zhao Peng (pengzhao@nwafu.edu.cn)
## Northwest A&F University
## This script briefly summarizes the codes in the analysis process, for reference only.

## Construct a phylogenies tree
vcftools --vcf snp.vcf --recode --recode-INFO-all --stdout --max-missing 0.5 --maf 0.05 > tree.vcf
raxmlHPC-PTHREADS-SSE3 -s tree.phy -f a -m GTRGAMMA -p 12346 -x 12346 -# 100 -T 24 -n m5f5
## PCA
plink --allow-extra-chr --out temp --recode --vcf SNP.vcf
plink --allow-extra-chr --file temp --noweb --make-bed --out SNP
gcta64 --bfile SNP --make-grm --autosome --out tmp
gcta64 --grm tmp --pca --out SNP.pca
## Structure
plink --bfile snp --indep-pairwise 50 10 0.1
plink --bfile snp --extract plink.prune.in --make-bed --out input
for K in 2 3; do admixture -j8 --cv input.bed $K | tee log${K}.out; done
