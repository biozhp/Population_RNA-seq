## SNP imputation and annotation
## Author:Zhao Peng (pengzhao@nwafu.edu.cn)
## Northwest A&F University
## This script briefly summarizes the codes in the analysis process, for reference only.

## SNP imputation
java -jar beagle.18May20.d20.jar gt=snp.vcf out=snp.impute.vcf nthreads=8
## SNP Annotation (https://pcingola.github.io/SnpEff/se_inputoutput/)
java -jar snpEff.jar build -gff3 -v wheat
java -jar snpEff.jar ann wheat snp.vcf > ann.vcf
