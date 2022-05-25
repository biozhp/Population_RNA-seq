## RNA-seq variant calling and quality control
## Author:Zhao Peng (pengzhao@nwafu.edu.cn)
## Northwest A&F University
## This script briefly summarizes the codes in the analysis process, for reference only.

## Make fasta index
java -Xmx10g -jar -Djava.io.tmpdir=temp_dir picard.jar CreateSequenceDictionary R=reference.fa O=reference.dict
samtools faidx reference.fa
## STAR make index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir index --genomeFastaFiles reference.fa --sjdbGTFfile gene.gtf --sjdbOverhang 149
## STAR mapping
STAR --twopassMode Basic --runThreadN 8 --genomeDir index_name --alignIntronMin 1 --alignIntronMax 67475 --outSAMtype BAM SortedByCoordinate --twopass1readsN -1 --sjdbOverhang 149 --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA --outFilterMismatchNmax 4 --outSJfilterReads Unique --outSAMmultNmax 1 --outFileNamePrefix sample_name --outSAMmapqUnique 60 --readFilesCommand gunzip -c --readFilesIn fq1.gz fq2.gz
## BWA make index
bwa index -p index_name seq.fasta
## BWA mapping
bwa aln -t 12 -n 4 index_name fq1.gz > fq1.sai
bwa aln -t 12 -n 4 index_name fq2.gz > fq2.sai
bwa sampe index_name fq1.sai fq2.sai fq1.gz fq2.gz | samtools view -Sbh > bwa.bam
## Select unique mapping reads
samtools view -h STAR.Aligned.sortedByCoord.out.bam | awk -F"\t" '$12=="NH:i:1"||$1~/@SQ/||$1~/@HD/||$1~/@PG/||$1~/@RG/||$1~/@CO/' | samtools view -S -b -o STAR.unique.bam
## Get the position of reads
samtools view STAR.unique.bam | awk -F"\t" {'print $1,$3,$4,$7,$8'} > STAR.txt
## Get the GeneID that reads were aligned (In the bed file, the second columnthe is StartPos-1)
bedtools intersect -a STAR.bed -b gene.bed -wa -wb | bedtools groupby -i - -g 1-4 -c 8 -o collapse > reads.gene.out
## Sort bam (reads in STAR.bwa.filter.bam were aligned to the same region by STAR and BWA at the same time)
java -Xmx20g -jar -Djava.io.tmpdir=temp_dir picard.jar SortSam I=STAR.bwa.filter.bam O=sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
## Make bam index
samtools index -b sort.bam
## Mark Duplicates
java -Xmx20g -jar -Djava.io.tmpdir=temp_dir picard.jar MarkDuplicates I=sort.bam O=markdup.bam METRICS_FILE=markdup.metrics REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2500000 MAX_FILE_HANDLES=1000
## Make bam index
samtools index -b markdup.bam
## Split Reads with N in Cigar
gatk --java-options "-Xmx20G" SplitNCigarReads -R reference.fa -I markdup.bam -O split.bam --create-output-bam-index true --TMP_DIR temp_dir
## Call germline SNPs and indels via local re-assembly of haplotypes
gatk --java-options "-Xmx30G" HaplotypeCaller -R reference.fa -I split.bam --dont-use-soft-clipped-bases -O gatk.raw.vcf.gz --TMP_DIR temp_dir --native-pair-hmm-threads 8 -ERC GVCF
## Merge vcf file names
all_gVCFs=""
for line in $(cat sample.txt)
do
all_gVCFs=${all_gVCFs}"-V "${line}".gatk.raw.vcf.gz "
done
## Import vcfs to GenomicsDB
gatk --java-options "-Xmx20G" GenomicsDBImport $all_gVCFs--genomicsdb-workspace-path chr_name --reference reference.fa --reader-threads 6 --intervals chr_name && \
## Joint genotyping
gatk --java-options "-Xmx20G" GenotypeGVCFs -V gendb://chr_name -R reference.fa -new-qual -O chr_name.vcf
## Merge vcf
gatk --java-options "-Xmx20G" GatherVcfsCloud -I chr.vcf.list -O merge.vcf
## Select a subset of SNPs from vcf
gatk --java-options "-Xmx20G" SelectVariants -V merge.vcf -O snp.vcf --select-type-to-include SNP -R reference.fa
## SNP quality control
gatk --java-options "-Xmx20G" VariantFiltration -R reference.fa -V snp.vcf -O snp.temp.vcf --filter-name "lowQUAL" --filter-expression "QUAL < 30.0" --filter-name "lowGQ" --filter-expression "GQ < 20.0" --filter-name "highFS" --filter-expression "FS > 60.0" --filter-name "lowReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "lowQD" --filter-expression "QD < 2.0" --filter-name "lowMQRankSum" --filter-expression "MQRankSum  < -12.5" --filter-name "lowMQ" --filter-expression "MQ < 40.0" --filter-name "highSOR" --filter-expression "SOR > 3.0"
awk '$1~/^#/ || $7~/PASS/' snp.temp.vcf > snp.filter.vcf
