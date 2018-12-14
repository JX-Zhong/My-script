#!/bin/bash
set  -e
set -u


#To prepare the reference for mapping you must first index it by typing the following command where <ref.fa> is the path to your reference file:
bwa index <ref.fa>
# it prepares the Burrows Wheeler Transform index for the reference, allowing the aligner to locate where your reads map within that reference.
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' <ref.fa> <read1.fa> <read1.fa> > lane.sam
#Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags:
samtools fixmate -O bam <lane.sam> <lane_fixmate.bam>
#To sort them from name order into coordinate order:
samtools sort -O bam -o <lane_sorted.bam> -T </tmp/lane_temp> <lane_fixmate.sam>

#In order to reduce the number of miscalls of INDELs in your data it is helpful to realign your raw gapped alignment with the Broad’s GATK Realigner.
java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R <ref.fa> -I <lane.bam> -o <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf>
java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R <ref.fa> -I <lane.bam> -targetIntervals <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf> -o <lane_realigned.bam>

#BQSR from the Broad’s GATK allows you to reduce the effects of analysis artefacts produced by your sequencing machines. It does this in two steps, the first analyses your data to detect covariates and the second compensates for those covariates by adjusting quality scores.
java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R <ref.fa> -knownSites >bundle/b38/dbsnp_142.b38.vcf> -I <lane.bam> -o <lane_recal.table>
java -Xmx2g -jar GenomeAnalysisTK.jar -T PrintReads -R <ref.fa> -I <lane.bam> --BSQR <lane_recal.table> -o <lane_recal.bam>


#It is helpful at this point to compile all of the reads from each library together into one BAM, which can be done at the same time as marking PCR and optical duplicates. To identify duplicates we currently recommend the use of either the Picard or biobambam’s mark duplicates tool.
java -Xmx2g -jar MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT INPUT=<lane_1.bam> INPUT=<lane_2.bam> INPUT=<lane_3.bam> OUTPUT=<library.bam>


#Once this is done you can perform another merge step to produce your sample BAM files.
samtools merge <sample.bam> <library1.bam> <library2.bam> <library3.bam>
samtools index <sample.bam>


#If you have the computational time and resources available it is helpful to realign your INDELS again:
java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R <ref.fa> -I <sample.bam> -o <sample.intervals> --known >bundle/b38/Mills1000G.b38.vcf>
java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R <ref.fa> -I <sample.bam> -targetIntervals <sample.intervals> --known >bundle/b38/Mills1000G.b38.vcf> -o <sample_realigned.bam>

#Lastly we index our BAM using samtools:

samtools index <sample_realigned.bam>

#To convert your BAM file into genomic positions we first use mpileup to produce a BCF file that contains all of the locations in the genome. We use this information to call genotypes and reduce our list of sites to those found to be variant by passing this file into bcftools call.

bcftools mpileup -Ou -f <ref.fa> <sample1.bam> <sample2.bam> <sample3.bam> | bcftools call -vmO z -o <study.vcf.gz>

#To prepare our VCF for querying we next index it using tabix:
tabix -p vcf <study.vcf.gz>

#Finally you will probably need to filter your data using commands such as:
bcftools filter -O z -o <study_filtered..vcf.gz> -s LOWQUAL -i'%QUAL>10' <study.vcf.gz>


#examples
bwa index yeast.fasta
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' yeast.fasta y1.fastq y2.fastq > yeast.sam
samtools sort -O bam -T /tmp -l 0 -o yeast.bam yeast.sam
samtools view -T yeast.fasta -C -o yeast.cram yeast.bam
samtools view yeast.cram
samtools mpileup -f yeast.fasta yeast.cram