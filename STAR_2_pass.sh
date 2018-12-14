#!/bin/bash
#
###############################
# Program: 
#	STAR: 2 pass
#	1. Mapping
#	2. count
#   3. merge
# Author: JerAx
# Email: zhongjx1993@gmail.com
######################################
############ Declaration #############
######################################


if [ $# -ne 5 ]
then
	echo "Usage: ./STAR.sh [fq1] [fq2] [inpu_dir] [work_dir] [threads]"
	exit 65 
fi	  
	 
	 
fq1=$1
fq2=$2
inpu_dir=$3
work_dir=$4	
threads=$5
genomeDir="/home/jiaxin/ref/mm_ref_150"
annotation="/home/jiaxin/ref/Mus_musculus.GRCm38.85.gtf"
sample_id=${fq1%%.*} 

################################################################################

#out_qc_dir=${work_dir}/qc/${sample_id}
runDir=${work_dir}/${sample_id}/runDir
clean_fastq=${work_dir}/${sample_id}/clean_fastq
#rsem=${work_dir}/rsem/${sample_id}
source activate py2
#mkdir -p ${out_qc_dir}
mkdir -p ${runDir}
mkdir -p ${clean_fastq}
#mkdir -p ${rsem}

#qc
fastp -i ${inpu_dir}/$fq1 -I ${inpu_dir}/$fq2  -o ${clean_fastq}/$fq1  -O ${clean_fastq}/$fq2


#1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows:

#mkdir $genomeDir
#STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /share/data0/reference/Genome/hg19/hg19.chr.fa --sjdbGTFfile /home/liukuai/single_cell/homo_annotation_hg19/gencode.v19.annotation.gtf --runThreadN 5 --sjdbOverhang # ReadLength-1


#2) Per-sample 2-pass mapping
STAR \
--genomeDir ${genomeDir} --sjdbGTFfile $annotation \
--readFilesIn ${clean_fastq}/$fq1 ${clean_fastq}/$fq2 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic \
--runThreadN $threads --outFileNamePrefix ${runDir}/${sample_id}. \
--sjdbOverhang 149 \
--readFilesCommand zcat

featureCounts -T $threads -p -t exon -g gene_id -a $annotation -o ${runDir}/../${sample_id}".featureCounts.txt" ${runDir}/${sample_id}".Aligned.sortedByCoord.out.bam"
cut -f 1,7 ${runDir}/../${sample_id}".featureCounts.txt"|grep -v '^#'>${runDir}/../${sample_id}".featureCounts.filter.tsv"

rm -r ${runDir}/*_STARtmp
rm -r ${runDir}/*_STARpass1
rm ${runDir}/*.toTranscriptome.out.bam

#merge_metaphlan_tables.py ${work_dir}/${sample_id}/*.featureCounts.filter.tsv >merge.counts.tsv
