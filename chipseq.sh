#!/bin/bash
#
###############################
# Program: 
#	CHIP-seq: 4 pass
#	1. Filtering
#	2. MApping
#  3. Marking duplicates
#  4. Analyzing peaks
#  need:MACS2 -----bowtie2 ------samtools-----bedtools-----deeptools
# Author: JerAx
# Email: zhongjx1993@gmail.com
######################################
############ Declaration #############
######################################






if [ $# -ne 5 ]
then
	echo "Usage: ./chipseq.sh [fq1] [fq2] [inpu_dir] [work_dir] [threads]"
	exit 65 
fi	 
source activate py2
fq1=$1
fq2=$2
inpu_dir=$3
work_dir=$4	
threads=$5
genomeDir="/home/jiaxin/ref/mm10/mm10"
annotation="/home/jiaxin/ref/mm10_refSeq.bed"
sample_id=${fq1%%.*} 



runDir=${work_dir}/${sample_id}/runDir
clean_fastq=${work_dir}/${sample_id}/clean_fastq



source activate py2
#mkdir -p ${out_qc_dir}
mkdir -p ${runDir}/${sample_id}/plot
mkdir -p ${clean_fastq}
START_TIME=`date +%s`
echo $START_TIME

#QC with fastp
fastp -i ${inpu_dir}/$fq1 -I ${inpu_dir}/$fq2  -o ${clean_fastq}/$fq1  -O ${clean_fastq}/$fq2



#Bowtie2 mapping
bowtie2 -p $threads -3 5 --local -x $genomeDir -U ${clean_fastq}/$fq1|samtools sort -O bam -@ $threads -o ${runDir}/${sample_id}/${sample_id}".bam"
samtools index ${runDir}/${sample_id}/${sample_id}".bam"
samtools markdup -r ${runDir}/${sample_id}/${sample_id}".bam" ${runDir}/${sample_id}/${sample_id}".markdu.bam"
samtools index ${runDir}/${sample_id}/${sample_id}".markdu.bam"
samtools flagstat ${runDir}/${sample_id}/${sample_id}".markdu.bam"  > ${runDir}/${sample_id}/${sample_id}".bam.stat"


#Call peaks by MAC2
macs2 callpeak -c Control.bam -t ${runDir}/${sample_id}/${sample_id}".markdu.bam" -q 0.01 -f BAM -g mm -n $sample_id --outdir ${runDir}/${sample_id} 2> ${runDir}/${sample_id}/${sample_id}.macs2.log

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ${runDir}/${sample_id}/${sample_id}"_summits.bed" >${runDir}/${sample_id}/${sample_id}"_summits.homer.tmp"



findMotifsGenome.pl ${runDir}/${sample_id}/${sample_id}"_summits.homer.tmp" mm10  ${runDir}/${sample_id}/Motif  -len 8,10,12 -p $threads
annotatePeaks.pl    ${runDir}/${sample_id}/${sample_id}"_summits.homer.tmp" mm10  -cpu $threads 1>${runDir}/${sample_id}/${sample_id}".homer.peakAnn.xls" 2>${runDir}/${sample_id}/${sample_id}".homer.annLog.txt"

#bam2bed and bam2bw
bamCoverage  -bs 10 ${runDir}/${sample_id}/${sample_id}".markdu.bam" -p $threads -o ${runDir}/${sample_id}/${sample_id}".bw" --normalizeUsing RPKM --smoothLength 300 --extendReads 200

# #不同peaks的交集并集处理
# bedtools intersect -a ap1_peaks.narrowpeaks -b ag_peaks.narrowpeaks -wa |wc -l

computeMatrix reference-point  \
       --referencePoint TSS  \
       -b 5000 -a 5000  \
       -R $annotation -p $threads \
       -S ${runDir}/${sample_id}/${sample_id}".bw" \
       --skipZeros  -o ${runDir}/${sample_id}/${sample_id}".TSS.matrix.gz" \
       --outFileSortedRegions ${runDir}/${sample_id}/${sample_id}".regions1_genes.bed"

plotHeatmap -m ${runDir}/${sample_id}/${sample_id}".TSS.matrix.gz" -out ${runDir}/${sample_id}/plot/${sample_id}".heatmap.png"
plotHeatmap -m ${runDir}/${sample_id}/${sample_id}".TSS.matrix.gz" -out ${runDir}/${sample_id}/plot/${sample_id}".heatmap.pdf" --plotFileFormat pdf  --dpi 720
plotProfile -m ${runDir}/${sample_id}/${sample_id}".TSS.matrix.gz"  -out ${runDir}/${sample_id}/plot/${sample_id}".profile.png"
plotProfile -m ${runDir}/${sample_id}/${sample_id}".TSS.matrix.gz"  -out ${runDir}/${sample_id}/plot/${sample_id}".profile.pdf" --plotFileFormat pdf --perGroup --dpi 720 



# # 统计reads在全基因组范围的情况
# multiBamSummary bins -bs 1000 --bamfiles 02-read-alignment/ap2_chip_rep1_1_sorted.bam 02-read-alignment/ap2_chip_rep1_2_sorted.bam 02-read-alignment/ap2_chip_rep1_3_sorted.bam 02-read-alignment/ap2_chip_rep2_1_sorted.bam 02-read-alignment/ap2_ctrl_rep1_1_sorted.bam 02-read-alignment/ap2_ctrl_rep1_2_sorted.bam 02-read-alignment/ap2_ctrl_rep2_1_sorted.bam --extendReads 130 -out treat_results.npz
# # 散点图
# plotCorrelation -in treat_results.npz -o treat_results.png --corMethod spearman -p scatterplot
# # 热图
# plotCorrelation -in treat_results.npz -o treat_results_heatmap.png --corMethod spearman -p heatmap
# # 主成分分析
# plotPCA -in treat_results.npz  -o pca.png


END_TIME=`date +%s`
EXECUTING_TIME=`expr $END_TIME - $START_TIME`
echo $EXECUTING_TIME