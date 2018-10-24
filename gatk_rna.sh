#!/bin/bash
#
# AUTHOR: JerAx
# Single cell RNA-Seq variant calling pieline accoring to GATK Best practices.
# https://www.broadinstitute.org/gatk/guide/article?id=3891
#database can be downloadde  at ftp://ftp.broadinstitute.org/bundle/hg19/
# Call with following arguments
# bash scrna_seq_variant_pipeline.sh <input_dir> <Input filename> <output_dir>
# 
# Assumes STAR aligner is under path
#
# GATK bundle set : one can obtain these from gatk ftp (knonw as gatk bundle)
#  ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
#

workdir=$1
i=$2
outputdir=$3


#Path to gatk and picard tools and ref
ref="/share/data0/reference/Genome/hg19/hg19.chr.fa"
gatk="/home/jiaxin/bio/gatk/GenomeAnalysisTK.jar"
picard="/home/jiaxin/bio/gatk/picard.jar"

#Path to gatk bundle set files
millsIndels="/home/jiaxin/ref/vcf_ref/vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
KGIndels="/home/jiaxin/ref/vcf_ref/vcf/1000G_phase1.indels.hg19.sites.vcf"
dbSNP138="/home/jiaxin/ref/vcf_ref/vcf/dbsnp_138.hg19.vcf"



#Create an output directory
mkdir -p $outputdir


#picard add tittle
echo -e "["$(date)"]\tAdd RGline.."
java -jar $picard AddOrReplaceReadGroups I=$workdir/$i".Aligned.sortedByCoord.out.bam" O=$outputdir/$i".Aligned.sortedAddRG.bam" SO=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$i 2>$outputdir/$i".AddRGline.log"

#picard mark duplicates
echo -e "["$(date)"]\tMarking duplicates.."
java -jar $picard MarkDuplicates I=$outputdir/$i".Aligned.sortedAddRG.bam" O=$outputdir/$i".dupMarked.bam" M=$outputdir/$i".dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$outputdir/$i".MarkDuplicates.log"

#SplitNCigarReads
echo -e "["$(date)"]\tSpliting reads.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T SplitNCigarReads -R $ref -I $outputdir/$i".dupMarked.bam" -o $outputdir/$i".split.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>$outputdir/$i".SplitNCigarReads.log"

rm $outputdir/$i".Aligned.sortedAddRG.bam"
rm $outputdir/$i".dupMarked.bam"
rm $outputdir/$i".dupMarked.bai"

#Create targets for indel realignment
echo -e "["$(date)"]\tCreating targets for indel realignment.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T RealignerTargetCreator -R $ref -I $outputdir/$i".split.bam" -o $outputdir/$i".intervals" -nt 50 -known $millsIndels -known $KGIndels 2>$outputdir/$i".indel.log"

#Perform indel realignment
echo -e "["$(date)"]\tPerforming Indel Realignment.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T IndelRealigner -R $ref -I $outputdir/$i".split.bam" -targetIntervals $outputdir/$i".intervals" -o $outputdir/$i".processed.bam"  -known $millsIndels -known $KGIndels 2>$outputdir/$i".indel2.log" 

rm $outputdir/$i".split.bam"
rm $outputdir/$i".split.bai"

#Perform BQSR
echo -e "["$(date)"]\tPerforming BQSR.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T BaseRecalibrator -I $outputdir/$i".processed.bam" -R $ref -knownSites $KGIndels -knownSites $millsIndels -knownSites $dbSNP138 -o $outputdir/$i".recal.table" 2>$outputdir/$i".BQSR1.log" 

#Print recalibrated reads
echo -e "["$(date)"]\tPrinting recalibrated reads.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T PrintReads -R $ref -I $outputdir/$i".processed.bam" -nct 4 -BQSR $outputdir/$i".recal.table" -o $outputdir/$i".recal.bam" 2>$outputdir/$i".BQSR2.log" 

rm $outputdir/$i".processed.bam"
rm $outputdir/$i".processed.bai"

#Run HaplotypeCaller
echo -e "["$(date)"]\tRunning HaplotypeCaller.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T HaplotypeCaller -R $ref  -I $outputdir/$i".recal.bam" -dontUseSoftClippedBases -stand_call_conf 20.0 -nct 64 -o $outputdir/$i".hp.vcf"  2>$outputdir/$i".HaplotypeCaller.log" 


#Filter variants
echo -e "["$(date)"]\tFiltering Variants.."
java -d64 -Djava.io.tmpdir=/home/jiaxin/tmp -jar $gatk -T VariantFiltration -R $ref -V $outputdir/$i".hp.vcf" -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0"  -filterName DP -filter "DP < 30.0" -o $outputdir/$i".filter.vcf" 2>$outputdir/$i".VariantFilter.log"

rm $outputdir/$i".recal.bam"
rm $outputdir/$i".recal.bai"
head -58 $outputdir/$i".filter.vcf" >$outputdir/$i".filtered.vcf"

awk 'length($4)==1&&length($5)==1&&$7=="PASS"' $outputdir/$i".filter.vcf">>$outputdir/$i".filtered.vcf"
rm $outputdir/$i".hp.vcf"
rm $outputdir/$i".intervals"
rm $outputdir/$i".filter.vcf"
mkdir -p $outputdir/final
mv $outputdir/$i".filtered.vcf" $outputdir/final/$i".filtered.vcf"

echo -e "["$(date)"]\tCalling Variants Successfully"