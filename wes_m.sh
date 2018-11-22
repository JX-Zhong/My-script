# ###### necessary documentation
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.64.alt &
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi & 
# nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.64.alt &
###### variable
## set sources path
# reference_path=/home/yudi/hg38
# ref_fasta=/home/yudi/hg38/hg38.fa
## set tools path
mypicard=/home/yudi/tool/picard/picard.jar
mygatk=/share/apps/gatk-4.0.3.0/gatk

## set config

export GATK_LOCAL_JAR=/share/apps/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar

############### reference_bwa_index
#bwa index ${reference_path}/Homo_sapiens_assembly38.fasta.gz

for i in `ls //share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/|cat`
do
############### change read1/2 F1/F2
    mkdir ${i}
    cd ${i}
    mkdir output
    mkdir log
    echo current sample ${i}--now-merge-fa-2-bam--`date`
    read1=`ls //share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/${i}|awk 'NR==1{print}'`
    read2=`ls //share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/${i}|awk 'NR==2{print}'`
    java -Xmx40G -XX:ParallelGCThreads=48 \
    -jar ${mypicard} FastqToSam \
    F1=//share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/${i}/${read1} \
    F2=//share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/${i}/${read2} \
    OUTPUT=./output/unmapped.bam \
    READ_GROUP_NAME=${i}\
    SAMPLE_NAME=mo \
    LB=WES \
    PL=illumina 1>./log/fastq2bam.log 2>&1

########################## bwa_mem ##################################
    echo current sample ${i}--now-bwa-mem--`date`
    bwa mem -t 48 -Y \
        -R"@RG\tID:${i}\tSM:mo\tLB:WES\tPL:illumina" \
        /home/yudi/hg38/hg38.fa \
        //share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/${i}/${read1} \
        //share/data3/174_WES_data/174_child-neuro-wes_fastq_0604/wes_m/${i}/${read2} > ./output/tmp.sam 2> ./log/bwa_mem.log
########################## sam2bam #################################
    echo current sample ${i}--now-sam2bam--`date`
    samtools view -@ 48 -b ./output/tmp.sam > ./output/aligned.bam 2> ./log/sam2bam.log
    rm ./output/tmp.sam

########################## Merge bam alignment with an unmaped bam
########################## add TMP_DIR ,or will get an error 'No Space left on device'
    echo current sample ${i}--now-merge-bam--`date`
    mkdir ./tmp
    java -XX:ParallelGCThreads=50 -Xms66666m -jar /home/yudi/tool/picard/picard.jar \
         MergeBamAlignment \
         ALIGNED_BAM=./output/aligned.bam \
         UNMAPPED_BAM=./output/unmapped.bam \
         OUTPUT=./output/aligned_merged.bam \
         REFERENCE_SEQUENCE=/home/yudi/hg38/hg38.fa \
         PAIRED_RUN=true \
         UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
         MAX_INSERTIONS_OR_DELETIONS=-1 \
         TMP_DIR=./tmp 2> ./log/merge.log
########################## MarkDuplicates
    echo current sample ${i}--now-MarkDuplicates--`date`
    java -XX:ParallelGCThreads=50 -Xms66666m -jar /home/yudi/tool/picard/picard.jar \
         MarkDuplicates \
         INPUT=./output/aligned_merged.bam \
         OUTPUT=./output/aligned_merged_dedup.bam \
         METRICS_FILE=./output/dedup.metrics 2> ./log/dedup.log

########################## SortAndFixTags
########################## Sort aggregated+deduped BAM file and fix tags
    echo current dict${i}--now SortBam--`date`
    set -o pipefail
    java -XX:ParallelGCThreads=50 -Xms66666m -jar /home/yudi/tool/picard/picard.jar \
      SortSam \
      INPUT=./output/aligned_merged_dedup.bam \
      OUTPUT=./output/sorted.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false 2> ./log/SortSam.log

    echo current dict${i}--now FixTags--`date`
    java -XX:ParallelGCThreads=50 -Xms66666m -jar /home/yudi/tool/picard/picard.jar \
      SetNmAndUqTags \
      INPUT=./output/sorted.bam \
      OUTPUT=./output/aligned_sort_dedup_fix.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      REFERENCE_SEQUENCE=/home/yudi/hg38/hg38.fa 2> ./log/fix.log

    rm ./output/aligned.bam
    rm ./output/aligned_merged.bam
    rm ./output/aligned_merged_dedup.bam
    rm ./output/sorted.bam
    rm ./output/unmapped.bam

    echo current dict${i}--now BaseRecalibrator--`date`
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    BaseRecalibrator \
    -R ~/hg38/hg38.fa \
    -I output/aligned_sort_dedup_fix.bam\
    --use-original-qualities \
    -O ./output/recal_data.csv \
    --known-sites ~/hg38/var/dbsnp_146.hg38.vcf \
    --known-sites ~/hg38/var/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    -L ~/hg38/seqtsv/exon_hg38.interval_list \
    --known-sites ~/hg38/var/1000G_phase1.snps.high_confidence.hg38.vcf 2> ./log/BaseRecalibrator.log

    # /share/apps/gatk-4.0.3.0/gatk \
    # --java-options '-Xms66666m' \
    # BaseRecalibratorSpark \
    # -R ~/hg38/hg38.2bit \
    # -I ./output/aligned_sort_dedup_fix.bam \
    # --use-original-qualities \
    # -O ./output/recal_data.csv \
    # --known-sites ~/hg38/var/dbsnp_146.hg38.vcf \
    # --known-sites ~/hg38/var/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    # -L ~/hg38/seqtsv/exon_hg38.interval_list \
    # --known-sites ~/hg38/var/1000G_phase1.snps.high_confidence.hg38.vcf \
    # --TMP_DIR ./tmp \
    # --spark-master local[8] 2> ./log/BaseRecalibrator.log

    echo cunrrent dict${i}--now GatherBQSRReports--`date`
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    GatherBQSRReports \
    -I ./output/recal_data.csv \
    -O ./output/GatherBQSRReports.res 2> ./log/GatherBQSR.log
    
    echo current dict${i}--now ApplyBQSR--`date`
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    ApplyBQSR \
    -R ~/hg38/hg38.fa \
    -I ./output/aligned_sort_dedup_fix.bam \
    -O ./output/bqsr.bam \
    -L ~/hg38/seqtsv/exon_hg38.interval_list \
    -bqsr ./output/GatherBQSRReports.res \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities 2> ./log/ApplyBQSR.log
    
    # /share/apps/gatk-4.0.3.0/gatk \
    # --java-options '-Xms66666m' \
    # ApplyBQSRSpark\
    # -R ~/hg38/hg38.2bit \
    # -I ./output/aligned_sort_dedup_fix.bam \
    # -O ./output/bqsr.bam \
    # -L ~/hg38/seqtsv/exon_hg38.interval_list \
    # -bqsr ./output/GatherBQSRReports.res \
    # --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    # --use-original-qualities \
    # --spark-master local[*] 2> ./log/ApplyBQSRSpark.log

    echo current dict${i}--now index bam--`date`
    cd ./output
    samtools index -@ 8 bqsr.bam
    cd ..
    echo data processing end `date`

    if [ -f "./${i}/output/bqsr.bam" ];then
    rm ./${i}/output/aligned_sort_dedup_fix.bai
    rm ./${i}/output/aligned_sort_dedup_fix.bam
    rm ./${i}/output/aligned_sort_dedup_fix.bam.md5
    fi

    echo current dict${i}--now HaplotypeCaller--`date`
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    HaplotypeCallerSpark \
    -L ~/hg38/seqtsv/exon_hg38.interval_list \
    --emit-ref-confidence GVCF \
    -R ~/hg38/hg38.2bit \
    -I ./output/bqsr.bam \
    -O ./output/Haplo.vcf \
    --TMP_DIR ./tmp \
    --spark-master local[*] 2> ./log/haplo.log

    #set -e
    echo current dict${i}--now indexfeaturefile--`date`
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xmx66666m' \
    IndexFeatureFile \
    -F ./output/Haplo.vcf 2> ./log/indexvcf.log

    echo current dict${i}--now GenotypeGVCFs--`date`
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    GenotypeGVCFs \
    -V ./output/Haplo.vcf\
    -R ~/hg38/hg38.fa \
    -L ~/hg38/seqtsv/exon_hg38.interval_list \
    -O ./output/GenotypeGVCFs.output.vcf \
    -D ~/hg38/var/dbsnp_146.hg38.vcf \
    -G StandardAnnotation 2> ./log/gvcf.log

    echo current dict${i}--now sel snp--`date`
    time \
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    SelectVariants \
    -select-type SNP \
    -V ./output/GenotypeGVCFs.output.vcf \
    -O ./output/snp.vcf 2> ./log/sel_snp.log

    echo current dict${i}--now filter snp--`date`
    time \
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    VariantFiltration \
    -V ./output/snp.vcf \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ./output/snp.filter.vcf 2> ./log/filter_snp.log

    echo current dict${i}--now sel InDel--`date`
    time \
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    SelectVariants \
    -select-type INDEL \
    -V ./output/GenotypeGVCFs.output.vcf \
    -O ./output/indel.vcf 2> ./log/sel_indel.log

    echo current dict${i}--now filter InDel--`date`
    time \
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    VariantFiltration \
    -V ./output/indel.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ./output/indel.filter.vcf 2> ./log/filter_indel.log

    echo current dict${i}--now merge vcf--`date`
    time \
    /share/apps/gatk-4.0.3.0/gatk \
    --java-options '-Xms66666m' \
    MergeVcfs \
    -I ./output/snp.filter.vcf \
    -I ./output/indel.filter.vcf \
    -O ./output/filter.vcf 2> ./log/merge.log
    
    echo current dict${i}--now snpEff ann--`date`
    java -Xmx6g -jar ~/tool/snpEff/snpEff/snpEff.jar \
    hg38 ./output/filter.vcf > snpEff_ann.vcf 2> ./log/snpEff_ann.log
    
    cd ..
done

