#!/bin/bash
# - - - - - - - - - - - - - - - - - - - - -
# @ Jiaxin
# Workflow for processing fastq files from Metagenome sequencing into special table and functiona abundance table
##load conda enviroment py2

######################################
############ Declaration #################
######################################

if [ $# -ne 3 ]; then
	echo "Usage: ./metagenome.sh [path to fastq] [path to work dir] [threads(>1)]"
	exit 65
fi

input=$1
work_dir=$2

mkdir -p ${work_dir}
mkdir -p ${work_dir}/kraken/map
mkdir -p ${work_dir}/kraken/result
mkdir -p ${work_dir}/humann2
mkdir -p ${work_dir}/metaphlan
mkdir -p ${work_dir}/metaphlan2
mkdir -p ${work_dir}/clean_fastq
for i in $(ls ${input} | cut -d "." -f 1 | uniq); do
	fastp -i ${input}/$i".R1.fastq" -I ${input}/$i".R2.fastq" -o ${work_dir}/clean_fastq/ -O ${work_dir}/clean_fastq/

	kraken --preload --db /share/data0/reference/megtagenome_database/minikraken_20141208 \
		--threads $3 ${input}/$i".R1.fastq" ${input}/$i".R2.fastq" >${work_dir}/kraken/$i".kraken.out"
	/share/apps/kraken/kraken-translate \
		--db /share/data0/reference/megtagenome_database/minikraken_20141208 \
		--mpa-format ${work_dir}/kraken/$i".kraken.out" >${work_dir}/kraken/map/$i".kraken.map"
	/share/apps/kraken/kraken-mpa-report \
		--db /share/data0/reference/megtagenome_database/minikraken_20141208 \
		${work_dir}/kraken/$i".kraken.out" >${work_dir}/kraken/$i".kraken.tsv"

	/share/apps/metaphlan2/metaphlan2.py demo.fastq \
		--bowtie2out metaphlan2.bz2 \ 
		-s metaphlan2.sam.bz2 \ 
		--nproc 10 \ 
		--input_type fastq \ 
		--mpa_pkl /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200.pkl \
		--bowtie2db /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200 >metaphlan2.tsv
done

#MetaPhlan
/home/jiaxin/metagnome/metaphlan/metaphlan.py demo.fastq --input_type multifastq --nproc 8 --bowtie2db /share/data0/reference/megtagenome_database/bowtie2db/mpa --bowtie2out decont.bt2.txt -o decont.metaphlan.tsv

/home/jiaxin/metagnome/metaphlan/metaphlan.py decont.bt2.txt --input_type bowtie2out --nproc 2 -t reads_map -o decont.metaphlan.map

#Kraken
kraken --preload --db /share/data0/reference/megtagenome_database/minikraken_20141208 --threads 10 demo.fastq >kraken/kraken.out
/share/apps/kraken/kraken-translate --db /share/data0/reference/megtagenome_database/minikraken_20141208 --mpa-format kraken/kraken.out >kraken/kraken.map
/share/apps/kraken/kraken-mpa-report --db /share/data0/reference/megtagenome_database/minikraken_20141208 kraken/kraken.out >kraken/kraken.tsv

#humann2
humann2 -i demo.fastq -o humann/ --input-format fastq --remove-temp-output --nucleotide-database /share/data0/reference/megtagenome_database/humann2/chocophlan/chocophlan --protein-database /share/data0/reference/megtagenome_database/humann2/uniref --threads 8 --metaphlan-options "--mpa_pkl /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200.pkl  --bowtie2db /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200"
#Use the HUMAnN2 tool renorm table, to compute the normalized abundances
humann2_renorm_table -i humann/demo_genefamilies.tsv -o humann/humann2_genefamilies.relab.tsv --units relab
humann2_renorm_table -i humann/demo_pathabundance.tsv -o humann/humann2_pathabundance.relab.tsv --units relab
