#!/bin/sh
# Add cluster information here (if needed)


# - - - - - - - - - - - - - - - - - - - - -
# @ Jiaxin
# Workflow for processing fastq files from
# 16s rRNA sequencing into an OTU table in
# .biom format. This workflow assumes the
# data has a forward/reverse and index from
# an Illumina Miseq run.
# - - - - - - - - - - - - - - - - - - - - -

# Load QIIME 1.9.1

######################################
############ Declaration #############
######################################
if [ $# -ne 3 ]
then
	echo "Usage: ./16s_pipeline.sh [path to merged:eg./home/peng/ocean_result/merged] [path to work dir] [threads(>1)]"
	exit 65
fi

input=$1
work_dir=$2

# - - - - - - - - - - - - - - - - - - - - -
# Make directories to store the data
# - - - - - - - - - - - - - - - - - - - - -

mkdir -p ${work_dir}
mkdir -p ${work_dir}/report
mkdir -p ${work_dir}/fastqc
#export PATH=/share/apps/microbiomeutil-r20110519/ChimeraSlayer:$PATH
source /share/apps/qiime_software/activate.sh

#clean data with hg19
#for fastq in fasta/*/*.fastq; do kneaddata --input $fastq --output per_sample_clean/$(basename "$fastq" .fastq) -db /home/jiaxin/huaman_reference/human_genome/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens --threads 2 --bypass-trim --output-prefix $(basename "$fastq" .fastq) ; done

# Enable reverse strand matching for uclust so that it double checks both sequence possibilities.
echo "pick_otus:enable_rev_strand_match True" >  ${work_dir}/pick_open_params.txt
# Change to program to RDP for assigning taxnomy of the OTU's instead of uclust 
echo "assign_taxonomy:assignment_method rdp" >> ${work_dir}/pick_open_params.txt
# Change to confidence to 0.5 in the RDP classifer
echo "assign_taxonomy:confidence 0.5" >> ${work_dir}/pick_open_params.txt
#rm ${input}/*/*un2.fastq
#rm ${input}/*/*un1.fastq
#fastqc ${input}/*/*.fastq -t $3 -o ${work_dir}/fastqc
#tar -cvzf ${work_dir}/fastqc.tar.gz ${work_dir}/fastqc
#split library & QC
multiple_split_libraries_fastq.py -i ${input} -o ${work_dir}/sl_out --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

#multiple_split_libraries_fastq.py -i ${input} -o ${work_dir}/sl_out --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name
#multiple_split_libraries_fastq.py -i ${input} -o ${work_dir}/sl_out --demultiplexing_method sampleid_by_file --sampleid_indicator .fastq -p /home/peng/ml_dev3/para.txt

# - - - - - - - - - - - - - - - - - - - - -
# Pick OTU's Workflow
# - - - - - - - - - - - - - - - - - - - - -
pick_open_reference_otus.py -i ${work_dir}/sl_out/seqs.fna -o ${work_dir}/open_otu -p ${work_dir}/pick_open_params.txt -a -O $3
echo "Picking open-reference otu's...completed."

# - - - - - - - - - - - - - - - - - - - - -
# Remove Chimeric Sequences
# You must have ChimeraSlayer + BLAST
# installed!
# - - - - - - - - - - - - - - - - - - - - -
identify_chimeric_seqs.py -m ChimeraSlayer -i ${work_dir}/open_otu/pynast_aligned_seqs/rep_set_aligned.fasta -a /share/apps/qiime_software/chimeraslayer-4.29.2010-release/ChimeraSlayer/../RESOURCES/rRNA16S.gold.NAST_ALIGNED.fasta -o ${work_dir}/open_otu/chimeric_check.txt
echo "Identifying chimeric sequences...completed."

#filter fasta
filter_fasta.py -f ${work_dir}/open_otu/pynast_aligned_seqs/rep_set_aligned.fasta -o ${work_dir}/open_otu/pynast_aligned_seqs/non_chimeric_rep_set_aligned.fasta -s ${work_dir}/open_otu/chimeric_check.txt -n
echo "Filtering fasta sequences...completed."

#filter alignment
filter_alignment.py -o ${work_dir}/open_otu/pynast_aligned_seqs_no_chimeric -i ${work_dir}/open_otu/pynast_aligned_seqs/non_chimeric_rep_set_aligned.fasta
echo "Filtering alignment...completed."

# Build phylogenetic tree command
make_phylogeny.py -i ${work_dir}/open_otu/pynast_aligned_seqs_no_chimeric/non_chimeric_rep_set_aligned_pfiltered.fasta -o ${work_dir}/open_otu/rep_set_no_chimeric.tre
echo "Making phylogenetic tree...completed."

#filter otu table chimeric
filter_otus_from_otu_table.py -i ${work_dir}/open_otu/otu_table_mc2_w_tax_no_pynast_failures.biom -o ${work_dir}/open_otu/otu_table_non_chimeric_final.biom -e ${work_dir}/open_otu/chimeric_check.txt

biom convert -i ${work_dir}/open_otu/otu_table_non_chimeric_final.biom -o ${work_dir}/open_otu/otu_table_non_chimeric_final.biom.consensuslineage.txt --to-tsv --header-key taxonomy --output-metadata-id "Consensus Lineage"
# copy important result to report
cp ${work_dir}/open_otu/otu_table_non_chimeric_final.biom.consensuslineage.txt ${work_dir}/report
cp ${work_dir}/open_otu/otu_table_non_chimeric_final.biom ${work_dir}/report
cp ${work_dir}/open_otu/rep_set_no_chimeric.tre ${work_dir}/report
tar -cvzf ${work_dir}/report.tar.gz ${work_dir}/report
