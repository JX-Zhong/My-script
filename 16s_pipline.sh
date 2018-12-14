#! /bin/bash
#############################################################
# Title: Pipeline of 16S amplicon by Hiseq2500-PE250
# Author: zhongjiaxin
# E-mail: zhongjiaxin1993@163.com
# Date: 3/9/2018
# Enviroment: Ubuntu 16.04x64, qiime 1.9.1, usearch10
# Description: Script for automatic from clean data and mapping file to otu table, taxonomy, alpha and beta raw result
# Requre file list:
# 1. clean data: R1/R2.fastq
# 2. mapping file: map.txt
#############################################################
# Hiseq2500 PE250 of bacterial 16S, data in clean_data/ and each library mapping file in doc/ 
source /share/apps/qiime_software/activate.sh

# 1. 下载实验数据双端测序数据PE250_1/2.fq.gz，和实验设计mappingfile.txt至当前工作目录
mkdir temp result

# # 2. 验证实验设计是否有错误
# validate_mapping_file.py -m mappingfile.txt

# # 3. 合并双端数据
# join_paired_ends.py -f PE250_1.fq -r PE250_2.fq -m fastq-join -o PE250_join

# # 4. 提取barcode
# extract_barcodes.py -f PE250_join/fastqjoin.join.fastq \
# 	-m mappingfile.txt \
# 	-o PE250_barcode \
# 	-c barcode_paired_stitched --bc1_len 0 --bc2_len 6 -a --rev_comp_bc2
# # 5. 拆分样品
# split_libraries_fastq.py -i PE250_barcode/reads.fastq \
# 	-b PE250_barcode/barcodes.fastq \
# 	-m mappingfile.txt \
# 	-o PE250_split/ \
# 	-q 20 --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 --barcode_type 6

# # 6. Remove adaptor
# cutadapt -g AACMGGATTAGATACCCKG -a GGAAGGTGGGGATGACGT -e 0.15 -m 300 --discard-untrimmed PE250_split/seqs.fna -o PE250_P5.fa


# # 7. format to Usearch
# sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' PE250_P5.fa > seqs_usearch.fa

 #依照实验设计批处理并合并
for i in `tail -n+2 map.txt | cut -f 1`;do
  usearch --fastq_mergepairs seq/${i}_1.fq --reverse seq/${i}_2.fq \
  --fastqout temp/${i}.merged.fq --relabel ${i}.
done 

# 合并所有样品至同一文件
cat temp/*.merged.fq > temp/merge.fq
ls -l temp/all.fq
# 质量控制fastq filter, keep reads error rates less than 1%
usearch -fastq_filter merge.fq -fastq_maxee 1.0 \
  -fastaout filtered.fa -relabel Filt
# 8. Dereplication
usearch -fastx_uniques  filtered.fa \
	-fastaout seqs_unique.fa \
	-minuniquesize 2 -sizeout

# 9. Cluster OTU
usearch -cluster_otus seqs_unique.fa -otus otus.fa -uparseout uparse.txt -relabel Otu

# 10. Remove chimeras
usearch -uchime2_ref otus.fa \
	-db /home/jiaxin/USEARCH/rdp_gold.fa \
	-chimeras otus_chimeras.fa \
	-notmatched otus_rdp.fa \
	-uchimeout otus_rdp.uchime \
	-strand plus -mode sensitive -threads 40

grep '>' otus_chimeras.fa|sed 's/>//g' > otus_chimeras.id

filter_fasta.py -f otus.fa -o otus_non_chimera.fa -s otus_chimeras.id -n
grep '>' -c otus_non_chimera.fa 

# 11. 去除非细菌序列
# http://greengenes.secondgenome.com/downloads/database/13_5
## 此处我们用基于reference的去嵌合，下载rdp_gold.fa作
#为reference数据库
#wget http://drive5.com/uchime/rdp_gold.fa

# QIIME提供的下载链接
# 下载Greengene最新数据库，320MB
#wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# 解压数据包后大小3.4G
#tar xvzf gg_13_8_otus.tar.gz
# 将OTU与97%相似聚类的代表性序列多序列比对，大约8min
time align_seqs.py -i otus_non_chimera.fa -t gg_13_8_otus/rep_set_aligned/97_otus.fasta -o aligned/
# 无法比对细菌的数量
grep -c '>' aligned/otus_non_chimera_failures.fasta 
# 获得不像细菌的OTU ID
grep '>' aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > aligned/otus_non_chimera_failures.id
# 过滤非细菌序列
filter_fasta.py -f otus_non_chimera.fa -o otus_rdp_align.fa -s aligned/otus_non_chimera_failures.id -n
# 看我们现在还有多少OTU:975
grep '>' -c otus_rdp_align.fa 


# 12. Generate representitive sequences(rep seqs) and OTU table, remove low abundance samples
# 重命名OTU，这就是最终版的代表性序列，即Reference
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' otus_rdp_align.fa > result/rep_seqs.fa
# 生成OTU表
usearch -usearch_global seqs_usearch.fa -db result/rep_seqs.fa -otutabout otu_table.txt -biomout otu_table.biom -strand plus -id 0.97 -threads 10
# 默认10线程，用时1分20秒，有32.3%的序列匹配到OTU上；用30线程反而用时3分04秒


# 13. 物种注释Taxonomy assignment
assign_taxonomy.py -i result/rep_seqs.fa \
	-r /home/jiaxin/USEARCH/gg_13_8_otus/rep_set/97_otus.fasta \
	-t /home/jiaxin/USEARCH/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt \
	-m rdp -o result

# 14. OTU格式转换，添加信息

biom convert -i otu_table.txt \
	-o result/otu_table.biom \
	--table-type="OTU table" --to-json

biom add-metadata -i result/otu_table.biom \
	--observation-metadata-fp result/rep_seqs_tax_assignments.txt \
	-o result/otu_table_tax.biom \
	--sc-separated taxonomy --observation-header OTUID,taxonomy 

# 15. OTU表数据筛选

# 按样品数据量过滤：选择counts>3000的样品
filter_samples_from_otu_table.py -i result/otu_table_tax.biom -o result/otu_table2.biom -n 3000
# 查看过滤后结果：975个OTU
biom summarize-table -i result/otu_table2.biom

# 按OTU丰度过滤：选择相对丰度均值大于万分之一的OTU
filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i result/otu_table2.biom -o result/otu_table3.biom
# 查看过滤后结果：只有25个样品，346个OTU
biom summarize-table -i result/otu_table3.biom

# # 按物种筛选OTU表：去除p__Chloroflexi菌门
# filter_taxa_from_otu_table.py -i result/otu_table3.biom -o result/otu_table4.biom -n p__Chloroflexi
# # 查看过滤后结果：只有25个样品，307个OTU
# biom summarize-table -i result/otu_table4.biom

# 转换最终biom格式OTU表为文本OTU表格
biom convert -i result/otu_table4.biom -o result/otu_table4.txt --table-type="OTU table" --to-tsv

biom convert -i result/otu_table4.biom -o result/otu_table4.tab --header-key taxonomy --output-metadata-id "ConesensusLinneage" --table-type="OTU table" --to-tsv

# OTU表格式调整方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt
# 筛选最终OTU表中对应的OTU序列
filter_fasta.py -f result/rep_seqs.fa -o tax_rep4.fa -b result/otu_table4.biom


# 16. 进化
clustalo -i tax_rep4.fa -o rep_seqs_align.fa --seqtype=DNA --full --force --threads=30
filter_alignment.py -i rep_seqs_align.fa -o   # rep_seqs_align_pfiltered.fa, only very short conserved region saved
make_phylogeny.py -i rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree # generate tree by FastTree


#Beta/alpha diversity
echo 'alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species,goods_coverage,simpson' >>paramater.txt
echo 'beta_diversity:metrics bray_curtis,euclidean,unweighted_unifrac,weighted_unifrac' >>paramater.txt

core_diversity_analyses.py -i result/otu_table4.biom -o core/ -m map.txt -t  result/rep_seqs.tree -e 40000 -p paramater.txt


normalize_table.py -i result/otu_table4.biom -o temp/otu_table_css.biom -a CSS
# 将mappingfile转换为R可读的实验设计
sed 's/#//' mappingfile.txt > result/design.txt
# 转换文本otu_table格式为R可读
sed '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt > result/otu_table.txt
# 转换物种注释信息为制表符分隔，方便R读取
sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt
