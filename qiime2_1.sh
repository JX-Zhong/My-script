#Import data
qiime tools import --input-path clean_otu_table.biom --type 'FeatureTable[Frequency]' --source-format BIOMV210Format --output-path qiime2/table.qza
qiime tools import --input-path rep_seqs.fa --type 'FeatureData[Sequence]' --output-path qiime2/sequences.qza
qiime tools import --input-path rep_seqs.tree --type 'Phylogeny[Unrooted]' --output-path qiime2/unrooted-tree.qza


qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   \
--input-path manifest.csv  \
--output-path paired-end-demux.qza   \
--input-format PairedEndFastqManifestPhred33


qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

#注意：需要20bp以上的overlap才能使用dada2拼接，否则会报错。
qiime dada2 denoise-paired  \
 --i-demultiplexed-seqs paired-end-demux.qza \
 --o-table fastq/table.qza \
 --o-representative-sequences fastq/rep-seqs.qza \
 --p-trim-left-f 0 \
 --p-trim-left-r 0 \
 --p-trunc-len-f 180 \
 --p-trunc-len-r 180 \
 --o-denoising-stats denoising-stats-dada2.qza \
 --verbose  \
 --p-n-threads


qiime metadata tabulate \
  --m-input-file denoising-stats-dada2.qza \
  --o-visualization denoising-stats.qzv




# qiime dada2 denoise-single --i-demultiplexed-seqs single-end-demux.qza --p-trunc-len 380 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats stats-dada2.qza --verbose --p-n-threads 48

qiime feature-table summarize \
  --i-table fastq/table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.txt

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv


#建树用于多样性分析
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


#Export data
qiime tools export \
  feature-table.qza \
  --output-dir exported-feature-table

qiime tools export \
  unrooted-tree.qza \
  --output-dir exported-tree


qiime diversity core-metrics-phylogenetic  \
--i-phylogeny rooted-tree.qza  \
--i-table table.qza  \
--p-sampling-depth 35000 \
--m-metadata-file sample-metadata.tsv \
--output-dir core-metrics-results

qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/faith_pd_vector.qza   --m-metadata-file sample-metadata.tsv   --o-visualization core-metrics-results/faith-pd-group-significance.qzv
#Alpha diversity plot 
qiime diversity alpha-group-significance  \
                                       --i-alpha-diversity core-metrics-results/evenness_vector.qza \
                                         --m-metadata-file sample-metadata.tsv \
                                         --o-visualization core-metrics-results/evenness-group-significance.qzv



#rarefaction 
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 55000 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
#beta diversity plot

qiime diversity beta-group-significance   \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv  \
  --m-metadata-column Description  \
  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv  \
 --p-pairwise

qiime emperor plot   --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza   --m-metadata-file sample-metadata.tsv  --o-visualization core-metrics-results/unweighted-unifrac-emperor-DaysSinceExperimentStart.qzv



#训练分类器,使用rep_set文件中的99_otus.fasta数据和taxonomy中的99_OTU_taxonomy.txt数据

#Note: the classifier has already been trained for you ,the V4 region, bound by the 515F/806R primer pair 
qiime feature-classifier classify-sklearn \
 --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


# 导入参考序列
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 99_otus.fasta \
  --output-path 99_otus.qza
# 导入物种分类信息
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza

# Extract reference reads
qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \      #341F引物
  --p-r-primer GACTACHVGGGTATCTAATCC \  #805R引物
  --o-reads ref-seqs.qza
# Train the classifier（训练分类器）
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza


qiime metadata tabulate   --m-input-file taxonomy.qza   --o-visualization taxonomy.qzv


qiime taxa barplot   --i-table table.qza  \
 --i-taxonomy taxonomy.qza   \
 --m-metadata-file sample-metadata.tsv  \
 --o-visualization taxa-bar-plots.qzv

##Differential abundance testing with ANCOM

qiime composition add-pseudocount   --i-table table.qza --o-composition-table comp-gut-table.qza

qiime composition ancom   --i-table comp-gut-table.qza   --m-metadata-file sample-metadata.tsv --m-metadata-column Description --o-visualization ancom-Subject.qzv

#eg:level4 difference
qiime taxa collapse   --i-table table.qza --i-taxonomy taxonomy.qza  --p-level 4  --o-collapsed-table gut-table-l4.qza
qiime composition add-pseudocount   --i-table gut-table-l4.qza   --o-composition-table comp-gut-table-l4.qza
qiime composition ancom   --i-table comp-gut-table-l4.qza   --m-metadata-file sample-metadata.tsv  --m-metadata-column Description   --o-visualization l4-ancom-Subject.qzv



#Practicum: gneiss Analysis ,ilr and Balance Trees ,Like ANCOM, gneiss uses logs, so requires a pseudocount 
qiime gneiss add-pseudocount \
                --i-table table.qza \
                --p-pseudocount 1 \
                --o-composition-table composition.qza


qiime gneiss correlation-clustering \
               --i-table composition.qza \
                --o-clustering hierarchy.qza  


qiime gneiss dendrogram-heatmap \
               --i-table composition.qza \
                --i-tree hierarchy.qza \
               --m-metadata-file sample-metadata.tsv \
                --m-metadata-column BodySite \
                --p-color-map seismic \
               --o-visualization tree_heatmap_by_bodysite.qzv

#Calculate the ilrtransforms 
qiime gneiss ilr-transform \
         --i-table composition.qza \
          --i-tree hierarchy.qza \
           --o-balances balances.qza


#Fit the linear regression model ,As this is a time-course, could instead uselme-regressiongrouped by Subject 
qiime gneiss ols-regression \
         --p-formula "Subject+BodySite+DaysSinceExperimentStart" \
          --i-table balances.qza \
           --i-tree hierarchy.qza \
            --m-metadata-file sample-metadata.tsv \
             --o-visualization regression_summary.qzv


#Well, great, but we seem to have strayed a long way from actual taxa
qiime gneiss balance-taxonomy \
         --i-table composition.qza \
          --i-tree hierarchy.qza \
           --i-taxonomy taxonomy.qza \
            --p-taxa-level 2 \
             --p-balance-name 'y0' \
              --m-metadata-file sample-metadata.tsv \
               --m-metadata-column BodySite \
                --o-visualization y0_taxa_summary.qzv




####q2-picrust2 Running,The required inputs are --i-table and --i-tree  :::FeatureTable[Frequency] and Phylogeny[Rooted],The Feature Table needs to contain the abundances of ASVs (i.e. a BIOM table) and the tree file needs to contain the ASVs placed into the PICRUSt2 reference tree.


wget http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/picrust2_tutorial_files/reference.fna.qza

wget http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/picrust2_tutorial_files/reference.tre.qza




qiime fragment-insertion sepp --i-representative-sequences rep_seqs.qza \
                                                  --p-threads 20 --i-reference-alignment reference.fna.qza \
                                                  --i-reference-phylogeny reference.tre.qza \
                                                  --output-dir tutorial_placed_out


qiime picrust2 custom-tree-pipeline --i-table mammal_biom.qza \
                                    --i-tree tutorial_placed_out/tree.qza \
                                    --output-dir q2-picrust2_output \
                                    --p-threads 46 --p-hsp-method mp \
                                    --p-max-nsti 2

qiime feature-table summarize --i-table q2-picrust2_output/pathway_abundance.qza --o-visualization q2-picrust2_output/pathway_abundance.qzv




qiime tools export --input-path taxonomy.qza --output-path taxtable/


biom convert -i  feature-table.biom  -o table-with-function.tsv --to-tsv
#Change the first line of biom-taxonomy.tsv (i.e. the header) to this:    #OTUID	taxonomy	confidence
biom add-metadata -i  feature-table.biom -o table-with-taxonomy.biom  --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy

#Annotating file named metacyc_pathways_info_prokaryotes.txt from humann2  from 2018.12.2




#filtering
qiime feature-table filter-samples --i-table table.qza --m-metadata-file  ../sample-metadata.tsv  --p-where   "Description IN ('WT','KO')" --o-filtered-table filter_table.qza
qiime feature-table summarize --i-table filter_table.qza --o-visualization filter_table.qzv --m-sample-metadata-file filter_sample-metadata.tsv


#	How	to	export	a	feature	(OTU)	table	and	convert	from	biom	to	.tsv	
(use	following	codes	in	QIIME2)
#	Step	1,	export	OTU	table	
qiime	tools	export	\
	--input-path table.qza	\
		--output-path	phyloseq
#	OTU	table	exports	as	feature-table.biom	so	convert	to	.tsv
biom	convert	-i	phyloseq/feature-table.biom	-o	phyloseq/otu_table.txt	--to-tsv
#	now	you	have	otu_table.txt
#	open	it	up	in	text	edit	and	change	#OTUID	to	OTUID
#	Step	2,	export	taxonomy	table
qiime	tools	export	\
--input-path taxonomy.qza\
		--output-path	phyloseq
#	now	you	have	taxonomy.tsv
#	open	it	up	in	text	edit	and	change	Feature	ID	to	OTUID ,both in taxonomy and otu table
#	Step	3,	export	tree	
qiime	tools	export	\
	--input-path	unrooted-tree.qza	\
		--output-path	phyloseq
#	Step	4,	if	you	filtered	out	any	sequences	(chloroplasts,	mitochondria,	etc)	then	your	taxonomy	and	OTU	tables	are	different	lengths.	QIIME2	doesn’t	filter	out	taxonomy,	so	you	have	to	merge	the	two	files	in	R	and	output	a	merged	file.	(use	following	codes	in	R)
#	set	working	directory	
#
#the first line of taxonomy.txt is "   #OTUID  taxonomy        confidence   " !!!!!!!Please pay attention to case sensitive
#biom add-metadata -i table.biom -o table.w_omd.biom --observation-metadata-fp taxonomy.txt
#convert biom to tsv with taxonomy message
#biom convert -i table.biom -o table.from_biom_w_consensuslineage.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"