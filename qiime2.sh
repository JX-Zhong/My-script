#Import data
qiime tools import --input-path clean_otu_table.biom --type 'FeatureTable[Frequency]' --source-format BIOMV210Format --output-path qiime2/table.qza
qiime tools import --input-path rep_seqs.fa --type 'FeatureData[Sequence]' --output-path qiime2/sequences.qza
qiime tools import --input-path rep_seqs.tree --type 'Phylogeny[Unrooted]' --output-path qiime2/unrooted-tree.qza


qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path manifest.csv  --output-path paired-end-demux.qza   --source-format PairedEndFastqManifestPhred33

qiime dada2 denoise-paired   --i-demultiplexed-seqs paired-end-demux.qza --o-table fastq/table.qza --o-representative-sequences fastq/rep-seqs.qza --p-trunc-len-f 180 --p-trunc-len-r 180 --o-denoising-stats fastq/stats-dada2.qza --p-n-reads-learn 15 --verbose
qiime dada2 denoise-single --i-demultiplexed-seqs single-end-demux.qza --p-trunc-len 380 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats stats-dada2.qza --verbose --p-n-threads 48

#Export data
qiime tools export \
  feature-table.qza \
  --output-dir exported-feature-table

qiime tools export \
  unrooted-tree.qza \
  --output-dir exported-tree


qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table.qza   --p-sampling-depth 35000 --m-metadata-file map.txt  --output-dir core-metrics-results

qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/faith_pd_vector.qza   --m-metadata-file map.txt   --o-visualization core-metrics-results/faith-pd-group-significance.qzv
#Alpha diversity plot 
qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/evenness_vector.qza   --m-metadata-file map.txt --o-visualization core-metrics-results/evenness-group-significance.qzv


#beta diversity plot

qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza   --m-metadata-file map.txt  --m-metadata-column Description  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv   --p-pairwise

qiime emperor plot   --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza   --m-metadata-file map.txt  --o-visualization core-metrics-results/unweighted-unifrac-emperor-DaysSinceExperimentStart.qzv



qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads sequences.qza --o-classification taxonomy.qza



qiime metadata tabulate   --m-input-file taxonomy.qza   --o-visualization taxonomy.qzv


qiime taxa barplot   --i-table table.qza   --i-taxonomy taxonomy.qza   --m-metadata-file map.txt  --o-visualization taxa-bar-plots.qzv

##Differential abundance testing with ANCOM

qiime composition add-pseudocount   --i-table table.qza --o-composition-table comp-gut-table.qza

qiime composition ancom   --i-table comp-gut-table.qza   --m-metadata-file map.txt --m-metadata-column Description --o-visualization ancom-Subject.qzv

#eg:level4 difference
qiime taxa collapse   --i-table table.qza --i-taxonomy taxonomy.qza  --p-level 4  --o-collapsed-table gut-table-l4.qza
qiime composition add-pseudocount   --i-table gut-table-l4.qza   --o-composition-table comp-gut-table-l4.qza
qiime composition ancom   --i-table comp-gut-table-l4.qza   --m-metadata-file map.txt  --m-metadata-column Description   --o-visualization l4-ancom-Subject.qzv

for i in `ls S*.vcf`