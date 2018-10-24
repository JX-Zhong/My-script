samtools mpileup -BQ0 -d10000 -f ref.fa -q 40 -b filenames.txt | monovar.py -p 0.002 -a 0.2 -t 0.05 -m 50 -f ref.fa -b filenames.txt -o output.vcf
for i in `cat list`; do samtools mpileup -B -d10000 -f ref/hg19.chr.fa -q 40 GSE69405/bam/$i/runDir/$i".Aligned.sortedByCoord.out.bam"|monovar.py -p 0.002 -a 0.2 -t 0.05 -m 50 -f ref.fa -o GSE69405/bam/$i/runDir/$i".vcf"; done
for i in `cat list`; do echo "GSE69405/bam/$i/runDir/"$i".Aligned.sortedByCoord.out.bam">a.b; samtools mpileup -B -d10000 -f ref/hg19.chr.fa -q 40 GSE69405/bam/$i/runDir/$i".Aligned.sortedByCoord.out.bam"|monovar.py -p 0.002 -a 0.2 -t 0.05 -b a.b -m 50 -f ref.fa -o GSE69405/bam/$i/runDir/$i".vcf"; done
for i in `cat list`; do echo "GSE69405/bam/$i/runDir/"$i".Aligned.sortedByCoord.out.bam">a.b; samtools mpileup -B -d10000 -f ref/hg19.chr.fa -q 40 GSE69405/bam/$i/runDir/$i".Aligned.sortedByCoord.out.bam"|monovar.py -p 0.002 -a 0.2 -t 0.05 -b a.b -m 50 -f ref.fa -o GSE69405/bam/$i/runDir/$i".vcf"; done
; done; done; done

#############################Removing Human Contamination#######################################

fastp -i r1.fastq -I r2.fastq -o r1_clean.fastq -O r2_clean.fastq

for fastq in per_sample_fastq/*.fastq
do
  kneaddata \
  --input $fastq \
  --output per_sample_fastq_clean/$(basename "$fastq" .fastq) \
  -db /home/jiaxin/huaman_reference/human_genome/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
  --threads 2 \
  --bypass-trim \
  --output-prefix $(basename "$fastq" .fastq)
done



#MetaPhlan2
/share/apps/metaphlan2/metaphlan2.py demo.fastq --bowtie2out metaphlan2.bz2 -s metaphlan2.sam.bz2 --nproc 10 --input_type fastq --mpa_pkl /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200.pkl --bowtie2db /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200>metaphlan2.tsv


#MetaPhlan
/home/jiaxin/metagnome/metaphlan/metaphlan.py demo.fastq --input_type multifastq --nproc 8   --bowtie2db /share/data0/reference/megtagenome_database/bowtie2db/mpa   --bowtie2out decont.bt2.txt  -o decont.metaphlan.tsv

/home/jiaxin/metagnome/metaphlan/metaphlan.py decont.bt2.txt --input_type bowtie2out   --nproc 2   -t reads_map -o decont.metaphlan.map

#Kraken
kraken --preload --db /share/data0/reference/megtagenome_database/minikraken_20141208 --threads 10 demo.fastq > kraken/kraken.out ; /share/apps/kraken/kraken-translate --db /share/data0/reference/megtagenome_database/minikraken_20141208   --mpa-format kraken/kraken.out > kraken/kraken.map; /share/apps/kraken/kraken-mpa-report --db /share/data0/reference/megtagenome_database/minikraken_20141208   kraken/kraken.out >kraken/kraken.tsv


#humann2
humann2 -i demo.fastq -o humann/ --input-format fastq --remove-temp-output --nucleotide-database /share/data0/reference/megtagenome_database/humann2/chocophlan/chocophlan --protein-database /share/data0/reference/megtagenome_database/humann2/uniref  --threads 8 --metaphlan-options "--mpa_pkl /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200.pkl  --bowtie2db /share/data0/reference/megtagenome_database/db_v20/mpa_v20_m200"
#Use the HUMAnN2 tool renorm table, to compute the normalized abundances
humann2_renorm_table -i humann/demo_genefamilies.tsv -o humann/humann2_genefamilies.relab.tsv --units relab
humann2_renorm_table -i humann/demo_pathabundance.tsv -o humann/humann2_pathabundance.relab.tsv --units relab

#Create a genus level gene families file
humann2_gene_families_genus_level --input $SAMPLE_genefamilies.tsv --output $SAMPLE_genefamilies_genus_level.tsv
#Run HUMAnN2, with the genus level gene families file as input, to get genus level pathways output files
humann2 --input $SAMPLE_genefamilies_genus_level.tsv --output humann2_genus_level_output
#Total genefamily
cut -f1 demo_fastq/demo_genefamilies.tsv | tail -n +3 | grep -v "|" | sort -u | wc -l

 #regroup our CPM-normalized gene family abundance values to enzyme commission (EC) categories
humann2_regroup_table --input demo_fastq/demo_genefamilies-cpm.tsv --output demo_fastq/level4ec-cpm.tsv --groups uniref90_level4ec

#each input file into a single output directory
for f in *.fasta; do humann2 -i $f -o 763577454; done
humann2_join_tables -i 763577454 -o 763577454_genefamilies.tsv --file_name genefamilies
humann2_renorm_table -i 763577454_genefamilies.tsv -o 763577454_genefamilies_cpm.tsv --units cpm



