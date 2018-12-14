#!/bin/bash
#
# AUTHOR: JerAx

workdir=$1
id=$2
outputdir=$3
mkdir -p $outputdir



kraken2 $workdir/$id".fastq"  --db /share/data0/reference/megtagenome_database/kraken2_Normal --threads 3 --output $outputdir/$id".kmer"
j=1

awk -F "\t" '$3!=9606&&$3!=0&&$3!=28384' $outputdir/$id".kmer"   >$outputdir/$id".kmer.filter"
for i in `cut -f3 $outputdir/$id".kmer.filter"`;do cut -f5 $outputdir/$id".kmer.filter" |sed -n ''"$j"'p'|grep -o "$i:[0-9]*"|cut -d":" -f2>$outputdir/tmp;for a in `cat $outputdir/tmp`;do if [ $a -gt "8" ]; then sed -n ''"$j"'p' $outputdir/$id".kmer.filter" >>$outputdir/filter.id; fi; done;j=$[j+1];done
uniq $outputdir/filter.id>$outputdir/filter.id
for i in `cut -f2 $outputdir/filter.id`;do grep "$i" -A3  $workdir/$id".fastq">>$workdir/$id".fastq.uniq"
rm $outputdir/tmp
rm $outputdir/filter.id
rm $outputdir/$id".kmer"
kraken2 $workdir/$id".fastq.uniq" --db /share/data0/reference/megtagenome_database/kraken2_Normal --threads 3 --output $outputdir/$id".final.kmer" --report $outputdir/$id".taxonomy" --use-mpa-style
