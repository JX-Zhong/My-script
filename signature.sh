for i in `ls haplotype/SRR*|cut -d"." -f1|uniq|cut -d"/" -f2`; do java -d64 -jar $gatk -T VariantFiltration -R $ref -V haplotype/$i".hp.vcf" -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o filter/$i".filtered.vcf"; done


for i in `ls haplotype/SRR*|cut -d"." -f1|uniq|cut -d"/" -f2`; do java -d64 -jar $gatk -T VariantFiltration -R $ref -V haplotype/$i".hp.vcf" -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -filterName AD -filter "AD < 5.0" -filterName AF -filter "AF < 0.1"  -o filter/$i".filtered.vcf"; done


java -d64 -jar $gatk -T VariantFiltration -R $ref -V haplotype/SRR2049347.hp.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -filterName DP -filter "DP < 6" -o filter/snp.filtered.vcf


head -58 filter/SRR2049347.filtered.vcf
awk -F "\t" '{OFS="\t";if($7=="PASS");print$0}' filter/SRR2049347.filtered.vcf >tt/ss.vcf



for i in `ls filter/SRR2049*|cut -d "." -f1|cut -d "/" -f2 |uniq`; do head -58 filter/$i".filtered.vcf" >tt/$i".dp.filter.vcf"; awk -F "\t" '{OFS="\t";if($7=="PASS")print$0}' filter/$i".filtered.vcf" >>tt/$i".dp.filter.vcf"; done


awk -F "\t" '{OFS="\t";if(length($4)==1 && length($5)==1)print$0}' tt/SRR2049348.dp.filter.vcf 

awk -F "\t" '{OFS="\t";if(length($4)==1 && length($5)==1 && $1!="chrMT" )print$0}' 


awk -F "\t" '$363>99 && $363<1001 {OFS="\t";print}' merge_clean.maf > 100_1000.maf

tail -n +33 SRR*.dp.filter.vcf.maf>>H358.maf
sort -k5,5 -k6nr,6 H358.filter.maf>H358.filter1.maf

awk '!a[$5 $6]++' H358.filter1.maf>abc.maf


awk -F "\t" '{OFS="\t";$16=$17="H358";print}' cbc.maf >H358.maf


sed '1,32d' SRR*.dp.filter.vcf.maf >>a.maf

sort -k5,5 -k6nr,6  a.maf |awk '!a[$5 $6]++' >b.maf

awk -F "\t" '{OFS="\t";$16=$17="LC-MBT-15";print}' b.maf >c.maf
