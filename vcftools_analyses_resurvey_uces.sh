#!/bin/bash
conda activate seq_analysis


vcf=$1
poplist=$2
outfolder=popgen_results

while read -r line;
do
cat sample_list.txt | grep $line > sample_lists/${line}_samples.txt
done < $poplist

#Calculate pi over 1 kbp windows
out=${outfolder}/winpi
mkdir -p $out
while read -r line;
do
vcftools --vcf $vcf --keep sample_lists/resurvey_all/${line}_samples.txt \
--window-pi 1000 \
--out ${out}/${line}
done < $poplist

#Calculate site-pi

while read -r line;
do
vcftools --vcf $vcf \
--keep sample_lists/resurvey_all/${line}_samples.txt \
--site-pi \
--out ${out}/${line}
done < $poplist



#Concatenate window pi vcfs into single tsv
out=1kbp_pi_all_samples
while read -r line;
do
sample_name=${line}
awk -v sample_name="$sample_name" '{print sample_name "\t" $0}' ${out}/${line}.windowed.pi >> win_pi_1kb_noEPOW_250717.tsv
done <pop_list.txt


while read -r line;
do
sample_name=${line}
awk -v sample_name="$sample_name" '{print sample_name "\t" $0}' ${out}/${line}.sites.pi >> site_pi_stats_allpops.tsv
done <pop_list.txt


#Calculate Tajima's D with 1kbp windows
out=TajimaD_all_samples
mkdir -p $out
while read -r line;
do
vcftools --vcf $in --keep  ../sample_lists/resurvey_all/${line}_samples.txt --TajimaD 1000 --out ${line}
done <../pop_list.txt

#Calculate heterozygosity by pop
out=het
mkdir -p $out
while read -r line;
do
vcftools --vcf $in --keep ../sample_lists/resurvey_all/${line}_samples.txt --het --out ${out}/${line}
done <../pop_list.txt


vcftools --vcf ../${in} --indv WL2201 --indv WL2202 --het --out het/WL22


#cat sample_list.txt | grep $line > ${line}_samples.txt
vcftools --vcf ../${in} --keep sample_list/${line}_samples.txt --site-pi --out het/WL22
# vcftools --vcf $in --keep ${line}_samples.txt --hardy --out ${out}/${line}
# vcftools --vcf $in --keep ${line}_samples.txt --weir-fst-pop --out ${out}/${line}
# vcftools --vcf $in --keep ${line}_samples.txt --het --out ${out}/${line}
#vcftools --vcf $in --keep ${line}_samples.txt --TajimaD --out ${out}/${line}
vcftools --vcf $in --keep ${line}_samples.txt --indv-freq-burden --out ${out}/${line}
done < pop_list.txt