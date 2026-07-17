#!/bin/bash
conda activate seq_analysis


vcf=$1
poplist=$2
outfolder=popgen_results
samplelist_dir=sample_lists/resurvey_all

mkdir -p $samplelist_dir

while read -r line;
do
cat sample_list.txt | grep $line > ${samplelist_dir}/${line}_samples.txt
done < $poplist

#Calculate pi over 1 kbp windows
pi_dir=${outfolder}/winpi
mkdir -p $pi_dir
while read -r line;
do
vcftools --vcf $vcf --keep ${samplelist_dir}/${line}_samples.txt \
--window-pi 1000 \
--out ${pi_dir}/${line}
done < $poplist

#Calculate site-pi

while read -r line;
do
vcftools --vcf $vcf \
--keep ${samplelist_dir}/${line}_samples.txt \
--site-pi \
--out ${pi_dir}/${line}
done < $poplist



#Concatenate window pi files into single tsv
while read -r line;
do
sample_name=${line}
awk -v sample_name="$sample_name" '{print sample_name "\t" $0}' ${pi_dir}/${line}.windowed.pi >> win_pi_1kb_noEPOW_250717.tsv
done < $poplist

#Concatenate site pi
while read -r line;
do
sample_name=${line}
awk -v sample_name="$sample_name" '{print sample_name "\t" $0}' ${pi_dir}/${line}.sites.pi >> site_pi_stats_allpops.tsv
done < $poplist


#Calculate Tajima's D with 1kbp windows
tajima_dir=TajimaD_all_samples
mkdir -p $tajima_dir
while read -r line;
do
vcftools --vcf $vcf --keep ${samplelist_dir}/${line}_samples.txt --TajimaD 1000 --out ${tajima_dir}/${line}
done < $poplist
