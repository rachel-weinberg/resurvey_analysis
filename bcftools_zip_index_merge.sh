module load anaconda3
module load bio/samtools
module load bio/bcftools

source activate seq_analysis

#index
while read -r line;
do
bgzip $line
bcftools index ${line}.gz
echo ${line}.gz | tee -a bgzip_vcf_list.txt
done <$1

out=$2

bcftools merge -l bgzip_vcf_list.txt -m snps -Ov -o $out