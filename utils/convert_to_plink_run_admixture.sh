conda activate seq_analysis
#Convert vcf to plink bed file
vcf=$1 #read in vcf containing historical and resurvey samples. Published analysis uses combined resurvey + historical vcf

bgzip $vcf
tabix -p vcf ${vcf}.gz

#remove multiallelic sites (should already be removed but just in case)
bcftools view -m2 -M2 -v snps ${vcf}.gz -Ov -o uces_filtered_biallelic.vcf

conda deactivate
conda activate plink-1.9

#set bed prefix
output_prefix=uce_historical #substitue with resurvey or historial depending on samples included

plink --allow-extra-chr --vcf uces_filtered_biallelic.vcf --make-bed --out $output_prefix

#Replace all chromosome codes with "1"
sed 's/^[^\t ]*\([[:space:]]\)/1\1/g' ${output_prefix}.bim > temp.bim
mv temp.bim ${output_prefix}.bim


#navigate to folder with admixture installation


#Cross-validation
for k in 1 2 3 4 5;
do admixture --cv --seed=123 ${output_prefix}.bed $k | tee log${k}.out; done

for i in log*.out;
do cat $i | grep 'CV error'; done