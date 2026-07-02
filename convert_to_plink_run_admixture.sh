conda activate seq_analysis
#Convert vcf to plink bed file
vcf=$1 #read in vcf containing historical samples. Published analysis uses combined resurvey + historical vcf

bgzip $vcf
tabix -p vcf ${vcf}.gz

#remove multiallelic sites (should already be removed but just in case)
bcftools view -m2 -M2 -v snps ${vcf}.gz -Ov -o uces_filtered_biallelic.vcf
conda deactivate seq_analysis
conda activate plink-1.9

#set bed prefix
output_prefix=uce_historical #substitue with resurvey or historial depending on samples included

plink --allow-extra-chr --vcf uces_filtered_biallelic.vcf --make-bed --out $output_prefix


#navigate to folder with admixture installation


#Cross-validation
for k in 1 2 3 4 5 10 20;
do admixture --cv ${output_prefix}.bed $k | tee log${k}.out; done

for i in log*.out;
do cat $i | grep 'CV error`

#After inspecting CV results, run supervised analysis with k=3

admixture --supervised $output_prefix.bed 3