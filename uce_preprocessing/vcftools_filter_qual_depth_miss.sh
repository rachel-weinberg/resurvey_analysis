module load bio/samtools
module load bio/vcftools
module load bio/bcftools
module load anaconda3
source activate seq_analysis

#Get missingness by individual
vcf=$1
ref=$2
pop=${vcf%.vcf}

#Get missing individuals
vcftools --vcf $vcf --missing-indv 

#Get list of individuals with >10% missing data
awk '$5 > 0.1' out.imiss | cut -f1 > lowDP.indv

vcftools --vcf $vcf --remove lowDP.indv --out ${pop}_90p --recode --recode-INFO-all 




out=${vcf%.vcf}_q30_dp5_mac2
outdir=variants
statdir=stats
mkdir -p $outdir
name=${vcf%.recode.vcf}


vcftools --vcf ${pop}_90p.recode.vcf --minQ 30 --min-meanDP 5 --mac 2 --recode --recode-INFO-all --out $out

#index the output vcf
bgzip ${out}.recode.vcf
tabix -p vcf ${out}.recode.vcf.gz

#Normalize indels
bcftools norm -f $ref -o ${out}.norm.vcf ${out}.recode.vcf.gz