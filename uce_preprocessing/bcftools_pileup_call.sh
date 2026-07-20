#!/bin/bash
#Job name:
#SBATCH --job-name=call_variants_uce_samples_by_pop
#Account:
#SBATCH --account=fc_tsutsuifca
#Partition
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rachel.weinberg@berkeley.edu
#SBATCH --time=12:00:00

module load anaconda3
module load bio/samtools
module load bio/bcftools
module load bio/bwa
source activate seq_analysis


ref=MP2303_uces.fasta
outdir=variants
statdir=stats
mkdir -p $outdir

bam_list=$1


name1=${bam_list#bamlist_}
name=${name1%.txt}

#Output all variant sites 
bcftools mpileup -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD -I -f $ref -b dupmarked_bamlist_all.txt -Ou | bcftools call -m -Ov -o EPOW2202_uce_vars_all_250525.vcf

