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


ref=cds_matching_uces.fasta
outdir=variants
statdir=stats
mkdir -p $outdir

bam_list=$1
#name=${bam_list#../sample_list_}
#name=${name%.txt}
name1=${bam_list#bamlist_}
name=${name1%.txt}
#command to output all variant sites (usually use this one)
#bcftools mpileup -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD -I -f $ref -b dupmarked_bamlist_all.txt -Ou | bcftools call -m -Ov -o EPOW2202_uce_vars_all_250525.vcf

#Command to output ALL sites
#bcftools mpileup -a FORMAT/AD,FORMAT/DP,FORMAT/SP -f $ref -b $bam_list -Ou | bcftools call -m -Ov -o ${outdir}/${name}_all_sites.vcf

#Command for specific regions only
bcftools mpileup -a FORMAT/AD,FORMAT/DP,FORMAT/SP -f $ref -b $bam_list -T contigs_with_three_changed_sites.txt -Ou | bcftools call -m -Ov -o ${outdir}/${name}_cds.vcf

#For single sample specific regions
#bcftools mpileup -a FORMAT/AD,FORMAT/DP,FORMAT/SP -f $ref NTKG002_AIW_3477_fmrgmd_uces.sorted.bam -T contigs_with_three_changed_sites.txt -Ou | bcftools call -m -Ov -o ${outdir}/${name}_cds.vcf

