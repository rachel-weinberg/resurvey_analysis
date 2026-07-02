module load anaconda3
module load bio/bcftools
module load bio/samtools
source activate seq_analysis


#historical vcf
vcf1=$1

#contemporary vcf for comparison
vcf2=$2 

pop=$3

out=$4


#Get list of transitions from contemporary file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'TYPE="snp" && QUAL>=30 && INFO/DP>=5 && ((REF="A" && ALT="G") || (REF="G" && ALT="A") || (REF="C" && ALT="T") || (REF="T" && ALT="C"))' $vcf2 > ref_transitions_${pop}.txt

#Get list of transitions in historical file and compare to contemporary list w grep
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'TYPE="snp" && ((REF="A" && ALT="G") || (REF="G" && ALT="A") || (REF="C" && ALT="T") || (REF="T" && ALT="C"))' $vcf1 | \
grep -v -F -f ref_transitions_${pop}.txt | \
cut -f1,2 > mask_positions_${pop}.txt

bcftools view --targets-file ^mask_positions_${pop}.txt $vcf1 -o $out


#all sites

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'TYPE="snp" && ((REF="A" && ALT="G") || (REF="G" && ALT="A") || (REF="C" && ALT="T") || (REF="T" && ALT="C"))' $vcf1 | \
grep -v -F -f ref_transitions_${pop}.txt | \
cut -f1,2 > mask_positions_${pop}.txt

bcftools view --targets-file ^mask_positions_${pop}.txt $vcf1 -o $out