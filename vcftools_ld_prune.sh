conda activate seq_analysis

vcf=UCEanalysis/resurvey_only_ab_100p.vcf.gz
ld_window=1000
ld_threshold=0.3

#Output LD by site in 1kbp window
vcftools --gzvcf ${vcf} \
    --geno-r2 \
    --ld-window-bp $ld_window \
    --out ld_stats/1kbp_ld_stats 

#Extract sites with r2 > 0.3
awk -v threshold=$ld_threshold '$5 > threshold {print $1"\t"$3}' ld_stats/1kbp_ld_stats.geno.ld | sort | uniq >ld_stats/1kbp_ld_sites_03.txt

#Note: may need to manually remove header row from LD sites file before running next command
sites=ld_stats/1kbp_ld_sites_03.txt
bcftools view -T ^${sites} -o resurvey_ld_pruned.vcf -Ov $vcf
