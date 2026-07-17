pacman::p_load(BiocManager, tidyverse, vcfR, SNPRelate)


vcf_path <- "UCEanalysis/all_pops_ts_filtered.recode.vcf"

vcf <- read.vcfR(
  vcf_path
)

#filter on allele balance
vcf_ab <- filter_allele_balance(vcf, min.ratio = 0.2, max.ratio = 0.8)

#Filter for biallelic and polymorphic sites
vcf_b <- vcf_ab[is.biallelic(vcf), ]
vcf_bp <- vcf_b[is.polymorphic(vcf_b), ]


#Exclude sites with missing data
vcf_ab_nomiss <- missing_by_snp(vcf_bp, cutoff = 1)


vcfR::write.vcf(
  vcf_ab_nomiss,
  file = "UCEanalysis/resurvey_only_ab_100p.vcf.gz"
)
