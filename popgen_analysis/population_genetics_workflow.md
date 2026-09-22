This workflow uses the following tools:
1. vcftools 0.1.17
2. plink 1.9
3. admixture 1.3.0

All R packages are loaded in the corresponding R scripts

1. Calculate Tajima's D, Fst, heterozygosity, and nucleotide diversity using vcftools. Requires sample_list.txt (all samples) and pop_list.txt (sample IDs for each pop) be present in the analysis folder

```
vcftools_analyses_resurvey_uces.sh <vcf_file> <pop_list.txt>
```

2. The bulk of the anlaysis is performed in the following script, which contains code to import and analyze vcftools outputs, calculate PCA, allelc richness, isolation by distance, allele sharing distances, heterozygosity, and heirarchical variance components (hierfstat)

```
Resurvey_UCE_popgen_analysis.R
```

3. Admixture cross-validation and unsupervised analysis
```
convert_to_plink_run_admixture.sh <vcf_file> 
```

4. Code to plot CV error for k=1-5 and admixture results for k=3 using 2 individuals per population

```
admixture_plots.R
```

5. Additional code to generate maps and plots of behavioral data b

```
behavior_maps_and_plots.R
```

