#!/bin/bash
#Job name:mapdamage_all
#SBATCH --job-name=dmg_profiler
#Account:
#SBATCH --account=fc_tsutsuifca
#Partition
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --mail-user=rachel.weinberg@berkeley.edu
#SBATCH --time=6:00:00


module load anaconda3/2024.02-1-11.4
module load gcc/10.5.0 
module load r/4.4.0
source activate damage-profiler #environment located in environments/damage-profiler.yml

batch=$1 #Samples to run DamageProfiler with
ref=$2 #Reference genome (MP2303 contigs)

#Must have sample_list_${batch}.txt somewhere in the working directory

mkdir -p stats/DamageProfiler_${batch}
while read -r line;
do
name=${line%_fmrgmd_uces.sorted.bam}
damageprofiler -i $line -r $ref -t 75 -yaxis_dp_max 0.1 -o stats/DamageProfiler_${batch}/${name};
done <sample_list_${batch}.txt


