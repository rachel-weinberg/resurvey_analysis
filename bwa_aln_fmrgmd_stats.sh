#!/bin/bash

# Load required modules
module load bio/samtools
module load bio/bcftools
module load bio/bwa


# Input validation
if [ $# -eq 0 ]; then
    echo "Error: No sample ID provided"
    echo "Usage: $0 <sample_id>"
    exit 1
fi

line=$1
id=$line
ref=/global/scratch/users/rachelweinberg/UCE_Data/alignment/cds_alignments/cds_matching_uces.fasta
r1=../../resurvey_clean/${id}/split-adapter-quality-trimmed/${id}-_READ1.fastp2.fastq.gz
r2=${r1/READ1/READ2}
outdir=.

echo "Processing started for sample: $id"
# Check if required files exist
if [ ! -f "$ref" ]; then
    echo "Error: Reference file $ref not found"
    exit 1
fi

if [ ! -f "$r1" ]; then
    echo "Error: Read 1 file $r1 not found"
    exit 1
fi

if [ ! -f "$r2" ]; then
    echo "Error: Read 2 file $r2 not found"
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "$outdir"
mkdir -p stats

# Check if reference is indexed
if [ ! -f "${ref}.bwt" ]; then
    echo "Indexing reference with BWA..."
    bwa index "$ref"
fi

# Alignment
echo "Aligning reads for sample: $id"
bwa mem "$ref" "$r1" "$r2" > "${outdir}/${id}_pe.sam"

# Convert to BAM
samtools view -bS "${outdir}/${id}_pe.sam" > "${outdir}/${id}_uces.bam"

# Sort BAM
samtools sort "${outdir}/${id}_uces.bam" -o "${outdir}/${id}_uces.sorted.bam"

# Remove unsorted SAN and BAM intermediate files
rm "${outdir}/${id}_uces.bam"
rm "${outdir}/${id}_pe.sam"  

# Change to output directory
cd "$outdir" || exit 1

# Add read groups
samtools addreplacerg -r "@RG\tID:${id}\tPG:samtools addreplacerg\tSM:${id}" --write-index -o "${id}_rg.bam" "${id}_uces.sorted.bam"

# Sort by name for fixmate
samtools sort -n -o "${id}_rg.namesort.bam" "${id}_rg.bam"

# Fixmate command
samtools fixmate -m -r "${id}_rg.namesort.bam" "${id}_fmrg.bam"

# Sort by coordinate for markdup command
samtools sort -o "${id}_fmrg_uces.sorted.bam" "${id}_fmrg.bam"

# Mark duplicates
samtools markdup -s -f "../stats/${id}_markdup_stats" "${id}_fmrg_uces.sorted.bam" "${id}_fmrgmd_uces.sorted.bam"

# Remove intermediate files
rm "${id}_rg.bam" "${id}_rg.namesort.bam" "${id}_fmrg.bam" "${id}_fmrg_uces.sorted.bam"

# Get depth and alignment stats for duplicate marked alignments
samtools depth "${id}_fmrgmd_uces.sorted.bam" > "../stats/${id}_depth"
samtools coverage "${id}_fmrgmd_uces.sorted.bam" > "../stats/${id}_cov"

echo "Processing complete for sample: $id"