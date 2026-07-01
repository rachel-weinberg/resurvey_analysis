module load anaconda3
source activate phyluce-1.7.2

cfile=$1
outdir=$2

phyluce_assembly_assemblo_spades \
    --conf $cfile \
    --output $outdir \
    --cores 24 \
    --memory 128

