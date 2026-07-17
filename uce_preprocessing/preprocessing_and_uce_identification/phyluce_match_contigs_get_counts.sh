
module load anaconda3
source activate phyluce-1.7.2


group=$1

probe_path=$2
contig_path=$3

#match contigs to probes
phyluce_assembly_match_contigs_to_probes \
 --contigs $contig_path \
 --probes $probe_path \
 --output uce-match-${group}

# #make data matrix config file by printing taxon names
echo '[all]' | tee taxon-set-${group}.conf
ls $contig_path | sed 's/\(.*\)\.contigs.fasta/\1/' | tee -a taxon-set-$group.conf

# #get match counts to extract UCE loci and create data matrix config file
# #make the taxon set directory
mkdir -p taxon-set-${group}/all

# #create config file for incomplete data matrix
phyluce_assembly_get_match_counts \
    --locus-db uce-match-${group}/probe.matches.sqlite \
    --taxon-list-config taxon-set-${group}.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-set-${group}/all/all-taxa-incomplete.conf
