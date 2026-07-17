module load anaconda3
source activate phyluce-1.7.2

#create config file for complete data matrix
# phyluce_assembly_get_match_counts \
#     --locus-db uce-match-$group/probe.matches.sqlite \
#     --taxon-list-config taxon-set-$group.conf \
#     --taxon-group 'all' \
#     --output taxon-set-$group/all/all-taxa-complete.conf


#extract fasta data for all loci in complete and incomplete data matrices
#go to taxon-set folder
cd taxon-set-${group}/all
# #make directory for log files
mkdir -p log

# #extract the FASTA data for this group

# #extract incomplete matrix data
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../${contig_path} \
    --locus-db ../../uce-match-${group}/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output ${group}-incomplete.fasta \
    --incomplete-matrix ${group}-incomplete.incomplete \
    --log-path log


#extract complete matrix data
# phyluce_assembly_get_fastas_from_match_counts \
#     --contigs ../../spades-assemblies-$group/contigs \
#     --locus-db ../../uce-match-$group/probe.matches.sqlite \
#     --match-count-output all-taxa-complete.conf \
#     --output ${group}-complete.fasta \
#     --log-path log



# explode the monolithic FASTA by taxon, you can also do by locus
#incomplete matrix
phyluce_assembly_explode_get_fastas_file \
    --input ${group}-incomplete.fasta \
    --output exploded-fastas-${group}-incomplete \
    --by-taxon

# #complete matrix
# # phyluce_assembly_explode_get_fastas_file \
# #     --input ${group}-complete.fasta \
# #     --output exploded-fastas-${group}-complete \
# #     --by-taxon

# # get summary stats on the FASTAS and save as .csv
# #incomplete matrix
for i in exploded-fastas-$group-incomplete/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv | tee -a phyluce_assembly_fasta_lengths_${group}_incomplete.csv;
done

#complete matrix
# for i in exploded-fastas-$group-complete/*.fasta;
# do
#     phyluce_assembly_get_fasta_lengths --input $i --csv | tee -a phyluce_assembly_fasta_lengths_$group_complete.csv;
# done

#get number of taxa for completeness assessment
# ntaxa=$(echo $(wc taxon-set-$group.conf -l) | awk '{print$1}')

# phyluce_align_seqcap_align \
    # --input ${group}-incomplete.fasta \
    # --output mafft-nexus-edge-trimmed-${group} \
    # --taxa $ntaxa \
    # --aligner mafft \
    # --cores 24 \
    # --incomplete-matrix \
    # --log-path log


# phyluce_align_get_only_loci_with_min_taxa \
    # --alignments  mafft-nexus-edge-trimmed-${group} \
    # --taxa $ntaxa \
    # --percent 0.98 \
    # --output mafft-nexus-edge-trimmed-$group-98p \
    # --cores 24 \
    # --log-path log