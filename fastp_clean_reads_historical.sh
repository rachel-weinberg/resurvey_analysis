
module load anaconda3
source activate fastp

#fastp: -c flag enables base correction for overlapping PE reads, cut-tail trims low-quality reads from 3' ends using sliding window of size 4 
#Use more aggressive trimming for the historical reads

# while read -r line;
# do
R1=WL0304/split-adapter-quality-trimmed/WL0304-READ1.fastq.gz
R2=${R1/READ1/READ2}
name=${R1%READ1.fastq.gz}
fastp -i $R1 -I $R2 -o ${name}_READ1.fastp2.fastq.gz -O ${name}_READ2.fastp2.fastq.gz --detect_adapter_for_pe --qualified_quality_phred 15 --length_required 30 --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 15;
# done<resurvey_clean/fastp_read_list_historical.txt 
