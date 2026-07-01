module load anaconda3
source activate fastp

#input txt list of contemporary reads
read_list=$1

while read -r line;
do
R1=$line
R2=${R1/READ1/READ2}
name=${R1%READ1.fastq.gz}
fastp -i $R1 -I $R2 -o ${name}_READ1.fastp2.fastq.gz -O ${name}_READ2.fastp2.fastq.gz --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 40
done<${read_list}