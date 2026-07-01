import csv
#adapter sequence are illumina TruSeq-style dual-indexed
#csv headers must be ordered as follows: plate	well	i7	i5	oligoname	libraryname
#will usually want to update names section manually to fit the following convention: <site><year><two digit sample number>
with open('<path_to_csv>', 'r') as f:
    #python3 requires this "with" statement instead of just naming and reading in the file
    d_reader = csv.DictReader(f)
    #get fieldnames from DictReader object and store in list
    headers = d_reader.fieldnames
    tag_seqs = list()
    tag_map = list()
    names = list()

    for line in d_reader:
        tag_seqs.append(headers[2]+"-"+line['oligoname']+":"+line['i7'])
        tag_seqs.append(headers[3]+"-"+line['oligoname']+":"+line['i5'])
        tag_map.append(line['libraryname']+":"+headers[2]+"-"+line['oligoname']+","+headers[3]+"-"+line['oligoname'])
        names.append(line['libraryname'])

with open("illumiprocessor_config.conf", "w") as g: #changed .txt to .conf, maybe generate the correct filetype the first time around?
    
g.write("[adapters]\ni7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG\ni5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT\n")

    g.write("[tag sequences]\n")
    for tag in tag_seqs:
        g.write(tag+"\n")
    g.write("\n[tag map]\n")

    for tag in tag_map:
        g.write(tag+"\n")
    g.write("[names]\n")    
    for name in names:
        g.write(name+":"+name+"\n")