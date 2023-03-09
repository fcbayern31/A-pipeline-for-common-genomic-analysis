export PATH="/Your/Path/To/fastp/":$PATH
export PATH="/Your/Path/To/fastuniq/":$PATH
export PATH="/Your/Path/To/gce/":$PATH

fq1=
fq2=
pre=
cpu=

fastp -q 20 -l 50 -i $fq1 -o ${pre}.clean_1.fq -I $fq2 -O ${pre}.clean_2.fq --adapter_sequence GATCGGAAGAGCACACGTCTGAACTCCAGTCACTCGTTAGTGCATCTCGTATGCCGTCTTCTGCTTG --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGA --thread ${cpu}

echo "${pre}.clean_1.fq \n ${pre}.clean_2.fq" > IN.lst

fastuniq -i IN.lst -o ${pre}.clean.uniq_1.fq -p ${pre}.clean.uniq_2.fq

peak=
size=  #genome size
H=     #gce haplotype


gce -f Stat.colum -c $peak -g $size -H $H > kmer.table 2>  kmer_gce.log

