export PATH="/Your/Path/To/hisat2":$PATH
export PATH="/Your/Path/To/TransDecoder":$PATH
export PATH="/Your/Path/To/stringtie":$PATH
export PATH="/Your/Path/To/samtools":$PATH
export PATH="/Your/Path/To/augustus":$PATH
export PATH="/Your/Path/To/maker":$PATH
export PATH="/Your/Path/To/snap":$PATH
export AUGUSTUS_CONFIG_PATH="/Your/Path/To/augustus/config"

fa=
pre=
cpu=

#step1 RNA-seq
fq1=
fq2=

hisat2-build $fa ${pre}.lib
hisat2 --dta -p ${cpu} -x ${pre}.lib -1 $fq1 -2 $fq2 | samtools sort -@ 32 > hisat2.sorted.bam
stringtie -p ${cpu} -o hisat2.sorted.gtf hisat2.sorted.bam

gtf_genome_to_cdna_fasta.pl hisat2.sorted.gtf $fa > transcripts.fasta
gtf_to_alignment_gff3.pl hisat2.sorted.gtf > transcripts.gff3

TransDecoder.LongOrfs -t transcripts.fasta
TransDecoder.Predict -t transcripts.fasta
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > ${pre}.transdecoder.gff3

#step2 De novo
#augustus
autoAugTrain.pl --genome=$fa --trainingset=${pre}.transdecoder.gff3 --species=${pre}

#SNAP #####running maker before
maker2zff ${pre}.maker.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl snap . > snap.hmm

#step3 maker : vim the opt.ctl file 
#est_gff=${pre}.transdecoder.gff3
#protein=Your homo species protein sequences
#rmlib=repeat.lib.fa
#snaphmm=snap.hmm
#augustus_species=${pre}

mpirun -np ${cpu} maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl --ignore_nfs_tmp 


