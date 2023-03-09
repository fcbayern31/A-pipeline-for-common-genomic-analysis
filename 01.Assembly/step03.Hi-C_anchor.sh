export PATH="/Your/Path/To/juicer/":$PATH
export PATH="/Your/Path/To/3D-DNA/":$PATH
export PATH="/Your/Path/To/samtools/":$PATH

#Put your contig.fa in ref/
#Put your Hi-C.reads_1/2.fq in fastq/

genome=
pre=
enz=
cpu=

python generate_site_positions.py $enz $pre $genome
samtools faidx $genome
cut -f1,2 ${genome}.fai > ${pre}.chrom.sizes

juicer.sh -g $pre -s $enz -z ./ref/$genome -y ./ref/${pre}_${enz}.txt -p ./ref/${pre}.chrom.sizes -t ${cpu} -D /Your/Path/To/juicer

/bin/bash run-asm-pipeline.sh $genome aligned/merged_nodups.txt

#Anchor your contig order and strand mannually by using JuicerBox!

#Hi-C interactions figure 

assem_file=   #*.assembly 
chr_num= 

grep ">" $assem_file |sed 's/>//' > draw.cprops
grep -v ">" $assem_file |head -n $chr_num > draw.asm

/bin/bash /Path/To/Your/3d-DNA/edit/edit-mnd-according-to-new-cprops.sh draw.cprops merged_nodups.txt > new.mnd.txt
/bin/bash /Path/To/Your/3d-DNA/visualize/run-asm-visualizer.sh -p false -c -r 1000000,500000,250000 draw.cprops draw.asm new.mnd.txt
python straw_matrix.py draw.hic  #hic-straw module https://pypi.org/project/hic-straw/
Rscript plot_hic.r draw.matrix

