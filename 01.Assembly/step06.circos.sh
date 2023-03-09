export PATH="/Your/Path/To/BLAST":$PATH
export PATH="/Your/Path/To/MCScanX":$PATH
export PATH="/Your/Path/To/circos":$PATH

#step1 collinearity (Optional)

ref=  #pep or cds
qry=  #pep or cds
cpu=
pre=
ref_gff=   #gff format
qry_gff=   #gff format

makeblastdb -dbtype prot -in ${ref}  #pep--prot , cds---nucl
blastp -db ${ref} -query ${qry} -evalue 1e-6 -num_threads $cpu -outfmt 6 -out ${pre}.blast

cat ${ref}.gene.gff| grep mRNA  | awk '{split($9,b,";");print $1 "\t" b[1] "\t" $4 "\t" $5}' | sed 's/ID=//g' >> ${pre}.gff
cat ${qry}.gene.gff| grep mRNA  | awk '{split($9,b,";");print $1 "\t" b[1] "\t" $4 "\t" $5}' | sed 's/ID=//g' >> ${pre}.gff

MCScanX ${pre}

grep -v "#" ${pre}.collinearity > tmp1
python transID.py tmp1 tmp2 ${pre}.gff
python add.color.py tmp2 ${pre}.collinearity.final colors.txt

#setp2 config file
bin_size=
fa=     #only one genome as input
python circos.py -fa $fa -b $bin_size -c prepare_config_file.txt

#step3 draw 
circos -conf circos.conf
#You can refer to Circos's official instructions [10.1101/gr.092759.109]
#Python scripts refered to [https://github.com/ponnhide/pyCircos]


