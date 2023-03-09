export PATH="/Your/Path/To/Minimap2/":$PATH
export PATH="/Your/Path/To/Assemblytics":$PATH

ref=
qry=
cpu=
pre=
anchor=
min=
max=

minimap2 -ax asm5 --cs -t $cpu  $ref $qry > Qry2Ref.sam
#sam2delta.py is from [https://github.com/malonge/RagTag]
python sam2delta.py  Qry2Ref.sam 

Assemblytics Qry2Ref.sam.delta $pre $anchor $min $max

#Treat your SV.bed and gene.gff as a and b

bedtools intersect -a SV.bed -b gene.gff -wa -wb > intersect.txt

