#You need to make the index for the following databases:
#1.NR
#2.Trembl
#3.Swissprot
#4.Interproscan
#5.KEGG

export PATH="/Your/Path/To/diamond":$PATH

cpu=
lib=
qry=    #Your pep file from gene annotation
out=

diamond blastp --threads ${cpu} --db $lib -q $qry -o $out
