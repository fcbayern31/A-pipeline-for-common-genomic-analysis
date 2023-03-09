export PATH="/Your/Path/To/BUSCO":$PATH

seq=    #Your genome.fa
mode=genome
out=busco
lineage=/Your/Path/To/embryophyta_odb10/  #Select a lineage of suitable species
cpu=

busco -i $seq -o $out -m $mode -l $lineage -c $cpu

