export PATH="/Your/Path/To/RepeatMasker":$PATH
export PATH="/Your/Path/To/RepeatModeler":$PATH
export PATH="/Your/Path/To/genometools/":$PATH
export PATH="/Your/Path/To/hmmer/":$PATH
export PATH="/Your/Path/To/blast":$PATH
export PATH="/Your/Path/To/LTR_retriever":$PATH
export PATH="/Your/Path/To/cd-hit":$PATH


genome=
pre=
cpu=

#repeatmodeler
BuildDatabase -name $pre -engine rmblast $genome
RepeatModeler -pa $cpu -database $pre 
dir=`ls |grep RM`
RepeatModeler -pa $cpu -database $pre -LTRStruct -recoverDir ./${dir} 
RepeatMasker -e rmblast -lib ${pre}-families.fa -pa $cpu $genome

#ltr_finder
gt suffixerator -db $genome -indexname $pre -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index $pre -similar 90 -vic 10 -seed 20 -seqids yes  -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1  > ${pre}.harvest.scn

ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 $genome > ${pre}.finder.scn

LTR_retriever -genome $genome -inharvest ${pre}.harvest.scn -infinder ${pre}.finder.scn -threads $cpu

#Merge
cat *LTRlib.fa *families.fa > ${pre}.repeat.lib.fa

RepeatMasker -nolow -no_is -norna -pa $cpu -lib ${pre}.repeat.lib.fa $genome

