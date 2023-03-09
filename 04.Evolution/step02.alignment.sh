export PATH="/Your/Path/To/Muscle":$PATH

#The perl script, run_msa.pl, was referenced from [https://github.com/BaconKwan/Perl_programme/blob/6095b224422f7430626b97b0db68fe3427b3f515/kaks_analysis/run_msa.pl]

pre=   
all_pep_seq=   #put all pep used in orthofinder in a fasta file
all_cds_seq=   #put all cds, which translated to the pep above, in a fasta file

ln -s OrthoFinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
ln -s OrthoFinder/Results_*/Orthogroups/Orthogroups.txt
grep -f Orthogroups_SingleCopyOrthologues.txt Orthogroups.txt | sed 's/ /\t/g' |sed 's/://' > Single_copy.family.txt

perl run_msa.pl --outpre $pre --outdir ./ Single_copy.family.txt $all_pep_seq $all_cds_seq 

G_lst=     #2 colunm : Species\tGene ID
S_num=     #number of the species used in orthofinder

python msa2SuperGene.py -l $G_lst -i Muscle.cds.msa -n $S_num -o Single2SuperGene.fasta
awk '{tmp=$0;getline;print tmp"    "$0}' Single2SuperGene.fasta |sed 's/>//g' > Single2SuperGene.phylip

G_len=   #length of the phylip
sed -i '1i\ '$S_num'   '$G_len'' Single2SuperGene.phylip

