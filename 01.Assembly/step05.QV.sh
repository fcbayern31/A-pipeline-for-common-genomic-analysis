export MERQURY=/Your/Path/To/merqury
export PATH="/Your/Path/To/merqury/bin":$PATH

fq1=
fq2=
kmer=
pre=
genome=

echo "$fq1\n$fq2" > fq.lst

cat fq.lst | while read id;do
meryl k=$kmer count output $id.meryl $id.fq.gz
done

union_sum.sh $kmer meryl.list $pre
merqury.sh ${pre}.K${kmer}.meryl $genome $pre

