export PATH="/Your/Path/To/MUMmer":$PATH

ref=
qry=
pre=
cpu=

nucmer -t $cpu --prefix $pre $ref $qry
delta-filter -l 10000 ${pre}.delta > ${pre}.delta.best
show-coords -r -T ${pre}.delta.best > ${pre}.delta.best.coords
mummerplot -p ${pre} ${pre}.delta.best -t postscript --color

show-snps -q ${pre}.delta.best > ${pre}.snp.txt


