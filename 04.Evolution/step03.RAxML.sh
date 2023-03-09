export PPATH="/Your/Path/To/RAxML":$PATH

pre=
cpu=

raxmlHPC-PTHREADS-AVX -f a -x 12345 -s Single2SuperGene.phylip -N 1000 -p 12345 -m GTRGAMMA -n $pre -T $cpu

