export PATH="/Your/Path/To/hifiasm":$PATH

fq= #hifi reads
cpu=

hifiasm -o hifiasm.asm -t $cpu $fq

awk '/^S/{print ">"$2;print $3}' hifiasm.asm.bp.p_ctg.gfa > hifiasm.asm.bp.p_ctg.fasta

