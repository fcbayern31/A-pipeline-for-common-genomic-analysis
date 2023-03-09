import argparse
import re
import os

parser=argparse.ArgumentParser(usage = '\npython msa2SuperGene.py -l [spe-geneID.list] -i [Single_copy_gene.msa] -o [out.fasta] ')
parser.add_argument('-l', '--lst', required = True, help="A list contain speice and gene ID.")
parser.add_argument('-i', '--input', required= True, help="MSA format file from mega/prank/muscle and so on.")
parser.add_argument('-n', '--number', required= True, help="How many species in your msa?")
parser.add_argument('-o', '--out', required= True , help="Output the table file.")


args=parser.parse_args()
outf=open(args.out,"w")
step=int(args.number)+1

dic_lst={}
dic_SuperGene={}
with open(args.lst,"r") as IN_lst:
	for line in  IN_lst :
		spe=line.strip().split("\t")[0]
		ID=line.strip().split("\t")[1]
		dic_lst[ID]=spe
		if spe not in dic_SuperGene.keys() :
			dic_SuperGene[spe]=""

with open (args.input,"r") as msa:
	lines=msa.readlines()
lines_cluster=[lines[i:i+step] for i in range(0,len(lines),step)]

for each in lines_cluster :
	for gene in each[1:] :
		clean_space=" ".join(gene.split())
		Spe=dic_lst[clean_space.split(" ")[0]]
		Seq=clean_space.split(" ")[1]
		dic_SuperGene[Spe]+=Seq

for k,v in dic_SuperGene.items():
	outf.write(">"+k+"\n"+v+"\n")
	




