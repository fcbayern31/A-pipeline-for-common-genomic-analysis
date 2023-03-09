from Bio import SeqIO
import sys
import argparse
import os
import re
import gzip
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-fa',required=True)
parser.add_argument('-b',required=True)
parser.add_argument('-c',required=False)
argv = vars(parser.parse_args())
fa = open(argv["fa"].strip(),'r')
BIN =int(argv["b"].strip())
if argv["c"]:
	conf = open(argv["c"].strip(),'r')
	CON={}
	for line in conf:
		L=line.strip().split("=")
		CON[L[0]]=L[1]

def gc(fa,BIN):
	outname1 =  'first.circos'
	outname2 =  'gc.circos'
	outname3 = 'gene_background_genome.txt'
	output1 = open(outname1, 'w')
	output2 = open(outname2, 'w')
	output3 = open(outname3, 'w')
	CO=[]
	with open("colors.conf","r") as color :
		for line in color:
			if line.startswith("continuous"):
				CO.append(line.strip().split("= ")[1])
	n=0
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm = seq_record.id
		seq=str(seq_record.seq).strip()
		if "unanchor" in chrm: continue
		output1.write("chr\t-\t"+chrm+"\t"+chrm+"\t0\t"+str(len(seq))+"\t"+CO[n]+"\n")
		output3.write(chrm+"\t0\t"+str(len(seq))+"\n")
		n+=1
		if n==50:
			n=0
		if len(seq)%BIN==0:
			B=int(len(seq)/BIN)
		else:
			B=int(len(seq)/BIN)+1
		for i in range(B):
			if (i+1)*BIN<=len(seq):
				SEQ=seq[i*BIN:(i+1)*BIN]
			else:
				SEQ=seq[i*BIN:len(seq)]
			l_SEQ=len(SEQ)
			GC=round(100 * len(re.findall('[GCgc]', SEQ)) / l_SEQ, 2)
			output2.write("{0}\t{1}\t{2}\t{3}\n".format(chrm,i*BIN+1, (i+1)*BIN, GC))
num_list=["gene","tRNA","rRNA","snRNA","miRNA"]
	
def number(fa,BIN,gff,TYPE):
	dic={}
	outname = TYPE+'.circos'
	output = open(outname, 'w')
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm = seq_record.id
		l_seq=len(str(seq_record.seq).strip())
		if "unanchor" in chrm: continue
		if l_seq%BIN==0:
			B=int(l_seq/BIN)
		else:
			B=int(l_seq/BIN)+1
		for i in range(B):
			start=i*BIN+1
			end=(i+1)*BIN
			ID=chrm+"\t"+str(start)+"\t"+str(end)
			dic[ID]=0
	for line in gff:
		if line.startswith("#"): continue
		if "unanchor" in line: continue
		if "Pseudo" in line: continue
		L=line.strip().split("\t")
		if TYPE=="gene" :
			if L[2]=="mRNA":
				for key in dic.keys():
					k=key.split("\t")
					if L[0]==k[0] and int(L[3]) >int(k[1]) and int(L[3]) < int(k[2]) :
						dic[key]+=1
		else:
			output.write("{0}\t{1}\t{2}\n".format(L[0],L[3],L[4]))
			#for key in dic.keys():
				#k=key.split("\t")
				#output.write("{0}\t{1}\t{2}\n".format(k[0],k[1], k[2],dic[key]))
				#if int(L[3]) >int(k[1]) and int(L[3]) < int(k[2]) :
					#dic[key]+=1
	if TYPE=="gene" :	
		for key in dic.keys(): 
			k=key.split("\t")
			output.write("{0}\t{1}\t{2}\t{3}\n".format(k[0],k[1], k[2],dic[key]))
def repeat_old(fa,BIN,gff):
	dic={}
	outname = 'repeat.circos'
	output = open(outname, 'w')
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm = seq_record.id
		l_seq=len(str(seq_record.seq).strip())
		if "unanchor" in chrm: continue
		if l_seq%BIN==0:
			B=int(l_seq/BIN)
		else:
			B=int(l_seq/BIN)+1
		for i in range(B):
			start=i*BIN+1
			end=(i+1)*BIN
			ID=chrm+"\t"+str(start)+"\t"+str(end)
			lst = [0] * BIN
			dic[ID]=lst
	for line in gff:
		if "unanchor" in line: continue
		L=line.strip().split("\t")
		for key in dic.keys():
			k=key.split("\t")
			if  int(L[3])> int(k[1]) and  int(L[3]) < int(k[2]) and int(L[4]) > int(k[1]) and  int(L[4]) < int(k[2]):
				LST=[1]*(int(L[4])-int(L[3])+1)
				dic[key][(int(L[3])%BIN-1):int(L[4])%BIN]=LST
			elif int(L[3])>int(k[1]) and  int(L[3]) < int(k[2]) and int(L[4]) > int(k[2]):
				LST=[1]*(int(k[2])-int(L[3])+1)
				dic[key][(int(L[3])%BIN-1):BIN]=LST
			elif int(L[3]) < int(k[1]) and int(L[4]) > int(k[1]) and  int(L[4]) < int(k[2]) :
				LST=[1]*(int(L[4])%BIN)
				dic[key][0:int(L[4])%BIN]=LST
	for key in dic.keys():
		k=key.split("\t")
		output.write("{0}\t{1}\t{2}\t{3}\n".format(k[0],k[1], k[2],dic[key].count(1)))
def repeat(fa,fa_m,BIN):
	dic={}
	outname = 'repeat.circos'
	output = open(outname, 'w')
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm1 = seq_record.id
		seq1=str(seq_record.seq).strip()
		l_seq=len(str(seq_record.seq).strip())
		if "unanchor" in chrm1: continue
		if l_seq%BIN==0:
			B=int(l_seq/BIN)
		else:
			B=int(l_seq/BIN)+1
		for i in range(B):
			num1=seq1[i*BIN:(i+1)*BIN].count("N")
			ID1=chrm1+"\t"+str(i*BIN)+"\t"+str((i+1)*BIN)
			dic[ID1]=[num1,0]
	for seq_record in SeqIO.parse(fa_m,"fasta"):
		chrm2 = seq_record.id
		seq2=str(seq_record.seq).strip()
		l_seq=len(str(seq_record.seq).strip())
		if "unanchor" in chrm2: continue
		if l_seq%BIN==0:
			B=int(l_seq/BIN)
		else:
			B=int(l_seq/BIN)+1
		for i in range(B):
			ID2=chrm2+"\t"+str(i*BIN)+"\t"+str((i+1)*BIN)
			num2=seq2[i*BIN:(i+1)*BIN].count("N")
			if ID2 in dic.keys():
				dic[ID2][1]=num2
	for key,v in dic.items():
		k=key.split("\t")
		num=round((int(v[1])-int(v[0]))/BIN*100,2)
		output.write("{0}\t{1}\t{2}\t{3}\n".format(k[0],k[1],k[2],num))
			
def snp_density(vcf, fa):
	chroms=[]
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm = seq_record.id
		if "unanchor" not in chrm:
			chroms.append(chrm)
	out_dic = {}
	for line in vcf:
		if line.startswith("#"): continue
		if "PASS" not in line: continue
		if "unanchor" in line: continue
		llst = line.split("\t")
		chrom, loc = llst[0], int(llst[1])
		if chrom not in chroms: continue
		if chrom not in out_dic:
			out_dic[chrom]=[]
		out_dic[chrom].append(loc)
	return out_dic
def indel_density(vcf, fa):
	chroms=[]
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm = seq_record.id
		if "unanchor" not in chrm:
			chroms.append(chrm)
	ins_dic, del_dic = {}, {}
	for line in vcf:
		if line.startswith("#") or "PASS"  not in line: continue
		if "unanchor" in line: continue
		llst = line.split("\t")
		chrom, loc, ref, alt = llst[0], int(llst[1]), llst[2], llst[3]
		if chrom not in chroms: continue
		if len(ref) < len(alt):
			ins_dic.setdefault(chrom, []).append(loc)
		else:
			del_dic.setdefault(chrom, []).append(loc)
	return ins_dic, del_dic
def out_vcf(out_dic, BIN, var, fa):
	outname = var + '.circos'
	output = open(outname, 'w')
	result = {}
	len_dic= {}
	for seq_record in SeqIO.parse(fa,"fasta"):
		chrm = seq_record.id
		l_seq=len(str(seq_record.seq).strip())
		#print (chrm)
		len_dic[chrm]=l_seq
	for k, v in out_dic.items():
		block_num = int(len_dic[k]/BIN)
		if len_dic[k]%BIN==0:
			block_num=int(len_dic[k]/BIN)
		else:
			block_num=int(len_dic[k]/BIN)+1
		lst = [0] * block_num
		for i in v:
			pos = int(i / BIN)
			lst[pos] += 1
		result[k] = lst
	for k, v in result.items():
		count = 0
		for i in v:
			output.write("{0}\t{1}\t{2}\t{3}\n".format(k,count*BIN+1, (count+1)*BIN, i))
			count +=1
			

gc(fa,BIN)
if argv["c"]:
	conf = open(argv["c"].strip(),'r')
	CON={}
	for line in conf:
		L=line.strip().split("=")
		CON[L[0]]=L[1]
	for i in num_list:
		if i in CON.keys():
			fa = open(argv["fa"].strip(),'r')
			BIN =int(argv["b"].strip())
			gff=open(CON[i],"r")
			TYPE=i
			number(fa,BIN,gff,TYPE)
	if "repeat" in CON.keys():
		fa = open(argv["fa"].strip(),'r')
		BIN =int(argv["b"].strip())
		fa_m=open(CON["repeat"],"r")	
		repeat(fa,fa_m,BIN)
					

	if "snp" in  CON.keys():
		fa = open(argv["fa"].strip(),'r')
		BIN =int(argv["b"].strip())
		vcf=open(CON["snp"],"r")
		out_dic=snp_density(vcf, fa)
		fa = open(argv["fa"].strip(),'r')
		out_vcf(out_dic, BIN, "snp", fa)
	if "indel" in CON.keys():
		fa = open(argv["fa"].strip(),'r')
		BIN =int(argv["b"].strip())
		vcf=open(CON["indel"],"r")
		ins_dic, del_dic = indel_density(vcf, fa)
		fa = open(argv["fa"].strip(),'r')
		out_vcf(del_dic, BIN, "indel", fa)
		fa = open(argv["fa"].strip(),'r')
		out_vcf(ins_dic, BIN, "insert", fa)
