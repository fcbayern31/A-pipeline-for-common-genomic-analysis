import sys
li=[]
dic={}
out=open(sys.argv[2],"w")
with open(sys.argv[3])as ID :
	for i in ID:
		i=i.strip().split("\t")
		dic[i[1]]=[i[0]]+i[2:]

with open(sys.argv[1],'r') as col :
	for i in col :
		i=i.strip().split("\t")
		outlist=dic[i[1]]+dic[i[2]]
		outline='\t'.join(outlist)
		out.write(outline+"\n")
