import sys
dic={}
li=[]
out=open(sys.argv[2],"w")
with open(sys.argv[3]) as cl :
	for i in cl :
		i=i.strip().split("\t")
		dic[i[0]]=i[1]

with open(sys.argv[1],"r") as col :
	for i in col:
		i=i.strip().split("\t")
		outline="\t".join(i)+"\t"+"color="+dic[i[0]]+",thickness=4"
		out.write(outline+"\n")
