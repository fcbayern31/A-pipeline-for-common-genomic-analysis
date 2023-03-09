library(pheatmap)


df = read.delim("draw.matrix", header=False ,sep = "\t")
pdf("Hi-C_inter.pdf")
pheatmap(df, color = colorRampPalette(c("#FFEFD5","#8B0000")),
	 annotation_legend=TRUE, 
	 border_color=NA, 
	 scale="none" )
dev.off()
