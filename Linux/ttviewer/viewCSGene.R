#!/usr/bin/env Rscript

library(gggenes)
library(ggplot2)
options (warn = -1)

Args <- commandArgs(trailingOnly=TRUE)
cbbPalette <-c("#000000","#FFFFFF")
geneData = read.csv(Args[1])
subgeneData = read.csv(Args[2])

#out = Args[3]
#name1 = gsub(".gbf","",Args[5])
#name1 = gsub(".gb","",name1)
#resu = paste(out, "/" ,Args[4], "(" ,name1, ")", "_csg.pdf", sep = "")

len_data <- nrow(geneData)
pdf(Args[3], width=6.27,height=0.632*len_data)

geneData$direction <- ifelse(geneData$strand == "forward", 1, -1)
subgeneData$direction <- ifelse(subgeneData$strand == "forward", 1, -1)
gp <- ggplot(geneData, aes(xmin = start, xmax = end, y = reorder(Gene_label,Gene),forward = direction)) + facet_wrap(~ Gene, scales = "free", ncol = 1) + geom_gene_arrow()+ geom_subgene_arrow(data = subgeneData, aes(xmin = start, xmax = end, fill = Subgene, xsubmin = from, xsubmax = to,y = reorder(Gene_label,Gene),forward = direction), color="black") + scale_x_continuous(breaks=c(subgeneData$from,subgeneData$to),labels=c(subgeneData$label_from,subgeneData$label_to)) + theme(axis.ticks.y=element_blank()) + theme(strip.text.x = element_blank()) +  theme(plot.title = element_text(hjust = 0.5,lineheight=0.8)) + scale_fill_manual(values = cbbPalette) + theme(axis.text.x=element_text(size = 7, face="bold.italic")) + theme(legend.title = element_text(colour='white'))+ theme(axis.text.y =element_text(size = 9,color = "black",face="bold")) + theme(panel.grid.major.y =element_blank(),panel.grid.minor=element_blank())


if(len_data == 1){gp + ylab("  ")}else if(len_data == 2 | len_data == 3){gp + scale_y_discrete(name="\nCis-splicing Genes\n")+ theme(axis.title.y=element_text(size=8))}else{gp + scale_y_discrete(name="\nCis-splicing Genes\n")+ theme(axis.title.y=element_text(size=12))}


dev.off()
print("---  Pdf file has been generated, please check the folder under the current file!  ---")
