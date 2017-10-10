
library(edgeR)

#input files, labels and groups
htseqFilesswap<-c("sample6_sortedByname.out.txt","sample15_sortedByname.out.txt","sample17_sortedByname.out.txt","sample8_sortedByname.out.txt","sample16_sortedByname.out.txt","sample18_sortedByname.out.txt")
varClass_swap<-c("Benign", "Benign","Benign","Cancer","Cancer","Cancer")
varLabel_swap<-c("Bening6", "Bening15","Benign17","Cancer8","Cancer16","Cancer18")
cbind(htseqFilesswap,varClass_swap,varLabel_swap)
datahtseq_swap<-readDGE(files=htseqFilesswap,group=varClass_swap,labels=varLabel_swap)

#Filter tags
filter_swap<-grepl("_",rownames(datahtseq_swap),)
filteredDatahtseq_swap<-datahtseq_swap[!filter_swap,]

#Remove median zeros across samples
filteredDatahtseqMedianZeros_swap<-apply(filteredDatahtseq_swap,1,median)==0
filteredDatahtseqRemovedMedianZeros_swap<-filteredDatahtseq_swap[!filteredDatahtseqMedianZeros_swap,]

#Normalization for the library sizes by default it used TMM method
filteredDatahtseqRemovedMedianZerosNorm_swap<-calcNormFactors(filteredDatahtseqRemovedMedianZeros_swap)
col.cell<-c("green","red")[filteredDatahtseqRemovedMedianZerosNorm_swap$samples$group]

# MDS PLOT
plotMDS(filteredDatahtseqRemovedMedianZerosNorm_swap,col=col.cell,cex=0.5)
dev.off()

#Pathway analysis using SSGSEA and MSigDB resource
library(GSVA)
library(GSEABase)
msigDB52<-getGmt("../../../data/CancerHallmarks50Pathwaysh.v6.0WithAR-NEC-GeneSets.symbols.gmt",sep="\t")
filteredDatahtseqRemovedMedianZerosNorm_EnrichmentScore52<-gsva(filteredDatahtseqRemovedMedianZerosNorm_swap$counts,msigDB52,method="ssgsea",rnaseq=TRUE,mx.diff=TRUE,abs.ranking=TRUE)
#Generate heatmap using pheatmap 
library(pheatmap)
library(RColorBrewer)
colr<-colorRampPalette(c("blue","white","red"))(256)
pheatmap(filteredDatahtseqRemovedMedianZerosNorm_EnrichmentScore52,fontsize=8,color=colr)


