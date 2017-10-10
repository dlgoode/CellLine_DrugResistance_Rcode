### Splines analysis for all genes, with 778 design.
### Start with 778-R-CDK4i vs 778-Parental
preTrimmed.gene.fit.778.CDK4i <- lmFit(preTrimmed.TimeCourse.Rlogs [,1:8], design.778.CDK4i)
preTrimmed.gene.fit.778.CDK4i <- eBayes(preTrimmed.gene.fit.778.CDK4i)
all.gene.topTable.778.CDK4i <- topTable( preTrimmed.gene.fit.778.CDK4i, num=39528 )

preTrimmed.gene.fit.778.Nutlin <- lmFit(preTrimmed.TimeCourse.Rlogs[,c(1:4,9:13)], design.778.Nutlin)
preTrimmed.gene.fit.778.Nutlin <- eBayes(preTrimmed.gene.fit.778.Nutlin)
all.gene.topTable.778.Nutlin <- topTable( preTrimmed.gene.fit.778.Nutlin, num=39528 )

preTrimmed.gene.fit.778.Tunicamycin <- lmFit(preTrimmed.TimeCourse.Rlogs[,c(1:4,14:18)], design.778.Tunicamycin)
preTrimmed.gene.fit.778.Tunicamycin <- eBayes(preTrimmed.gene.fit.778.Tunicamycin)
all.gene.topTable.778.Tunicamycin <- topTable( preTrimmed.gene.fit.778.Tunicamycin, num=39528 )

preTrimmed.gene.fit.H1975.Erlotinib <- lmFit(preTrimmed.TimeCourse.Rlogs[,c(23:31)], design.778.Erlotinib)
preTrimmed.gene.fit.H1975.Erlotinib <- eBayes(preTrimmed.gene.fit.H1975.Erlotinib)
all.gene.topTable.H1975.Erlotinib <- topTable( preTrimmed.gene.fit.H1975.Erlotinib, num=39528 )

### Get the names of all the gene sets:
names( Combined.gene.sets[[1]])

gene=rownames(all.gene.topTable.778.CDK4i)[2]

exprn.vals <- preTrimmed.TimeCourse.Rlogs[ rownames(preTrimmed.TimeCourse.Rlogs)==gene, 1:8]

plot( order.together[1:8], exprn.vals, frame=F,
      main=paste("Time course expression of ",gene,sep=""), col=colours.by.Expt[1:8], pch=CellLine.Shapes[1:8], 
      xaxt="n", las=2, ylab="Normalized (Rlog) Expression", xlab="Time point", ylim=c(min(exprn.vals)-0.1,max(exprn.vals)+0.1) 
)
lines( order.together[1:4], exprn.vals[1:4], type="l", col=colours.by.Expt[1] )
lines( order.together[5:8], exprn.vals[5:8], type="l", col=colours.by.Expt[5] )

### Plot 
for( i in 1:19 ) {
  
  Scores <- DNA.repair.GSVA.scores[i,]
  
  min.max <- summary(Scores)[c(1,6)]
  
  Set <- rownames(DNA.repair.GSVA.scores)[i]
  
  filename = paste(Set,"_grouped",".pdf",sep="")
  
  #  pdf(filename,width=10)
  
  plot( order.together, Scores[order.samples.by.date], ylab="Scaled ssGSEA score", frame=F, pch=CellLine.Shapes, xlim=c(0,27), las=2,
        xaxt="n", col=colours.by.Expt, ylim=c(min.max[1]-0.1, min.max[2]+0.1), xlab="", main=Set )
  
  abline(v=c(6,12,15,18),col="grey80", lty=2)
  
  legend( 21,min.max[2], legend=Expt.Names, col=unique(colours.by.Expt), pch=c(rep(15,4),rep(c(17,16,18,8),each=2) ), bty="n", cex=0.75 )
  
  axis(side=1,at=c(2.5,8.5,13.5,16.5,19.5),col="white",labels=c("778","H1975","AU565","SKBR3","SKMEL3"), cex=0.7, las=2 )
  
  axis(side=1,at=c(1:5,7:11,13:14,16:17,19:20),col="white", col.ticks="grey67", labels=F, cex=0.7, las=2 )
  
  #    dev.off()
  
}
dim(Combined.gene.sets[[1]][1])