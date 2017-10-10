### Script for running GSVA on all time points and plotting ssGSEA values
### across all time points for each gene set.

All.gene.sets <- list(read.table("All_Gene_Sets",sep="\t",header=T, fill=T ) )

library(GSVA)

All.sets.GSVA.scores <- gsva( Trimmed.CellLine.TimeCourse.Rlogs, All.gene.sets[[1]][1:28], method="ssgsea", rnaseq=T, abs.ranking=T )

All.sets.GSVA.scores.unnormed <- gsva( Trimmed.CellLine.TimeCourse.Rlogs, All.gene.sets[[1]][1:28], method="ssgsea", rnaseq=T, abs.ranking=T, ssgsea.norm=F )

for( i in 1:27 ) {
  
  Scores <- All.sets.GSVA.scores[i,]
  
  min.max <- summary(Scores)[c(1,6)]
  
  Set <- rownames(All.sets.GSVA.scores)[i]
  
  filename = paste(Set,"_grouped",".pdf",sep="")
  
  pdf(filename,width=10.5)
  
  plot( order.together, Scores[order.samples.by.date], ylab="Scaled ssGSEA score", frame=F, pch=CellLine.Shapes, xlim=c(0,28), las=2,
        xaxt="n", col=colours.by.Expt, ylim=c(min.max[1]-0.1, min.max[2]+0.1), xlab="", main=Set )
  
  abline(v=c(6,12,15,18),col="grey80", lty=2)
  
  legend( 21,min.max[2]+0.1, legend=Expt.Names, col=unique(colours.by.Expt), pch=c(rep(15,4),rep(c(17,16,18,8),each=2) ), bty="n", cex=0.75 )
  
  axis(side=1,at=c(2.5,8.5,13.5,16.5,19.5),col="white",labels=c("778","H1975","AU565","SKBR3","SKMEL3"), cex=0.7, las=2 )
  
  axis(side=1,at=c(1:5,7:11,13:14,16:17,19:20),col="white", col.ticks="grey67", labels=F, cex=0.7, las=2 )
  
  dev.off()
  
}

r <-  read.csv("DNA damage Signaling Pathway.csv", head=F, row.names=1)

l2 <- list(as.data.frame(t(r)))
> length(l2[[1]])
[1] 30

library(GSVA)
Main.sets.GSVA.scores <- gsva( Trimmed.CellLine.TimeCourse.Rlogs, l2[[1]][1:length(l2[[1]])], method="ssgsea", rnaseq=T, abs.ranking=T )
rownames(Main.sets.GSVA.scores) <-  gsub("/","-",rownames(Main.sets.GSVA.scores))

for( i in 1:length(l2[[1]]) ) {
  
  Scores <- Main.sets.GSVA.scores[i,]
  
  min.max <- summary(Scores)[c(1,6)]
  
  Set <- rownames(Main.sets.GSVA.scores)[i]
  
  filename = paste(Set,"_grouped",".pdf",sep="")
  
  pdf( paste( "Main Gene Sets/",filename, sep=""),width=10)
  
  plot( order.together, Scores[order.samples.by.date], ylab="Scaled ssGSEA score", frame=F, pch=CellLine.Shapes, xlim=c(0,28), las=2,
        xaxt="n", col=colours.by.Expt, ylim=c(min.max[1]-0.1, min.max[2]+0.1), xlab="", main=Set )
  
  abline(v=c(6,12,15,18),col="grey80", lty=2)
  
  legend( 21,min.max[2]+0.1, legend=Expt.Names, col=unique(colours.by.Expt), pch=c(rep(15,4),rep(c(17,16,18,8),each=2) ), bty="n", cex=0.75 )
  
  axis(side=1,at=c(2.5,8.5,13.5,16.5,19.5),col="white",labels=c("778","H1975","AU565","SKBR3","SKMEL3"), cex=0.7, las=2 )
  
  axis(side=1,at=c(1:5,7:11,13:14,16:17,19:20),col="white", col.ticks="grey67", labels=F, cex=0.7, las=2 )
  
  dev.off()
  
}

### Combine all gene sets into one list:
Combined.gene.sets <- list(read.table("Combined_Gene_Sets.txt",sep="\t",header=T, fill=T ) , as.data.frame(t(r)))
