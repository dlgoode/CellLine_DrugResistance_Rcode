
> getwd()
[1] "C:/Users/goode david/Documents/Peter Mac/Arcadi_htseq_data"

### create sample info data frame. 
Arcadi.sample.files <- list.files("htseq_Files/")

Arcadi.sample.names <- vapply( list.files("htseq_Files/"),FUN=function(lt) { strsplit(lt, ".s")[[1]][1] } , FUN.VAL="nothing" )

cell.type.Arcadi.sample <- vapply( list.files("htseq_Files/"),FUN=function(lt) { strsplit(lt, "_|-")[[1]][1] } , FUN.VAL="nothing" )


Arcadi.sample.status <- vapply( Arcadi.sample.files,FUN=function(lt) { ifelse ( grepl("parental",lt ), "Parental", "Resistant") }, FUN.VAL="nothing" )

Arcadi.sample.table <- cbind( "SampleName"=Arcadi.sample.names, "SampleFile"=Arcadi.sample.files, "CellLine"=cell.type.Arcadi.sample, "Status"=Arcadi.sample.status)

### Load data into DESeq2 object:
library(DESeq2)

Arcadi.DESeq2 <- DESeqDataSetFromHTSeqCount( Arcadi.sample.table, directory = "htseq_Files", ~1)

min.allgenes.unlogged <- apply( unnormed.Arcadi.timecourse.counts, 1, min)
max.allgenes.unlogged <- apply( unnormed.Arcadi.timecourse.counts, 1, max)

fraction.zero.counts <- apply( unnormed.Arcadi.timecourse.counts, 2, FUN=function(counts) { mean(counts>0)})

Rlog.Arcadi.DESeq2 <- rlog(Arcadi.DESeq2, blind=T )

> CellLine.TimeCourse.Rlogs <- assay( Rlog.Arcadi.DESeq2 ) 
> dim( CellLine.TimeCourse.Rlogs)
[1] 50775    39

min.allgenes <- apply( CellLine.TimeCourse.Rlogs, 1, min)
max.allgenes <- apply( CellLine.TimeCourse.Rlogs, 1, max)

### Trim out no count genes:
Trimmed.CellLine.TimeCourse.Rlogs <- CellLine.TimeCourse.Rlogs[ min.allgenes!=0 | max.allgenes!=0, ]
> dim( Trimmed.CellLine.TimeCourse.Rlogs)
[1] 39528    39

### Check by PCA:
Rlog.PCA <- prcomp( CellLine.TimeCourse.Rlogs )

plot( Rlog.PCA$rotation[order.samples.by.date,1], Rlog.PCA$rotation[order.samples.by.date,2], 
        pch=CellLine.Shapes, frame=F, cex=1.1, col=colours.by.Expt , xlab="PC1", ylab="PC2", xlim=c(0.155,0.167) )

legend( 0.163, 0.19, legend=Expt.Names, col=unique(colours.by.Expt), pch=c(rep(15,4),rep(c(17,16,18,8),each=2) ), bty="n", cex=0.75 )

### Run GSVA
library(GSVA)

### Test sets first
b<-list(c("ATM","ATR","APC"),c("TP53","BRCA1","BRCA2") )

gsva( Trimmed.CellLine.TimeCourse.Rlogs, b, method="ssgsea" )

### Try with one of Arcadi's lists.
### Imported from Excel in the format given, with one column per set,
### with lists of genes in the set. 

### Can convert to list on the fly.
DNA.repair.gene.sets <- list(read.table("DNA_repair_GeneSets",sep="\t",header=T, fill=T ) )

DNA.repair.GSVA.scores <- gsva( Trimmed.CellLine.TimeCourse.Rlogs, DNA.repair.gene.sets[[1]][1:19], method="ssgsea", rnaseq=T, abs.ranking=T )

DNA.repair.GSVA.scores.unnormed <- gsva( Trimmed.CellLine.TimeCourse.Rlogs, DNA.repair.gene.sets[[1]][1:19], method="ssgsea", rnaseq=T, abs.ranking=T, ssgsea.norm=F )

order.samples.by.date <- c(1,4,2,3,5,6,8,7,9,10,13,12,11,14,15,18,17,16,23,26,24,25,27,28,31,30,29,19:22,32:39)
plot( DNA.repair.GSVA.scores.unnormed[1,order.samples.by.date], ylab="ssGSEA score", frame=F)

### Factor by cell line type for point shape:
CellLine.Factor <- as.numeric( as.factor( Arcadi.sample.table[order.samples.by.date,3] ) )

> summary( as.factor( Arcadi.sample.table[order.samples.by.date,3] ) )
778   AU565   H1975   SKBR3 SKMEL28 
18       4       9       4       4 

shapes.by.CellLine <- c( rep(0,18),rep(1,9),rep(c(2,5,6),each=4))

colours.by.Expt <- c( rep(c("skyblue","dodgerblue"),each=4), rep(c("blue","navyblue"),each=5),
                    rep("green",each=4),rep("darkgreen",each=5),rep("red",each=2),rep("darkred",each=2),
                    rep("violet",2), rep("purple",2), rep("orange",2), rep("orangered",2) )
                    

plot( DNA.repair.GSVA.scores.unnormed[1,order.samples.by.date], ylab="ssGSEA score", frame=F, pch=shapes.by.CellLine)

plot( DNA.repair.GSVA.scores[1,order.samples.by.date], ylab="ssGSEA score", frame=F, pch=shapes.by.CellLine,xaxt="n")

Expt.Names=c("778-Parental","778-CDK4i","778-Nutlin","778-Tunicamycin","H1975-Parental","H1975-Erlotinib", "AU565-Parental", "AU565-Lapatinib", "SKBR3-Parental", "SKBR3-Lapatinib", "SKMEL-Parental", "SKMEL-Vemirafinib")

CellLine.Shapes <- c(CellLine.Factor[-c(36:39)]+14,rep(8,4) )

## Put altogether in function can use for each set:
for( i in 1:19 ) {
  
  Scores <- DNA.repair.GSVA.scores[i,]

  min.max <- summary(Scores)[c(1,6)]

  Set <- rownames(DNA.repair.GSVA.scores)[i]

  filename = paste(Set,"pdf",sep=".")

  pdf(filename,width=10)

  plot( Scores[order.samples.by.date], ylab="ssGSEA score", frame=F, pch=CellLine.Shapes, xlim=c(0,49),
         xaxt="n", col=colours.by.Expt, ylim=c(min.max[1]-0.1, min.max[2]+0.1), xlab="", main=Set )

  abline(v=c(18.5,27.5,31.5,35.5),col="grey80", lty=2)

  legend( 40,min.max[2]+0.1, legend=Expt.Names, col=unique(colours.by.Expt), pch=c(rep(15,4),rep(c(17,16,18,8),each=2) ), bty="n", cex=0.75 )

  axis(side=1,at=c(9,22.5,29.5,33.5,37.5),col="white",labels=c("778","H1975","AU565","SKBR3","SKMEL3"), cex=0.7, las=2 )

  dev.off()
}

### Group different cell lines together.
order.together <- c( rep(c(1:4),2), rep(1:5,2), c(7:10), c(7:11), rep(c(13:14),2),rep(c(16:17),2),rep(c(19:20),2) )

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

### Test for trend over time in an exptl group:

plot( order.together[1:18], Scores[order.samples.by.date][1:18], ylab="Scaled ssGSEA score", frame=F, pch=CellLine.Shapes, xlim=c(0,8), las=2,
              xaxt="n", col=colours.by.Expt, ylim=c(min.max[1]-0.1, min.max[2]), xlab="", main=Set )

text( 5.5, Scores[order.samples.by.date][4], cor.test( order.together[1:4], Scores[order.samples.by.date][1:4], method="sp" )$est, col=colours.by.Expt[1] )
lines( 1:4, lm( Scores[order.samples.by.date][1:4] ~ c(1:4) )$fitted.values, col=colours.by.Expt[1], lwd=2 )
rd.med <- round( median( Scores[order.samples.by.date][1:4]), 2)
text( 0.2, rd.med, rd.med, col=colours.by.Expt[1] )
lines( c(0.5,4.5), rep( rd.med , 2 ), col=colours.by.Expt[1], lty=3 )

text( 5.5, Scores[order.samples.by.date][8], cor.test( order.together[1:4+4], Scores[order.samples.by.date][1:4+4], method="sp" )$est, col=colours.by.Expt[5] )
lines( 1:4, lm( Scores[order.samples.by.date][1:4+4] ~ c(1:4) )$fitted.values, col=colours.by.Expt[5], lwd=2 )
rd.med <- round( median( Scores[order.samples.by.date][5:8]), 2)
text( 0.2, rd.med, rd.med, col=colours.by.Expt[5] )
lines( c(0.5,4.5), rep( rd.med , 2 ), col=colours.by.Expt[5], lty=3 )


text( 5.5, Scores[order.samples.by.date][13], cor.test( order.together[9:13], Scores[order.samples.by.date][9:13], method="sp" )$est, col=colours.by.Expt[9] )
lines( 1:5, lm( Scores[order.samples.by.date][9:13] ~ c(1:5) )$fitted.values, col=colours.by.Expt[9], lwd=2 )
rd.med <- round( median( Scores[order.samples.by.date][9:13]), 2)
text( 0.2, rd.med, rd.med, col=colours.by.Expt[9] )
lines( c(0.5,4.5), rep( rd.med , 2 ), col=colours.by.Expt[9], lty=3 )

text( 5.5, Scores[order.samples.by.date][18], labels=cor.test( order.together[9:13+5], Scores[order.samples.by.date][9:13+5], method="sp" )$est, col=colours.by.Expt[14] )
lines( 1:5, lm( Scores[order.samples.by.date][14:18] ~ c(1:5) )$fitted.values, col=colours.by.Expt[14], lwd=2 )
rd.med <- round( median( Scores[order.samples.by.date][9:13+5]), 2)
text( 0.2, rd.med, rd.med, col=colours.by.Expt[14] )
lines( c(0.5,4.5), rep( rd.med , 2 ), col=colours.by.Expt[14], lty=3 )

legend( 6, 1.8, legend=Expt.Names[1:4], col=unique(colours.by.Expt), pch=rep(15,4), bty="n", cex=0.75 )
