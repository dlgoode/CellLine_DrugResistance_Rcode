### Retry the time course analysis 
### 1) Zero count genes removed *prior* to taking Rlogs
### 2) genes ranked according to sign

> trimmed.unnormed.Arcadi.timecourse.counts <- unnormed.Arcadi.timecourse.counts[  max.allgenes.unlogged > 0 , ]
> dim( trimmed.unnormed.Arcadi.timecourse.counts )
[1] 39528    39

trimmed.Arcadi.DESeq2 <- DESeqDataSetFromMatrix( countData = trimmed.unnormed.Arcadi.timecourse.counts,
                                                 colData=data.frame(colnames(trimmed.unnormed.Arcadi.timecourse.counts)), ~1)

Rlog.trimmed.Arcadi.DESeq2 <- rlog(trimmed.Arcadi.DESeq2, blind=T )

preTrimmed.TimeCourse.Rlogs <- assay( Rlog.trimmed.Arcadi.DESeq2 ) 

### Check by PCA:
preTrimmed.Rlog.PCA <- prcomp( preTrimmed.TimeCourse.Rlogs )

plot( preTrimmed.Rlog.PCA$rotation[order.samples.by.date,1], preTrimmed.Rlog.PCA$rotation[order.samples.by.date,2], 
      pch=CellLine.Shapes, frame=F, cex=1.1, col=colours.by.Expt , xlab="PC1", ylab="PC2", xlim=c(0.155,0.167) )

library( GSVA )
library(splines)

preTrimmed.Combined.sets.GSVA.scores <- gsva( preTrimmed.TimeCourse.Rlogs, Combined.gene.sets[[1]][1:length(Combined.gene.sets[[1]])], method="ssgsea", rnaseq=T, abs.ranking=T )

### Visualize selected sets to see if AU565 reversed now:
i <- grep("CIN",rownames(preTrimmed.Combined.sets.GSVA.scores))

Scores <- preTrimmed.Combined.sets.GSVA.scores[i,]

min.max <- summary(Scores)[c(1,6)]

Set <- rownames(preTrimmed.Combined.sets.GSVA.scores)[i]

plot( order.together, Scores[order.samples.by.date], ylab="Scaled ssGSEA score", frame=F, pch=CellLine.Shapes, xlim=c(0,28), las=2,
      xaxt="n", col=colours.by.Expt, ylim=c(min.max[1]-0.1, min.max[2]+0.1), xlab="", main=Set )

### Still looks the same as before!
### Distributions of ssGSEA scores still very similar
boxplot( Combined.sets.GSVA.scores )
boxplot( preTrimmed.Combined.sets.GSVA.scores )

######
### Repeat spline analysis with new sets anyway 
preTrimmed.fit.778.CDK4i <- lmFit(preTrimmed.Combined.sets.GSVA.scores[,1:8], design.778.CDK4i)

preTrimmed.fit.778.CDK4i <- eBayes(preTrimmed.fit.778.CDK4i)

X.778.Nutlin <- ns(splines.targets$Time[c(1:4,9:13)], df=2 )
Group.778.Nutlin  <- factor(splines.targets$Group[c(1:4,9:13)] )

design.778.Nutlin  <- model.matrix(~Group.778.Nutlin*X.778.Nutlin)

preTrimmed.fit.778.Nutlin <- lmFit(preTrimmed.Combined.sets.GSVA.scores[,c(1:4,9:13)], design.778.Nutlin)
preTrimmed.fit.778.Nutlin <- eBayes(preTrimmed.fit.778.Nutlin)

preTrimmed.fit.778.Tunicamycin <- lmFit(preTrimmed.Combined.sets.GSVA.scores[,c(1:4,14:18)], design.778.Tunicamycin)
preTrimmed.fit.778.Tunicamycin <- eBayes(preTrimmed.fit.778.Tunicamycin)

preTrimmed.fit.778.Erlotinib <- lmFit(preTrimmed.Combined.sets.GSVA.scores[,c(23:31)], design.778.Erlotinib)
preTrimmed.fit.778.Erlotinib <- eBayes(preTrimmed.fit.778.Erlotinib)

######

### Save results to workspace
CDK4i.preTrimmed.spline.results <- topTable( preTrimmed.fit.778.CDK4i, num=dim(preTrimmed.Combined.sets.GSVA.scores)[1] )
Nutlin.preTrimmed.spline.results <- topTable( preTrimmed.fit.778.Nutlin, num=dim(preTrimmed.Combined.sets.GSVA.scores)[1] )
Tunicamycin.preTrimmed.spline.results <- topTable( preTrimmed.fit.778.Tunicamycin, num=dim(preTrimmed.Combined.sets.GSVA.scores)[1] )
Erlotinib.preTrimmed.spline.results <- topTable( preTrimmed.fit.778.Erlotinib, num=dim(preTrimmed.Combined.sets.GSVA.scores)[1] )

### Now plot the results:
library( pheatmap )

preTrimmed.spline.adj.pvalues <- 
  cbind(CDK4i.preTrimmed.spline.results[order(rownames(CDK4i.preTrimmed.spline.results)),9], 
                Nutlin.preTrimmed.spline.results[ order(rownames(Nutlin.preTrimmed.spline.results)),9],
                
                Tunicamycin.preTrimmed.spline.results[order(rownames(Tunicamycin.preTrimmed.spline.results)),9],
                Erlotinib.preTrimmed.spline.results[order(rownames(Erlotinib.preTrimmed.spline.results)),9]
  )

rownames(preTrimmed.spline.adj.pvalues) <- ordered.gene.set.names2 
colnames(preTrimmed.spline.adj.pvalues) <- c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib")


preTrimmed.spline.log10.adj.pvalues <- 
  -log10( cbind(CDK4i.preTrimmed.spline.results[order(rownames(CDK4i.preTrimmed.spline.results)),9], 
                Nutlin.preTrimmed.spline.results[ order(rownames(Nutlin.preTrimmed.spline.results)),9],
                
                Tunicamycin.preTrimmed.spline.results[order(rownames(Tunicamycin.preTrimmed.spline.results)),9],
                Erlotinib.preTrimmed.spline.results[order(rownames(Erlotinib.preTrimmed.spline.results)),9]
  ) )

rownames(preTrimmed.spline.log10.adj.pvalues) <- ordered.gene.set.names2 
colnames(preTrimmed.spline.log10.adj.pvalues) <- c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib")

ordered.gene.set.names2 <- rownames(CDK4i.preTrimmed.spline.results)[ order(rownames(CDK4i.preTrimmed.spline.results)) ]

jpeg("pheatmap.clustered.adj.pvalues.Sept7.jpg", height=790, width=600)

pheatmap( preTrimmed.spline.log10.adj.pvalues, labels_row = ordered.gene.set.names2, cluster_cols=F, cluster_rows=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"), fontsize_row=8, fontsize_col = 9 )

dev.off()

### Nominal p-values
preTrimmed.spline.log10.nom.pvalues <- 
  -log10( cbind(CDK4i.preTrimmed.spline.results[order(rownames(CDK4i.preTrimmed.spline.results)),8], 
                Nutlin.preTrimmed.spline.results[ order(rownames(Nutlin.preTrimmed.spline.results)),8],
                
                Tunicamycin.preTrimmed.spline.results[order(rownames(Tunicamycin.preTrimmed.spline.results)),8],
                Erlotinib.preTrimmed.spline.results[order(rownames(Erlotinib.preTrimmed.spline.results)),8]
  ) )

pheatmap( preTrimmed.spline.log10.nom.pvalues, labels_row = ordered.gene.set.names2, cluster_cols=F, cluster_rows=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"), fontsize_row=7, fontsize_col = 9 )

### Compare the ranks of the sets overall:
### Order gene sets by combined rank of p-value across all samples:
preTrimmed.spline.log10.ranked.pvalues <- apply( preTrimmed.spline.log10.nom.pvalues, 2, rank )

pheatmap( preTrimmed.spline.log10.ranked.pvalues, labels_row = ordered.gene.set.names, cluster_cols=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"), cluster_rows =  F, fontsize_row=7 ) 


preTrimmed.sum.rank.per.set <- apply( preTrimmed.spline.log10.ranked.pvalues, 1, sum )

jpeg("pheatmap.order.by.cum.rank.nom.pvalues.Sept7.jpg", height=790, width=600)

pheatmap( preTrimmed.spline.log10.ranked.pvalues[ order(preTrimmed.sum.rank.per.set, decreasing=TRUE), ], labels_row = ordered.gene.set.names2[ order(preTrimmed.sum.rank.per.set, decreasing=TRUE) ], cluster_cols=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"),  fontsize_row=8, cluster_rows = F )

dev.off()

### Look at the change in expression (Resistant vs Parental) of
### all gene sets, ordered by p-value rank:
preTrimmed.spline.exprn.diffs <- 
  cbind(CDK4i.preTrimmed.spline.results[order(rownames(CDK4i.preTrimmed.spline.results)),1], 
        Nutlin.preTrimmed.spline.results[ order(rownames(Nutlin.preTrimmed.spline.results)),1],
        
        Tunicamycin.preTrimmed.spline.results[order(rownames(Tunicamycin.preTrimmed.spline.results)),1],
        Erlotinib.preTrimmed.spline.results[order(rownames(Erlotinib.preTrimmed.spline.results)),1]
  ) 

jpeg("pheatmap.diff.exprn.order.by.cum.rank.Sept7.jpg", height=790, width=600)

pheatmap( preTrimmed.spline.exprn.diffs[ order(preTrimmed.sum.rank.per.set, decreasing=TRUE), ], labels_row = ordered.gene.set.names2[ order(preTrimmed.sum.rank.per.set, decreasing=TRUE) ], cluster_cols=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"),  fontsize_row=8, cluster_rows = F )

dev.off()

### Compare expression differences:
preTrimmed.Median.Exprn.Diffs <- cbind( 
  "778-R-CDK4i"= apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "CDK4i", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ) - apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "778-par", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ),
  "778-R-Nutlin"= apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "Nutlin", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ) - apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "778-par", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ),
  "778-R-Tunicamycin"= apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "Tunicamycin", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ) - apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "778-par", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ) )

preTrimmed.Median.Exprn.Diffs <- cbind(preTrimmed.Median.Exprn.Diffs,
                            "H1975-R-Erlotinib"=apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "Erlot", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ) - apply( preTrimmed.Combined.sets.GSVA.scores[,grep( "H1975-par", colnames( preTrimmed.Combined.sets.GSVA.scores) )], 1, median ) 
)