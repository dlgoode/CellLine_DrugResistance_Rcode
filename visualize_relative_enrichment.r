plot(c( rep(c(1:4),2), rep(c(1:5),2)), Trimmed.CellLine.TimeCourse.Rlogs[ rownames(Trimmed.CellLine.TimeCourse.Rlogs)=="TP53I11",1:18][ order.samples.by.date[1:18]], frame=F )

ordered.gene.set.names <- rownames(CDK4i.spline.results)[ order(rownames(CDK4i.spline.results)) ]

pheatmap( -log10( cbind(CDK4i.spline.results[order(rownames(CDK4i.spline.results)),9], 
                        Nutlin.spline.results[ order(rownames(Nutlin.spline.results)),9],
                        
                  Tunicamycin.spline.results[order(rownames(Tunicamycin.spline.results)),9],
                  Erlotinib.spline.results[order(rownames(Erlotinib.spline.results)),9]
                  ) ),
          labels_row = ordered.gene.set.names, scale='column', cluster_rows =  F, fontsize_row=7 )

annotation_row = data.frame(ordered.gene.set.names)[,1],

### combine adjusted p-values.
spline.log10.adj.pvalues <- 
  -log10( cbind(CDK4i.spline.results[order(rownames(CDK4i.spline.results)),9], 
                Nutlin.spline.results[ order(rownames(Nutlin.spline.results)),9],
              
              Tunicamycin.spline.results[order(rownames(Tunicamycin.spline.results)),9],
              Erlotinib.spline.results[order(rownames(Erlotinib.spline.results)),9]
) )

jpeg("pheatmap.clustered.adj.pvalues.jpg", height=740, width=600)
pheatmap( spline.log10.adj.pvalues, labels_row = ordered.gene.set.names, cluster_cols=F,
                         labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"), fontsize_row=7, fontsize_col = 9 )
dev.off()

spline.log10.nom.pvalues <- 
  -log10( cbind(CDK4i.spline.results[order(rownames(CDK4i.spline.results)),8], 
                Nutlin.spline.results[ order(rownames(Nutlin.spline.results)),8],
                
                Tunicamycin.spline.results[order(rownames(Tunicamycin.spline.results)),8],
                Erlotinib.spline.results[order(rownames(Erlotinib.spline.results)),8]
  ) )

pheatmap( spline.log10.nom.pvalues, labels_row = ordered.gene.set.names, cluster_cols=F,
            labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"), scale='column', cluster_rows =  F, fontsize_row=7 ) 

spline.log10.ranked.pvalues <- apply( spline.log10.nom.pvalues, 2, rank )

pheatmap( spline.log10.ranked.pvalues, labels_row = ordered.gene.set.names, cluster_cols=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"), cluster_rows =  F, fontsize_row=7 ) 

### Order gene sets by combined rank of p-value across all samples:
sum.rank.per.set <- apply( spline.log10.ranked.pvalues, 1, sum )

pheatmap( spline.log10.ranked.pvalues[ order(sum.rank.per.set, decreasing=TRUE), ], labels_row = ordered.gene.set.names[ order(sum.rank.per.set, decreasing=TRUE) ], cluster_cols=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"),  fontsize_row=7, cluster_rows = F )

### Look at the change in expression (Resistant vs Parental) of
### all gene sets, ordered by p-value rank:
spline.exprn.diffs <- 
  cbind(CDK4i.spline.results[order(rownames(CDK4i.spline.results)),1], 
              Nutlin.spline.results[ order(rownames(Nutlin.spline.results)),1],
              
              Tunicamycin.spline.results[order(rownames(Tunicamycin.spline.results)),1],
              Erlotinib.spline.results[order(rownames(Erlotinib.spline.results)),1]
) 

pheatmap( spline.exprn.diffs[ order(sum.rank.per.set, decreasing=TRUE), ], labels_row = ordered.gene.set.names[ order(sum.rank.per.set, decreasing=TRUE) ], cluster_cols=F,
          labels_col=c("778_CDK4i","778_Nutlin","778_Tunicamycin", "H1975_Erlotinib"),  fontsize_row=7, cluster_rows = F )

### Compare median ssGSEA difference for AU565 and other samples with only 2 time points.
summary( apply( All.sets.GSVA.scores[,grep( "AU565_R", colnames( All.sets.GSVA.scores) )], 1, mean ) 
         - apply( All.sets.GSVA.scores[,grep( "AU565_p", colnames( All.sets.GSVA.scores) )], 1, mean ) )

### Collate into one data frame for plotting
Median.Exprn.Diffs <- cbind( 
     "778-R-CDK4i"= apply( Combined.sets.GSVA.scores[,grep( "CDK4i", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "778-par", colnames( Combined.sets.GSVA.scores) )], 1, median ),
     "778-R-Nutlin"= apply( Combined.sets.GSVA.scores[,grep( "Nutlin", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "778-par", colnames( Combined.sets.GSVA.scores) )], 1, median ),
     "778-R-Tunicamycin"= apply( Combined.sets.GSVA.scores[,grep( "Tunicamycin", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "778-par", colnames( Combined.sets.GSVA.scores) )], 1, median ) )

Median.Exprn.Diffs <- cbind(Median.Exprn.Diffs,
   "H1975-R-Erlotinib"=apply( Combined.sets.GSVA.scores[,grep( "Erlot", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "H1975-par", colnames( Combined.sets.GSVA.scores) )], 1, median ) 
)

Median.Exprn.Diffs  <- cbind(Median.Exprn.Diffs,
                            "AU565_R_lapatinib"=apply( Combined.sets.GSVA.scores[,grep( "AU565_R", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "AU565_par", colnames( Combined.sets.GSVA.scores) )], 1, median ), 
                            "SKBR3_R_lapatinib"=apply( Combined.sets.GSVA.scores[,grep( "SKBR3_R", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "SKBR3_par", colnames( Combined.sets.GSVA.scores) )], 1, median ), 
                            "SKMEL28_R_vem"=apply( Combined.sets.GSVA.scores[,grep( "SKMEL28_R", colnames( Combined.sets.GSVA.scores) )], 1, median ) - apply( Combined.sets.GSVA.scores[,grep( "SKMEL28_par", colnames( Combined.sets.GSVA.scores) )], 1, median ) 
                            
)

summary( Median.Exprn.Diffs)

pheatmap( Median.Exprn.Diffs, cluster_rows=F )

### Want to order the Expression Diffs in same order as ranked by sum
### Do this by finding position of each rowname in the row names vector ordered by rank sum.
### Need to remove redundant 'mTOR.signaling' entries first.

Median.Exprn.Diffs <- Median.Exprn.Diffs[-43,]

reordered.names.for.median.exprn <- vapply(rownames(Median.Exprn.Diffs),
                                     FUN=function(X) { grep( X, ordered.gene.set.names[ order(sum.rank.per.set, decreasing=TRUE) ][-51])}, FUN.VAL=0)

pheatmap( Median.Exprn.Diffs[order(reordered.names.for.median.exprn),], cluster_rows=F, cluster_col=F, fontsize_row = 7 )
