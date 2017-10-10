num.Main.sets <- dim( Main.sets.GSVA.scores)[1]

cor.btwn.gene.sets <- matrix( nrow=num.Main.sets, ncol=num.Main.sets )

for( i in 1:num.Main.sets ) { 
  
  for( j in 1:num.Main.sets ) {
    
    cor.btwn.gene.sets[i,j] <- cor.test( Main.sets.GSVA.scores[i,], Main.sets.GSVA.scores[j,], method="sp")$est
  }
}

rownames( cor.btwn.gene.sets ) <- substr( rownames(Main.sets.GSVA.scores), 1, 10 )
colnames( cor.btwn.gene.sets ) <- substr( rownames(Main.sets.GSVA.scores), 1, 10 )

pheatmap( cor.btwn.gene.sets )

### Now do for the other gene sets

num.other.Gene.sets <- dim(All.sets.GSVA.scores)[1]

cor.btwn.other.gene.sets <- matrix( nrow=num.other.Gene.sets, ncol=num.other.Gene.sets )

for( i in 1:num.other.Gene.sets ) { 
  
  for( j in 1:num.other.Gene.sets ) {
    
    cor.btwn.other.gene.sets[i,j] <- cor.test( All.sets.GSVA.scores[i,], All.sets.GSVA.scores[j,], method="sp")$est
  }
}

rownames( cor.btwn.other.gene.sets ) <- substr( rownames(All.sets.GSVA.scores), 1, 19 )
colnames( cor.btwn.other.gene.sets ) <- substr( rownames(All.sets.GSVA.scores), 1, 19 )

pheatmap( cor.btwn.other.gene.sets )