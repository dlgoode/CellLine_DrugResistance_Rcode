### Redo splines analysis with
### 1) Zero count genes removed *prior* to taking Rlogs
### 2) genes ranked according to sign

library(GSEA)

All.sets.GSVA.signed.scores <- gsva( Trimmed.CellLine.TimeCourse.Rlogs, All.gene.sets[[1]][1:28], method="ssgsea", rnaseq=T, abs.ranking=F )

Nutlin.spline.results <- topTable( fit.778.Nutlin, num=dim(Combined.sets.GSVA.scores)[1] )
Tunicamycin.spline.results <- topTable( fit.778.Tunicamycin, num=dim(Combined.sets.GSVA.scores)[1] )
Erlotinib.spline.results <- topTable( fit.778.Erlotinib, num=dim(Combined.sets.GSVA.scores)[1] )

fit.778.Nutlin <- lmFit(Combined.sets.GSVA.scores[,c(1:4,9:13)], design.778.Nutlin)
fit.778.Tunicamycin <- lmFit(Combined.sets.GSVA.scores[,c(1:4,14:18)], design.778.Tunicamycin)
fit.778.Erlotinib <- lmFit(Combined.sets.GSVA.scores[,c(23:31)], design.778.Erlotinib)


### Function to color strings based on content.

color.by.string <- function( string, search, color1="black", color2="red") {

  if ( grepl( search, string ) ) {

    return( color2 )
  }

  return( color1 )
}
