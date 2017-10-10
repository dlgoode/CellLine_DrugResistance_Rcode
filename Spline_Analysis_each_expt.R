X.778.CDK4i <- ns(splines.targets$Time[1:8], df=2 )
Group.778.CDK4i  <- factor(splines.targets$Group[1:8] )

design.778.CDK4i  <- model.matrix(~Group.778.CDK4i*X.778.CDK4i)

fit.778.CDK4i <- lmFit(Combined.sets.GSVA.scores[,1:8], design.778.CDK4i)

fit.778.CDK4i <- eBayes(fit.778.CDK4i)

CDK4i.spline.results <- topTable( fit.778.CDK4i, num=dim(Combined.sets.GSVA.scores)[1] )

write.table( CDK4i.spline.results, file="CDK4i.spline.results.txt", sep="\t", quote=FALSE )

### Nutlin
X.778.Nutlin <- ns(splines.targets$Time[c(1:4,9:13)], df=2 )
Group.778.Nutlin  <- factor(splines.targets$Group[c(1:4,9:13)] )

design.778.Nutlin  <- model.matrix(~Group.778.Nutlin*X.778.Nutlin)

fit.778.Nutlin <- lmFit(Combined.sets.GSVA.scores[,c(1:4,9:13)], design.778.Nutlin)

fit.778.Nutlin <- eBayes(fit.778.Nutlin)

Nutlin.spline.results <- topTable( fit.778.Nutlin, num=dim(Combined.sets.GSVA.scores)[1] )

write.table( Nutlin.spline.results, file="Nutlin.spline.results.txt", sep="\t", quote=FALSE )

### Tunicamycin
X.778.Tunicamycin <- ns(splines.targets$Time[c(1:4,14:18)], df=2 )
Group.778.Tunicamycin  <- factor(splines.targets$Group[c(1:4,14:18)] )

design.778.Tunicamycin  <- model.matrix(~Group.778.Tunicamycin*X.778.Tunicamycin)

fit.778.Tunicamycin <- lmFit(Combined.sets.GSVA.scores[,c(1:4,14:18)], design.778.Tunicamycin)

fit.778.Tunicamycin <- eBayes(fit.778.Tunicamycin)

Tunicamycin.spline.results <- topTable( fit.778.Tunicamycin, num=dim(Combined.sets.GSVA.scores)[1] )

write.table( Tunicamycin.spline.results, file="Tunicamycin.spline.results.txt", sep="\t", quote=FALSE )

### Erlotinib (H1975)
X.778.Erlotinib <- ns(splines.targets$Time[c(23:31)], df=2 )
Group.778.Erlotinib  <- factor(splines.targets$Group[c(23:31)] )

design.778.Erlotinib  <- model.matrix(~Group.778.Erlotinib*X.778.Erlotinib)

fit.778.Erlotinib <- lmFit(Combined.sets.GSVA.scores[,c(23:31)], design.778.Erlotinib)

fit.778.Erlotinib <- eBayes(fit.778.Erlotinib)

Erlotinib.spline.results <- topTable( fit.778.Erlotinib, num=dim(Combined.sets.GSVA.scores)[1] )

write.table( Erlotinib.spline.results, file="Erlotinib.spline.results.txt", sep="\t", quote=FALSE )
