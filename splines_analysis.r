expt_names<-apply( Arcadi.sample.table, 1, FUN=function(x) { y=strsplit(x[1], "-[0-9]|_[0-9]", perl=T ); return(unlist(y)[1]) } )

### Get the experiment each sample belongs to
Expt_Groups <- as.vector( unlist( lapply( strsplit( Arcadi.sample.table[,2], split="_0|_1|_2|-0|-1|-2", perl=T ), FUN=function(x) { x[1] } ) ) )


lapply( strsplit( Arcadi.sample.table[1:4,1], "-"), FUN=function(x) { x[4] } )

### Order samples within experiment by date of collection:
rank( unlist( lapply( strsplit( Arcadi.sample.table[5:8,1], "-"), FUN=function(x) { x[5] } ) ) )
778-R-CDK4i-03-01-2014.sorted.txt 778-R-CDK4i-04-02-2014.sorted.txt 
                                1                                 2 
778-R-CDK4i-16-05-2014.sorted.txt 778-R-CDK4i-24-03-2014.sorted.txt 
                                4                                 3
                                
### Make a vector of relative order of each experimental set, for splines analysis
> for( i in 2:length(Expt_Groups) ) {  
+   if( Expt_Groups[i] != Expt_Groups[i-1] ) { print(i) }
+ }
[1] 5
[1] 9
[1] 14
[1] 19
[1] 21
[1] 23
[1] 27
[1] 32
[1] 34
[1] 36
[1] 38
                                
j=1                                
Relative.orders <- rank( unlist( lapply( strsplit( Arcadi.sample.table[j:(j+3),1], "-"), FUN=function(x) { x[4] } ) ) )

j=5
Relative.orders <- c( Relative.orders, rank( unlist( lapply( strsplit( Arcadi.sample.table[j:(j+3),1], "-"), FUN=function(x) { x[5] } ) ) ) )

j=9
Relative.orders <- c( Relative.orders, rank( unlist( lapply( strsplit( Arcadi.sample.table[j:(j+4),1], "-"), FUN=function(x) { x[5] } ) ) ) )

j=14
Relative.orders <- c( Relative.orders, rank( unlist( lapply( strsplit( Arcadi.sample.table[j:(j+4),1], "-"), FUN=function(x) { x[5] } ) ) ) )

### Next sets are AU565_parental & AU565_R_lapatinib, have date listed differently.
Relative.orders <- c( Relative.orders, c(2,1) )
Relative.orders <- c( Relative.orders, c(1,2) )

j=23 
Relative.orders <- c( Relative.orders, rank( unlist( lapply( strsplit( Arcadi.sample.table[j:(j+3),1], "-"), FUN=function(x) { x[4] } ) ) ) )

j=27
Relative.orders <- c( Relative.orders, rank( unlist( lapply( strsplit( Arcadi.sample.table[j:(j+4),1], "-"), FUN=function(x) { x[5] } ) ) ) )

Relative.orders <- c( Relative.orders, 2, 1, 1, 2, 2, 1, 1, 2)

### Try Splines analysis
### Following example on p 49 of limma user guide
library(splines)

## create the targets matrix:
#splines.targets <- as.data.frame( cbind( "FileName"=Arcadi.sample.table[,2], "Group"=Expt_Groups, "Time"=as.numeric(as.vector(Relative.orders) ) ) )
#splines.targets <- data.frame( "FileName"=Arcadi.sample.table[,2], "Group"=Expt_Groups, "Time"=Relative.orders )

X <- ns(splines.targets$Time, df=3)

Group <- factor(splines.targets$Group)

design <- model.matrix(~Group*X)

### But problems doing the analysis itself:
> fit <- lmFit(CellLine.TimeCourse.Rlogs, design)
Coefficients not estimable: GroupAU565_parental:X4 GroupAU565_R_lapatinib:X4 GroupH1975-parental:X4 GroupH1975-R-Erlotinib:X4 GroupSKBR3_parental:X4 GroupSKBR3_R_lapatinib:X4 GroupSKMEL28_parental:X4 GroupSKMEL28_R_vem:X4 Group778-R-CDK4i:X5 Group778-R-Nutlin:X5 Group778-R-Tunicamycin:X5 GroupAU565_parental:X5 GroupAU565_R_lapatinib:X5 GroupH1975-parental:X5 GroupH1975-R-Erlotinib:X5 GroupSKBR3_parental:X5 GroupSKBR3_R_lapatinib:X5 GroupSKMEL28_parental:X5 GroupSKMEL28_R_vem:X5 X5 GroupAU565_parental:X2 GroupAU565_R_lapatinib:X2 GroupSKBR3_parental:X2 GroupSKBR3_R_lapatinib:X2 GroupSKMEL28_parental:X2 GroupSKMEL28_R_vem:X2 GroupAU565_parental:X3 GroupAU565_R_lapatinib:X3 GroupSKBR3_parental:X3 GroupSKBR3_R_lapatinib:X3 GroupSKMEL28_parental:X3 GroupSKMEL28_R_vem:X3 Group778-R-CDK4i:X4

> fit <- eBayes(fit)
Error in ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  : 
                  No residual degrees of freedom in linear model fits
                
### I think the samples with data from only 2 time points are not valid for this analysis.

### Try with just parental and CDK4i 778 samples, as test
X <- ns(splines.targets$Time[1:8], df=2)

Group <- factor(splines.targets$Group[1:8])

design <- model.matrix(~Group*X)

fit <- lmFit(Trimmed.CellLine.TimeCourse.Rlogs[,1:8], design)

fit <- eBayes(fit)

### Works!

> topTable(fit)
Removing intercept from test coefficients
Group778.R.CDK4i          X1          X2 Group778.R.CDK4i.X1 Group778.R.CDK4i.X2   AveExpr        F      P.Value  adj.P.Val
UBD              6.117322  2.24281088 -0.91487839          -5.5012789          0.46617924  3.619338 6950.882 6.036977e-06 0.04580578
TP53I11          4.049904  0.01448683  0.02867828          -0.8465387         -1.88815714  5.373945 5963.541 7.478964e-06 0.04580578
JAKMIP2          6.142798  0.04384864  0.09382021           1.6115981         -0.09768742  3.541536 5738.971 7.891277e-06 0.04580578
SLC16A9          3.323312 -1.20496008  0.03160345           3.1808353         -0.27787398  4.891106 5719.249 7.929348e-06 0.04580578
ELSPBP1          2.463520  0.04359734  0.09423654           8.0533721         -4.51710401  0.821178 5697.974 7.970770e-06 0.04580578
PTPRT            4.062169  0.03727801  0.08057718           0.3364585          3.34334928  0.905696 5227.196 8.992060e-06 0.04580578
CYP24A1         -3.769114  0.26505643  0.02726145          -1.7309296          1.26899298 11.445900 5099.188 9.309210e-06 0.04580578
ADAMTS19         4.085897 -1.66242535  0.86911794           2.7216883         -2.76837062  2.396959 4296.092 1.182940e-05 0.04580578
SCUBE1           5.543535  0.04193620  0.09055874          -3.3856066         -3.18344885  1.558758 3431.251 1.619687e-05 0.04580578
PLIN4            2.019467  0.58589816  0.20985965           3.3034147         -2.55462882  3.904861 3212.182 1.776168e-05 0.04580578

### Just have to redo with all samples, and with all 778 samples.
### What about if there are differing number of time points?
### Also, degrees of freedom depends on how many different conditions there are? 

### Can I do for each experiment separately?

X.778 <- ns(splines.targets$Time[1:18], df=3)

Group.778 <- factor(splines.targets$Group[1:18])

design.778 <- model.matrix(~Group.778*X.778)
View(design.778)

fit.778 <- lmFit(Trimmed.CellLine.TimeCourse.Rlogs[,1:18], design.778)

fit.778 <- eBayes(fit.778)

topTable( fit.778 )


### Plot time course for all 778 experiments:
gene <- 'TP53'
plot(c( rep(c(1:4),2), rep(c(1:5),2)), Trimmed.CellLine.TimeCourse.Rlogs[ rownames(Trimmed.CellLine.TimeCourse.Rlogs)==gene,1:18][ order.samples.by.date[1:18]],
       frame=F, col=c(rep("blue",4),rep("red",4),rep("green",5),rep("purple",5) ), ylab="Rlog value", xlab="Time Point")


#########
#########
### Try with gene sets instead of genes, using ssGSEA scores
### (as described in All_gene_set_analysis.r )
###
#########

X.778.CDK4i <- ns(splines.targets$Time[1:8], df=2 )
Group.778.CDK4i  <- factor(splines.targets$Group[1:8] )

design.778.CDK4i  <- model.matrix(~Group.778.CDK4i*X.778.CDK4i)

fit.778.CDK4i <- lmFit(All.sets.GSVA.scores[,1:8], design.778.CDK4i)

fit.778.CDK4i <- eBayes(fit.778.CDK4i)

topTable( fit.778.CDK4i, num=27 )

set1 <- rownames( topTable( fit.778.CDK4i ) )[5]

## Fits for the other time courses:
fit.778.Nutlin <- lmFit(Combined.sets.GSVA.scores[,c(1:4,9:13)], design.778.Nutlin)
fit.778.Tunicamycin <- lmFit(Combined.sets.GSVA.scores[,c(1:4,14:18)], design.778.Tunicamycin)
fit.778.Erlotinib <- lmFit(Combined.sets.GSVA.scores[,c(23:31)], design.778.Erlotinib)


### Save results to workspace
CDK4i.spline.results <- topTable( fit.778.CDK4i, num=dim(Combined.sets.GSVA.scores)[1] )

### Likewise for other conditions:
Nutlin.spline.results <- topTable( fit.778.Nutlin, num=dim(Combined.sets.GSVA.scores)[1] )
Tunicamycin.spline.results <- topTable( fit.778.Tunicamycin, num=dim(Combined.sets.GSVA.scores)[1] )
Erlotinib.spline.results <- topTable( fit.778.Erlotinib, num=dim(Combined.sets.GSVA.scores)[1] )


### Plot the normalized ssGSEA value at each time point for the two conditions:
plot( rep(c(1:4),2) , Combined.sets.GSVA.scores[rownames(Combined.sets.GSVA.scores)==set1,1:8][ order.samples.by.date[1:8]],
     frame=F, col=c(rep("skyblue",4),rep("navyblue",4)), ylab="Rlog value", xlab="Time Point", 
       pch=c(rep(15,4),rep(17,4)), ylim=c(1.2,1.4), main=set1, xaxt="n") ##

axis( 1, at=c(1:4))
points( c(1:4), Combined.sets.GSVA.scores[rownames(Combined.sets.GSVA.scores)==set1,1:8][ order.samples.by.date[1:4]], col="skyblue", type="l" )
points( c(1:4), Combined.sets.GSVA.scores[rownames(Combined.sets.GSVA.scores)==set1,1:8][ order.samples.by.date[5:8]], col="navyblue", type="l" )

legend( 1, 1.4, bty="n", legend=c("Parental","CDK4i"), col=c("skyblue","navyblue"), pch=c(15,17) )

