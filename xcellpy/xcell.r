#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

library(GSEABase)
library(xCell)
rnaseq=as.logical(toupper(args[1]))
scale=as.logical(toupper(args[2]))
alpha=as.numeric(args[3])
nperm=as.numeric(args[4])
parallel.sz=as.numeric(args[5])
verbose=as.logical(toupper(args[6]))
tempdir=args[7]
beta_pval=as.logical(toupper(args[8]))
perm_pval=as.logical(toupper(args[9]))

mat = read.csv(file.path(tempdir,'expr.csv'),header=TRUE,row.names=1,check.names=FALSE)

output = xCellAnalysis(mat,
	    rnaseq=rnaseq,
	    scale=scale,
	    alpha=alpha,
        parallel.sz=parallel.sz
    )
ofile = file.path(tempdir,"pathways.csv")
write.csv(output,ofile)

if(beta_pval) {
    distroBeta = xCellSignifcanceBetaDist(output,rnaseq=rnaseq)
    ofileBeta = file.path(tempdir,"beta.csv")
    write.csv(distroBeta,ofileBeta)
}
if(perm_pval) {
    distroRandom = xCellSignifcanceRandomMatrix(scores = output, expr = mat, spill = NULL, alpha=alpha,nperm=nperm)
    ofileRandomP = file.path(tempdir,"randomP.csv")
    ofileRandomD = file.path(tempdir,"randomD.csv")
    write.csv(distroRandom$pvals,ofileRandomP)
    write.csv(distroRandom$shuff_xcell,ofileRandomD)
}
