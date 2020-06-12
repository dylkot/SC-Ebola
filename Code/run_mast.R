suppressMessages(library(MAST))

args <- commandArgs(trailingOnly=TRUE)
expressionfn <- args[1]
cdatfn <- args[2]
gdatfn <- args[3]
outfn <- args[4]
ncores <- args[5]
testvar <- args[6]
covars <- strsplit(args[7], ",")[[1]]
options(mc.cores=ncores)


## Prepare formula##
formula <- paste(covars, collapse=" + ")
formula <- paste(testvar, formula, sep=" + ")
formula <- paste("~ ", formula, sep="")
formula <- as.formula(formula)
print(formula)

## Load data and convert to SingleCellExperiment ##
expr <- read.csv(expressionfn,  row.names=1)
cdat <- read.csv(cdatfn, row.names=1)
gdat <- read.csv(gdatfn, row.names=1)
scaRaw <- FromMatrix(as.matrix(expr), cdat, gdat)

## Run MAST ##
zlm.output <- zlm(formula, scaRaw)
res = summary(zlm.output, doLRT=testvar)$datatable

## Collect output ##
ind <- res$contrast == testvar
logFC <- res[(res$component == "logFC") & ind,c("primerid", "coef")]
H <- res[(res$component == "H") & ind,c("primerid", "Pr(>Chisq)")]
C <- res[(res$component == "C") & ind,c("primerid", "Pr(>Chisq)")]
D <- res[(res$component == "D") & ind,c("primerid", "Pr(>Chisq)")]
colnames(logFC) <- c('primerid', 'logFC')
colnames(H) <- c('primerid', 'P_H')
colnames(C) <- c('primerid', 'P_C')
colnames(D) <- c('primerid', 'P_D')
merged = merge(logFC, H, by='primerid', all=T, suffixes=c('D', 'C'))
merged = merge(merged, C, by='primerid', all=T, suffixes=c('D', 'C'))
merged = merge(merged, D, by='primerid', all=T, suffixes=c('D', 'C'))

## Save to output file ##
write.table(merged, file=outfn, sep=",")

