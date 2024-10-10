#### Perform Stouffer's Z score meta-analysis on differential expression results

library(jaffelab)
library(collapse)
library(ggplot2)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(stringr)

genes = c("Akap11", "Dagla", "Gria3", "Grin2a", "Sp4", "Srrm2", "Zmym2")

#define region
reg = "STR"

#load all individual mutant DE tables, add Z score
exps = genes
genenames = list()
allde = list()
for (exp in exps){
  gene = exp

  geno = "HT"
  
  de = read.csv("path/to/DE/table.csv")
  de = de[!is.na(de$padj),]
  de$gene = exp

  #calculate Z from p-value
  de$Z = qnorm(de$pvalue/2) * sign(de$log2FoldChange) * -1
  if(sum(is.infinite(de$Z)) > 0){
    de$Z[is.infinite(de$Z) & de$Z > 0] = max(de$Z[is.finite(de$Z)])
    de$Z[is.infinite(de$Z) & de$Z < 0] = min(de$Z[is.finite(de$Z)])
  }
  allde[[exp]] = de
  genenames[[exp]] = de$X
}

#reduce to genes that are tested in all experiments
#make Z score matrix
int = Reduce(intersect, genenames)
zs = lapply(allde, function(x){
  dat = x[match(int, x$X),]
  return(dat$Z)
})

zs = do.call("cbind", zs)
rownames(zs) = int

#get sample size for all experiments (# HT + # WT)
sampsize = read.csv("/path/to/bulkSampleSize.csv")
rownames(sampsize) = sampsize$Gene
ns = sampsize[colnames(zs),]$n

#calculate weights
w = sqrt(ns)
denom = sqrt(sum(w^2))

#multiply weight through
weighted_zs = zs %r*% w
#calculate meta Z
stouffers = rowSums(weighted_zs)/denom

#convert meta Z to p-value
metap = 2*pnorm(abs(stouffers), lower.tail = FALSE)
dat = as.data.frame(cbind(stouffers, metap))
dat$FDR = p.adjust(dat$metap, "fdr")
dat = dat[order(dat$FDR),]

write.csv(dat, file = "path/to/output.csv")
