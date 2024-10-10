library(MatrixEQTL)
library(egg)
library(ggplot2)
library(cowplot)
library(stringr)

#### load data ####
reg1 = "TH"

exps = c("Akap11", "Dagla", "Gria3","Sp4", "Grin2a", "Srrm2", "Zmym2")

for (exp in exps){
  gene = exp
  geno = "HT"
  
  go = read.csv("path/to/GSEA/output.csv")
  if(gene == "Akap11"){
    dat = data.frame(go$NES)
    rownames(dat) = go$desc
    p = data.frame(go$pval < 0.05)
    rownames(dat) = go$desc
  } else {
    mat = match(rownames(dat), go$desc)
    dat[,gene] = go$NES[mat]
    p[,gene] = go$pval[mat] < 0.05
  }
}
colnames(dat)[1] = "Akap11"
p[is.na(p)] = FALSE
#filter to gene sets which are nominally significant in at least 3 mutants
dat = dat[rowSums(p) > 2,]

#filter to gene sets which are synapse related
dat = dat[rownames(dat) %in% c(unique(rownames(dat)[grepl("synap", rownames(dat))]),
                               unique(rownames(dat)[grepl("axon", rownames(dat))]),
                               unique(rownames(dat)[grepl("dendrit", rownames(dat))])),]

reg1_dat = dat

#define second region
reg2 = "SSC"

exps = c("Akap11", "Dagla", "Gria3","Sp4", "Grin2a", "Srrm2", "Zmym2")

for (exp in exps){
  gene = exp
  geno = "HT"
  
  go = read.csv("path/to/GSEA/output.csv")
  if(gene == "Akap11"){
    dat = data.frame(go$NES)
    rownames(dat) = go$desc
    p = data.frame(go$pval < 0.05)
    rownames(dat) = go$desc
  } else {
    mat = match(rownames(dat), go$desc)
    dat[,gene] = go$NES[mat]
    p[,gene] = go$pval[mat] < 0.05
  }
}
colnames(dat)[1] = "Akap11"
p[is.na(p)] = FALSE

#filter to gene sets which are nominally significant in at least 3 mutants
dat = dat[rowSums(p) > 2,]
reg2_dat = dat


#### correlation ####
library(jaffelab)
rownames(reg1_dat)= paste0(reg1, "_", rownames(reg1_dat))
rownames(reg2_dat)= paste0(reg2, "_", rownames(reg2_dat))

dat = rbind(reg1_dat, reg2_dat)
#correlate all rows to each other
x = cor(t(dat), method = "spearman") #or pearson
x[upper.tri(x, diag = TRUE)] <- NA
x <- na.omit(reshape::melt(t(x)))
x = as.data.frame(x)
x$X1 = as.character(x$X1)
x$X2 = as.character(x$X2)
x$reg1 = ss(x$X1, "\\_")
x$reg2 = ss(x$X2, "\\_")
#remove correlations within the same region
x$keep = x$reg1 != x$reg2
c = x[x$keep,]

c = c[order(abs(c$value), decreasing = TRUE),]

#calculate correlation p-val
spearmant <- function(r,n){r*sqrt((n-2)/(1-r^2))}
c$t = spearmant(c$value, 7)
c$p = 2*pt(q = abs(c$t), df=6, lower.tail=FALSE)
c$fdr = p.adjust(c$p, "fdr")
sum(c$fdr < 0.05)
c = c[c$fdr < 0.05,]

