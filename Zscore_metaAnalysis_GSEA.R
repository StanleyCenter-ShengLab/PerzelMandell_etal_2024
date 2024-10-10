#### Perform Stouffer's Z score meta-analysis on GSEA results

library(jaffelab)
library(collapse)
library(reshape2)
library(ggplot2)
library(rutils)
library(stringr)
library(cowplot)
library(RColorBrewer)

genes = c("Akap11", "Dagla", "Gria3", "Grin2a", "Sp4", "Srrm2", "Zmym2")

#define region
reg = "STR"

#load each mutant's GSEA results, add Z score
exps = genes
pathnames = list()
allgo = list()
for (exp in exps){
  gene = exp
  geno = "HT"
  
  go = read.csv("path/to/GSEA/output.csv")
  go$gene = exp

  #calculate Z from p-val
  go$Z = qnorm(go$pval/2) * sign(go$NES) * -1
  # if(sum(is.infinite(de$Z)) > 0){
  #   de$Z[is.infinite(de$Z) & de$Z > 0] = max(de$Z[is.finite(de$Z)])
  #   de$Z[is.infinite(de$Z) & de$Z < 0] = min(de$Z[is.finite(de$Z)])
  # }
  allgo[[exp]] = go
  pathnames[[exp]] = go$desc
}

#reduce to gene sets tested in all mutants, extract Z score matrix
int = Reduce(intersect, pathnames)
zs = lapply(allgo, function(x){
  dat = x[match(int, x$desc),]
  return(dat$Z)
})

zs = do.call("cbind", zs)
rownames(zs) = int

#load samples sizes (# HT + # WT)
sampsize = read.csv("path/to/bulkSampleSize.csv")
rownames(sampsize) = sampsize$Gene
ns = sampsize[colnames(zs),]$n

#calculate weights
w = sqrt(ns)
denom = sqrt(sum(w^2))

#multiply weight through
weighted_zs = zs %r*% w

#calculate stouffer's z
stouffers = rowSums(weighted_zs)/denom

#calculate meta-p from meta-z
metap = 2*pnorm(abs(stouffers), lower.tail = FALSE)
dat = as.data.frame(cbind(stouffers, metap))
dat$FDR = p.adjust(dat$metap, "fdr")
dat = dat[order(dat$FDR),]

write.csv(dat, file = "path/to/output.csv")


#### Plot top GSEA meta-analysis results
# select only signficant gene sets
dat = dat[dat$FDR < 0.05,]
dat$X = rownames(dat)
fullgo = do.call("rbind", allgo)

#collapse redundant GO terms
sub = fullgo[match(unique(fullgo$pathway), fullgo$pathway),]
ont = sub$ont[match(dat$X, sub$desc)]
id = sub$pathway[match(dat$X, sub$desc)]

pathway_df = as.data.frame(cbind(ont, id))
colnames(pathway_df) = c("go_type", "go_id")

#use meta-z as score to prioritize gene sets in collapsing
scores = abs(dat$stouffers)
names(scores) = pathway_df$go_id

red = go_reduce(pathway_df, "org.Mm.eg.db", threshold = 0.9, scores = scores)

#filter results to "parent terms"
plotdat = dat[dat$X %in% red$parent_term,]

#### PLOT ####
#top 20 gene sets
plotdat = plotdat[1:20,]
plotdat = plotdat[order(plotdat$stouffers),]
#get individual mutant Z scores
subzs = zs[plotdat$X,]
meltzs = melt(subzs)

colnames(meltzs)[2] = "Mutant"

p <- ggplot(plotdat, aes(x=stouffers, y = X)) +
  geom_bar(stat = "identity", fill = "gray75") +
  scale_y_discrete(limits=plotdat$X, labels = function(description) str_wrap(description, width = 40)) +
  theme_minimal_grid(12) +
  theme(axis.text=element_text(size=13), legend.text=element_text(size=13), legend.title=element_text(size=13)) +
  geom_point(data = meltzs, aes(x = value, y = Var1, fill =Mutant), color = "black", pch = 21, size = 2.4) +   
  scale_fill_manual(values=brewer.pal(7, "Accent")) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", color = "gray44") +
  ggtitle("Top meta-analysis gene sets") + xlab("Z") + ylab("")

p
