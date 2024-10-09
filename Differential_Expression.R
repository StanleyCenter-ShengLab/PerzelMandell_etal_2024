library(tximport)
library(edgeR)
library(jaffelab)
library(DESeq2)

#get samples and metadata
path = "path/to/Alignment/folders/"
samples = list.files(path)

region = ss(samples, "\\_")
genotype = ss(samples, "\\_", 2)
colData = as.data.frame(cbind(samples, region, genotype))

#select one brain region at a time to analyze
reg = "PFC"
colData = colData[colData$region == reg,] 

###read in Salmon
filePath <- "/path/to/Alignment/folders/"
files <- paste0(filePath, colData$samples, "/Salmon_out/quant.sf")
names(files) <- colData$samples
tx2gene = read.table(paste0(filePath, colData$samples[1], "/Salmon_out/trans2gene.txt"))
txi.salmon = tximport(files, type = "salmon", tx2gene = tx2gene)

txi.salmon$abundance <-
  txi.salmon$abundance[apply(txi.salmon$length,1,function(row) all(row > 0 )),]
txi.salmon$counts <-
  txi.salmon$counts[apply(txi.salmon$length,1,function(row) all(row > 0 )),]
txi.salmon$length <-
  txi.salmon$length[apply(txi.salmon$length,1,function(row) all(row > 0 )),]

### run DESeq2
dds <- DESeqDataSetFromTximport(txi.salmon, colData, ~ genotype)
dds$genotype <- relevel(dds$genotype, "WT")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

check <- estimateSizeFactors(dds)
norm <- counts(check, normalized=TRUE)

#check PCA for any issues
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="samples")
plotPCA(vsdata, intgroup="genotype")

#HT v WT
resHT <- results(dds, contrast=c("genotype", "HT", "WT"))
resHT <- lfcShrink(dds=dds, contrast=c("genotype", "HT", "WT"), res=resHT, type="normal")
sum(resHT$padj < 0.05, na.rm = TRUE)
write.csv(resHT, file = "path/to/output.csv")



