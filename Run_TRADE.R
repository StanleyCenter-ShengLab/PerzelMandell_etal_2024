#### Run TRADE analysis for each experiment
library(DESeq2)
library(sva)
library(ashr)
library(ggplot2)
library("tximport")
library(jaffelab)
library(TRADE)
`%notin%` <- Negate(`%in%`)

set.seed(333)
#load salmon alignment data
load_salmon_data = function(gene, age, reg) {
  path = paste0("/path/to/", gene, "/bulkRNAseq/Alignment/", age, "/")
  samples = list.files(path)

  #get meta data
  region = ss(samples, "\\_")
  genotype = ss(samples, "\\_", 2)
  colData = as.data.frame(cbind(samples, region, genotype))
  
  #rename HT and WT to perturb and control
  colData$condition = "control"
  colData$condition[colData$genotype == "HT"] = "perturb"
  colData = colData[colData$region == reg,] 
  
  filePath <- paste0("/path/to/", gene, "/bulkRNAseq/Alignment/", age, "/")
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
  
  counts_matrix = txi.salmon$counts
  
  filter_perturb = rowSums(counts_matrix[,colData$condition == "perturb"]) != 0
  filter_control = rowSums(counts_matrix[,colData$condition == "control"]) != 0
  counts_matrix = counts_matrix[filter_perturb & filter_control & sapply(1:nrow(counts_matrix), function(row) all(counts_matrix[row,] > 0)),]
  counts_matrix = round(counts_matrix)
  
  return(list(counts = counts_matrix,
              metadata = colData))
}


run_SVA_DESeq2 <- function(counts,metadata) {
  #Make the DESeq2 Object
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
  
  #This block of code runs SVA to identify the number of surrogate variables to include
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized = TRUE)
  idx <- rowMeans(dat) > 1
  dat <- dat[idx,]
  mod <- model.matrix(~ condition, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  n.sv = num.sv(dat,mod,method="be", seed = 12321902)
  
  if (n.sv > 0) {
    #Compute the SV values
    svseq <- svaseq(dat, mod = mod, mod0 = mod0, n.sv = n.sv)
    sv_matrix = svseq$sv
    colnames(sv_matrix) = sapply(1:n.sv, function(x) paste0("SV",x))
    metadata_sv = cbind(metadata, sv_matrix)
    
    #build the SVs into the DESeq2 formula
    dds_sva = DESeqDataSetFromMatrix(countData = counts, colData = metadata_sv, design = as.formula(paste("~",paste(c(colnames(sv_matrix),"condition"), collapse = "+"))))
    
    pre_analysis_dds = dds_sva
    
    #Run DESeq2
    dds_sva <- DESeq(dds_sva,betaPrior = FALSE)
    
    
    results_condition = results(dds_sva,contrast = c("condition","perturb","control"))
    
  } else {
    pre_analysis_dds = dds
    dds <- DESeq(dds,
                 betaPrior = FALSE)
    
    
    results_condition = results(dds,contrast = c("condition","perturb","control"))
    
  }
  
  return(list(dds = pre_analysis_dds,
              results = results_condition))
  
}

#This function does the permutation testing
get_perm_varmean_data <- function(dds, n_perm) {
  #Randomly resample columns of the counts matrix
  sample = sapply(1:n_perm, function(x) sample(1:ncol(dds)))
  
  perm_stats <- lapply(1:n_perm,
                       function(iter) {
                         print(iter)
                         dds_perm = dds
                         #This is where the permuting actually happens
                         dds_perm$condition = dds$condition[sample[,iter]]
                         #Print out the number of times when the permuted condition label is the same as the true one
                         #This was a helpful way to QC that permuting runs as expected
                         print(table(dds_perm$condition == dds$condition))
                         print(table(dds_perm$condition))
                         
                         #Run SVA/DESeq2
                         perm_output = run_SVA_DESeq2(assay(dds_perm,"counts"),colData(dds_perm))
                         results = perm_output$results
                         
                         #Print the mean squared effect size, useful for QC
                         print(sqrt(mean(results$log2FoldChange^2, na.rm = TRUE)))
                         
                         #Run TRADE, print the TRADE variance (i.e. transcriptome-wide impact)
                         TRADE_perm = TRADE(mode = "univariate", results1 = results)
                         print(TRADE_perm$distribution_summary$transcriptome_wide_impact)
                         return(TRADE_perm)
                       })
  return(list(samples = sample,
              perm_stats = perm_stats))
  
}

#This final wrapper function is basically a efficient way to call the functions above
#Avoids repeating many lines of code for the different regions/drugs
get_TRADE_output <- function(gene,
                                       age,
                                       reg,
                                       n_perm) {
  
  dataset <- load_salmon_data(gene, age, reg)
  
  #Run SVA/DESeq2 on the total dataset
  output <- run_SVA_DESeq2(dataset$counts,dataset$metadata)
  
  #Run TRADE on the total dataset
  TRADE_output = TRADE(mode = "univariate", results1 = output$results)
  
  #Permutation analysis
  permstats <- get_perm_varmean_data(output$dds,n_perm)
  
  #Grab the transcriptome-wide impact estimates from the permutation TRADE output
  perm_vars <- sapply(permstats$perm_stats,
                      function(output) {
                        output$distribution_summary$transcriptome_wide_impact
                      })
  return(list(deseq2_output = output,
              TRADE_output = TRADE_output,
              permstats = permstats,
              perm_vars = perm_vars))
  
}

MUTANT_AGE_REG = get_TRADE_output("[gene]", "[age]", "[reg]", 1000)
MUTANT_AGE_REG_TRADE = MUTANT_AGE_REG$TRADE_output
MUTANT_AGE_REG_TRADE_pval = mean(MUTANT_AGE_REG$perm_vars > MUTANT_AGE_REG$TRADE_output$distribution_summary$transcriptome_wide_impact)


