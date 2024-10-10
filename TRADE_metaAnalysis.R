#### Meta-Analyze TRADE results and plot

library(DESeq2)
library(sva)
library(ashr)
library(ggplot2)
library(mashr)
library(reshape2)
library(jaffelab)
library(gtools)

theme_bhr <- function(){
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = "None",
        legend.direction = "horizontal",
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(fill = "white"))
}

theme_bhr_legend <- function(){
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15),
        legend.title = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(fill = "white"))
}

compiled_schema_output_df <- read.csv("path/to/TRADE/allgene_TRADE1000perm_results.csv")
compiled_schema_output_df$p = compiled_schema_output_df$TI_pval
compiled_schema_output_df$exp = paste0(compiled_schema_output_df$gene, "_", compiled_schema_output_df$age)


#### Meta analyzing across regions ####

meta_TI_df <- data.frame(gene = c("Akap11",
                                  "Dagla",
                                  "Gria3",
                                  "Grin2a",
                                  "Sp4",
                                  "Srrm2",
                                  "Zmym2",
                                  "Grin2b"))

get_metabeta <- function(betas,ses) {
  we = 1/(ses^2)
  return(sum(betas * we)/sum(we))
}

get_metase <- function(ses) {
  we = 1/(ses^2)
  return(sqrt(1/sum(we)))
}

meta_TI_df$fisher.pval <- sapply(1:nrow(meta_TI_df),
                                 function(row) {
                                   gene = meta_TI_df$gene[row]
                                   
                                   gene_output = compiled_schema_output_df[compiled_schema_output_df$gene == gene,]
                                   gene_output$p[gene_output$p == 0] <- 0.001
                                   
                                   fishers_pval = poolr::fisher(gene_output$p)
                                   fishers_pval$p
                                 })

meta_TI_df$meta_TI <- sapply(1:nrow(meta_TI_df),
                             function(row) {
                               gene = meta_TI_df$gene[row]
                               
                               gene_output = compiled_schema_output_df[compiled_schema_output_df$gene == gene,]
                               gene_output$p[gene_output$p == 0] <- 0.001
                               gene_output$pseudo_se = gene_output$TI/qnorm(gene_output$p, lower.tail = FALSE)
                               
                               meta_TI = get_metabeta(gene_output$TI,gene_output$pseudo_se)
                             })

meta_TI_df$meta_TI_SE <- sapply(1:nrow(meta_TI_df),
                                function(row) {
                                  gene = meta_TI_df$gene[row]
                                  
                                  gene_output = compiled_schema_output_df[compiled_schema_output_df$gene == gene,]
                                  gene_output$p[gene_output$p == 0] <- 0.001
                                  gene_output$pseudo_se = gene_output$TI/qnorm(gene_output$p, lower.tail = FALSE)
                                  
                                  meta_TI = get_metase(gene_output$pseudo_se)
                                })


meta_TI_df$FDR = p.adjust(meta_TI_df$fisher.pval)
meta_TI_df$sig = stars.pval(meta_TI_df$FDR)
ggplot(data = meta_TI_df,
       mapping = aes(x = gene,
                     y = meta_TI,
                     ymin = meta_TI - meta_TI_SE,
                     ymax = meta_TI + meta_TI_SE,
                     color = gene,
                     #label = signif(FDR, 2)))+
                     label = sig )) +
  
  geom_pointrange()+
  geom_text(vjust = -2, col = "black", size = 6) +
  theme_bhr()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Gene", y = "Transcriptome-Wide Impact\nMeta-Analyzed Across Regions") 



#### Meta analyzing across mutants ####
compiled_schema_output_df <- read.csv("path/to/TRADE/allSCHEMA_TRADE1000perm_results.csv") #does not include Grin2b
compiled_schema_output_df$p = compiled_schema_output_df$TI_pval
compiled_schema_output_df$exp = paste0(compiled_schema_output_df$gene, "_", compiled_schema_output_df$age)


meta_TI_df <- data.frame(region = c("PFC", "SSC", "HP", "TH", "STR", "SN"))

get_metabeta <- function(betas,ses) {
  we = 1/(ses^2)
  return(sum(betas * we)/sum(we))
}

get_metase <- function(ses) {
  we = 1/(ses^2)
  return(sqrt(1/sum(we)))
}

meta_TI_df$fisher.pval <- sapply(1:nrow(meta_TI_df),
                                 function(row) {
                                   region = meta_TI_df$region[row]
                                   
                                   region_output = compiled_schema_output_df[compiled_schema_output_df$region == region,]
                                   region_output$p[region_output$p == 0] <- 0.001
                                   
                                   fishers_pval = poolr::fisher(region_output$p)
                                   fishers_pval$p
                                 })

meta_TI_df$meta_TI <- sapply(1:nrow(meta_TI_df),
                             function(row) {
                               region = meta_TI_df$region[row]
                               
                               region_output = compiled_schema_output_df[compiled_schema_output_df$region == region,]
                               region_output$p[region_output$p == 0] <- 0.001
                               region_output$pseudo_se = region_output$TI/qnorm(region_output$p, lower.tail = FALSE)
                               
                               meta_TI = get_metabeta(region_output$TI,region_output$pseudo_se)
                             })

meta_TI_df$meta_TI_SE <- sapply(1:nrow(meta_TI_df),
                                function(row) {
                                  region = meta_TI_df$region[row]
                                  
                                  region_output = compiled_schema_output_df[compiled_schema_output_df$region == region,]
                                  region_output$p[region_output$p == 0] <- 0.001
                                  region_output$pseudo_se = region_output$TI/qnorm(region_output$p, lower.tail = FALSE)
                                  
                                  meta_TI = get_metase(region_output$pseudo_se)
                                })

meta_TI_df$FDR = p.adjust(meta_TI_df$fisher.pval)
meta_TI_df$sig = stars.pval(meta_TI_df$FDR)

ggplot(data = meta_TI_df,
       mapping = aes(x = region,
                     y = meta_TI,
                     ymin = meta_TI - meta_TI_SE,
                     ymax = meta_TI + meta_TI_SE,
                     color = region,
                     #label = signif(fisher.pval, 2)))+
                     label = sig)) +
  
  geom_pointrange(size = 1)+
  geom_text(vjust = -2, col = "black", size = 6) +
  scale_color_manual(values = c("darkgoldenrod4", "yellow2", "navyblue", "seagreen3", "firebrick1", "deepskyblue3")) +
  theme_bhr()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "Meta-Transcriptome-Wide Impact") +
  theme(axis.text.x = element_text(size = 18, angle = 0), axis.text.y = element_text(size = 18), axis.title=element_text(size=18,face="bold"))

