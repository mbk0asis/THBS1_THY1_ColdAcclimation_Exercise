# GEOquery
# GSE156249
# Microarray data for mouse skeletal muscles before and after cold acclimation


library(GEOquery)
library(limma)
library(umap)
library(tidyverse)
library(ggrepel)
library(matrixStats)
library(ggpubr)


setwd("/home/sc/00--NGS/cold_acclimation_DrYang")

# load series and platform data from GEO
gset <- getGEO("GSE156249", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL11532", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# phenotypes
pheno <- pData(gset)
phenoData <- pheno %>% 
  mutate(group = case_when(str_detect(title,"before training") ~ "preEx", str_detect(title,"after training") ~ "postEx",
                           str_detect(title,"before cold") ~ "preCold", str_detect(title,"after 10 day") ~ "postCold")) %>%
  dplyr::select(group, subject=`subject id:ch1`)

# log2 transform
exp <- exprs(gset)
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
  exp[which(exp <= 0)] <- NaN
  exp <- log2(exp)}


# probe ID to gene symbols
p2s <- read.delim("probeIDs2Symbols", header = FALSE, row.names = 1)
colnames(p2s) <- c("refSeq","symbol")
head(p2s)
# secretory proteins
secProt <- read.delim("secreted_protein_Predicted_TheHumanProteinAtlas.tsv", header = TRUE)
head(secProt[,1:3])


## Exercise data =====================================================================================================================
exSamples <- phenoData %>% dplyr::filter(group %in% c("preEx","postEx")) 
exp_ex <- exp[, grep(paste(rownames(exSamples),collapse="|"), colnames(exp))]
# PCA
ntop=500
rv <- rowVars(exp_ex)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(exp_ex[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pcaDta <- pca$x[,1:2]
p1 <- merge(pcaDta, exSamples, by=0) %>%
  ggplot(aes(x=PC1, y=PC2, color=group)) + theme_classic() +
  geom_point(size=5) + coord_fixed() +
  scale_color_manual(values = c("magenta2","dodgerblue")) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100, digits = 2), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100, digits = 2), "% variance")) +
  ggtitle("Exercise")

# Differential expression analysis 
design <- exSamples %>% mutate(preEx = case_when(group %in% "preEx" ~ 1, TRUE ~ 0), postEx = case_when(group %in% "postEx" ~ 1, TRUE ~ 0)) %>%
  dplyr::select(preEx, postEx)

pval=0.05           
log2FoldChange=log2(1)

c <- makeContrasts("postEx-preEx", levels=design) 
fit <- lmFit(exp_ex, design) 
fit2<- contrasts.fit(fit,c) 
fit2<- eBayes(fit2) 
top <- topTable(fit2, coef=1, n=nrow(exp_ex), adjust = "fdr")
res1 <- merge(top, p2s, by=0) %>% mutate(group=case_when(logFC > log2FoldChange & P.Value < pval ~ "up", logFC < -log2FoldChange & P.Value < pval ~ "down", TRUE ~ "ns")) %>%
  mutate(thbs1 = case_when(symbol %in% "THBS1" ~ "THBS1", TRUE ~ ""))
nSamples <- res1 %>% dplyr::count(group)
p2 <- ggplot(res1, aes(x=logFC, y=-log10(P.Value), color=group, size=group)) + theme_classic() +
  geom_point(alpha=0.5) + geom_hline(yintercept = -log10(pval), linetype="dashed") + geom_vline(xintercept = 0) +
  geom_text(aes(label=thbs1), color="black", size=7) +
  scale_color_manual(values = c("dodgerblue","grey","magenta2")) +
  scale_size_manual(values = c(2,1,2)) + 
  ggtitle(paste("preEx (",nSamples$n[1],")  vs  postEx (", nSamples$n[3],")")) +
  coord_cartesian(xlim=c(-1.5,2.5),ylim=c(0,6))


## Cold acclimation data ======================================================================================================================
coldSamples <- phenoData %>% dplyr::filter(group %in% c("preCold","postCold")) 
exp_cold <- exp[, grep(paste(rownames(coldSamples),collapse="|"), colnames(exp))]

# PCA
ntop=500
rv <- rowVars(exp_cold)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(exp_cold[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pcaDta <- pca$x[,1:2]
p3 <- merge(pcaDta, coldSamples, by=0) %>%
  ggplot(aes(x=PC1, y=PC2, color=group)) + theme_classic() +
  geom_point(size=5) + coord_fixed() +
  scale_color_manual(values = c("magenta2","dodgerblue")) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100, digits = 2), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100, digits = 2), "% variance")) +
  ggtitle("Cold acclimation")
p3

# Differential expression analysis 
pval=0.05           
log2FoldChange=log2(1)
design <- coldSamples %>% mutate(preCold = case_when(group %in% "preCold" ~ 1, TRUE ~ 0), postCold = case_when(group %in% "postCold" ~ 1, TRUE ~ 0)) %>%
  dplyr::select(preCold, postCold)
c <- makeContrasts("postCold-preCold", levels=design) 
fit <- lmFit(exp_cold, design) 
fit2<- contrasts.fit(fit,c) 
fit2<- eBayes(fit2) 
top <- topTable(fit2, coef=1, n=nrow(exp_cold), adjust = "fdr")
res2 <- merge(top, p2s, by=0) %>% 
  mutate(thbs1 = case_when(symbol %in% "THBS1" ~ "THBS1", TRUE ~ "")) %>% 
  mutate(group=case_when(logFC > log2FoldChange & P.Value < pval~ "up",logFC < -log2FoldChange & P.Value < pval ~ "down", TRUE ~ "ns"))
write.csv(res2, "Raw_data_Fig.C_microarray_coldAccilimation.csv")

nSamples <- res2 %>% dplyr::count(group)
point2 <- res2 %>% dplyr::filter(thbs1 %in% "THBS1")
p4 <- ggplot(res2, aes(x=logFC, y=-log10(P.Value), color=group)) + theme_classic() +
  geom_point(alpha=0.7) +
  geom_point(data=point2, aes(x=logFC, y=-log10(P.Value)), color="red", size=3) +
  geom_hline(yintercept = -log10(pval), linetype="dashed") + 
  geom_vline(xintercept = 0) + geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linetype="dashed") +
  geom_text_repel(data=point2, aes(label=thbs1), color="black", size=8) +
  scale_color_manual(values = c("dodgerblue","grey","pink")) +
  scale_size_manual(values = c(3,1)) + 
  ggtitle(paste("preCold (",nSamples$n[1],")  vs  postCold (", nSamples$n[3],")")) +
  coord_cartesian(xlim=c(-1.5,2.5),ylim=c(0,6)) + xlab("log2FC") + ylab("-log10(pvalue)") +
  theme(axis.title = element_text(size=15),
        axis.text.x = element_text(size=11),axis.text.y = element_text(size=11))
p4

#
library(patchwork)
wrap_plots(p1, p3)   # PCA plots
wrap_plots(p2, p4, guides = 'collect')   # volcano plots

