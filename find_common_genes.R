# To find common genes 
#

library(readxl)
library(VennDiagram)
library(tidyverse)
library(gplots)
library(viridis)
library(magrittr)

setwd("/home/sc/00--NGS/cold_acclimation_DrYang/ColdAcclimation")

## Bicycle exercise data ===========================================================================================================
# Read an excel file (Bicycle exercise)
excel_file <- "sigDEGs_Bicycle_Exercise_secretory_proteins.xlsx"
sheets <- excel_sheets(excel_file)
sheet_list <- lapply(sheets, function(sheet) read_excel(excel_file, sheet = sheet))
names(sheet_list) <- c("A1_A2","A1_A3","B1_B2","B1_B3","proteomics")
common_genes <- Reduce(intersect, list(sheet_list[[1]]$Symbol, sheet_list[[2]]$Symbol, sheet_list[[3]]$Symbol, sheet_list[[4]]$Symbol))

list <- list()
for (n in 1:4) {
  list[[n]] <- sheet_list[[n]] %>% dplyr::filter(Symbol %in% common_genes) %>% arrange(Symbol)
}

list[[1]] %<>% mutate(group="A2A1")
list[[2]] %<>% mutate(group="A3A1")
list[[3]] %<>% mutate(group="B2B1")
list[[4]] %<>% mutate(group="B3B1")

exp_common_genes <- list[[1]] %>% inner_join(list[[2]], by="Symbol")%>% inner_join(list[[3]], by="Symbol")%>% inner_join(list[[4]], by="Symbol")
write.csv(exp_common_genes, "Raw_data_Fig.C_RNAseq_coldAcclimation.csv")

dta <- exp_common_genes %>% dplyr::select(1,contains("FC")) %>% column_to_rownames("Symbol")
colnames(dta) <- c("A2A1","A3A1", "B2B1","B3B1")
breaks <- unique(c(seq(-4, 0, 0.1),seq(0, 4, 0.1)))
heatmap.2(as.matrix(dta), trace="none",density="none",
               margins = c(5,15), breaks = breaks, Colv = FALSE, Rowv = TRUE,
               col=viridis(length(breaks)-1, option = "D"))


dta2 <- do.call(rbind, list) %>% mutate(fc2 = case_when(FC > 5 ~ 5, TRUE ~ FC))
dta2$Symbol <- factor(dta2$Symbol, levels = rev(unique(dta2$Symbol)))
ggplot(dta2, aes(x=group, y=Symbol, color=fc2, size=-log10(qval))) + theme_classic() +
  geom_point() + scale_color_viridis() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=12,angle=90, vjust = 0.5),
        axis.text.y = element_text(size=12, color="black"))
unique(dta2$Symbol)

# Create a list of the sets
set_list <- list(Set1 = sheet_list[[1]]$Symbol, Set2 = sheet_list[[2]]$Symbol, Set3 = sheet_list[[3]]$Symbol, Set4 = sheet_list[[4]]$Symbol)
set_list

library(ggVennDiagram)
ggVennDiagram(set_list, label_alpha = 0, label_size = 7, label="count", set_size=0,
              category.names = c("A2/A1","A3/A1","B2/B1", "B3/B1")) +
  #scale_fill_viridis_c(option="C")
  ggplot2::scale_fill_gradient(low="dodgerblue",high = "violet")





## PROTEOMICS (blood) =======================================================================
dta3 <- sheet_list[[5]] %>% mutate(group=case_when(fold.Y > log2(1.5) ~ "up", fold.Y < -log2(1.5) ~ "down", TRUE ~ "nc")) %>% rename(symbol=`Gene names`)
thbs1 <- dta3 %>% dplyr::filter(symbol %in% "Thbs1")
table(dta3$group)
library(ggrepel)
ggplot(dta3, aes(x=log2(mean), y= fold.Y, color=group)) + theme_classic() + xlab("Mean expression") + ylab("log2FC") +
  geom_point() + scale_color_manual(values = c("dodgerblue","grey","pink")) +
  geom_hline(yintercept = c(log2(1.5),-log2(1.5)), linetype="dashed") +
  geom_hline(yintercept = 0) +
  geom_point(data=thbs1, aes(x=log2(mean), y= fold.Y), color="red", size=2.5) +
  geom_text_repel(data=thbs1, aes(x=log2(mean), y= fold.Y, label=symbol), color="black", size=6) +
  theme(axis.text = element_text(size=10))




## Overlapping genes in all three studies =============================================================================
exercise <- common_genes[-grep("MMP28|CYR61", common_genes)]
proteomics <- na.omit(toupper(dta3[dta3$fold.Y > 0, ]$symbol))
proteomics[1:517]
cold <- res2 %>% dplyr::filter(logFC > 0 & P.Value < 0.05 & !(symbol %in% "")) %>% pull(symbol)
list <- list(exercise,proteomics[1:517],cold)
list

library(ggVennDiagram)
ggVennDiagram(list, label_alpha = 0, label_size = 7, label="count", set_size=7,
              category.names = c("Ex","Prot","Cold")) +
  #ggplot2::scale_fill_gradient2(low="#B6D3EA",mid="#BFE7C5",high="#FCDDB7")
  ggplot2::scale_fill_gradient(low="dodgerblue",high = "pink")
  

overlap <- calculate.overlap(list)
overlap2 <- data.frame(unlist(overlap)) %>% rownames_to_column("area") 
overlap2$area2 <- substr(overlap2$area, 1, 2) 
write.csv(overlap2, "overlapping_genes.csv")
