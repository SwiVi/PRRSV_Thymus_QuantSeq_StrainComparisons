lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(AUCell)
library(writexl)
library(readxl)
library(tidyverse)
library(gplots)
library(ggplot2)
library(scales)
library(report)

## Make Gu et al. gene sets ----
# Gene set enrichment analysis - Gu ----

## Make gene sets ----
ref <- read_xlsx('/project/nadc_prrsv/Wiarda/ReferenceData/Gu2022_scRNAseqPigThymus/Thy_5999_FindAllMarkers.xlsx')
ref <- subset(ref, p_val_adj < 0.05 & avg_log2FC > 0.25)

c0_DN_C <- subset(ref, cluster == '0_DN_C')
c1_DN_Q <- subset(ref, cluster == '1_DN_Q')
c2_DP_C_1 <- subset(ref, cluster == '2_DP_C_1')
c3_DP_C_2 <- subset(ref, cluster == '3_DP_C_2')
c4_DP_Q <- subset(ref, cluster == '4_DP_Q')
c5_CD2pos_gd_T <- subset(ref, cluster == '5_CD2pos_gd_T')
c6_CD2neg_gd_T <- subset(ref, cluster == '6_CD2neg_gd_T')
c7_T_entry <- subset(ref, cluster == '7_T_entry')
c8_Treg1 <- subset(ref, cluster == '8_Treg1')
c9_Treg2 <- subset(ref, cluster == '9_Treg2')
c10_CD8SP <- subset(ref, cluster == '10_CD8SP')
c11_CD4SP <- subset(ref, cluster == '11_CD4SP')
c12_ISG_CD8_T <- subset(ref, cluster == '12_ISG_CD8_T')
c13_cytotoxic_CD8_T <- subset(ref, cluster == '13_cytotoxic_CD8_T')
c14_CD8aa <- subset(ref, cluster == '14_CD8aa')
c15_B_cells <- subset(ref, cluster == '15_B_cells')

geneSets <- list(c0_DN_C = c0_DN_C$gene,
                 c1_DN_Q = c1_DN_Q$gene,
                 c2_DP_C_1 = c2_DP_C_1$gene,
                 c3_DP_C_2 = c3_DP_C_2$gene,
                 c4_DP_Q = c4_DP_Q$gene,
                 c5_CD2pos_gd_T = c5_CD2pos_gd_T$gene,
                 c6_CD2neg_gd_T = c6_CD2neg_gd_T$gene,
                 c7_T_entry = c7_T_entry$gene,
                 c8_Treg1 = c8_Treg1$gene,
                 c9_Treg2 = c9_Treg2$gene,
                 c10_CD8SP = c10_CD8SP$gene,
                 c11_CD4SP = c11_CD4SP$gene,
                 c12_ISG_CD8_T = c12_ISG_CD8_T$gene,
                 c13_cytotoxic_CD8_T = c13_cytotoxic_CD8_T$gene,
                 c14_CD8aa = c14_CD8aa$gene,
                 c15_B_cells = c15_B_cells$gene)

# Read in refined gene count table and target file ----
counts <- read.csv('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_GeneCounts_Refined.csv', check.names = FALSE, row.names = 1)
target <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_TargetFile_Refined.xlsx')
identical(colnames(counts), target$name) # make sure the sample IDs match (identical = TRUE)

# Filter genes ----
keep<- rowSums(counts)>=10
counts <- counts[keep,]

## Calculate AUC scores ----
exprMatrix <- as.matrix(counts)

set.seed(123)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # calculate cell signature AUC score for each gene set in each cell

AUC <- as.data.frame(getAUC(cells_AUC))
AUC <- data.frame(t(AUC))
AUC$sample <- rownames(AUC)
identical(AUC$sample, target$name) # make sure TRUE
AUC$group <- target$group
write_xlsx(AUC, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/GSEA_GuThymusGeneSets_AUCScores.xlsx')

## Make heatmaps (all samples hierarchical) ----
my_palette <- colorRampPalette(c("navy", "grey95", "red2"))(n = 399)
colors = seq(-4, 4,length=400)
heatmap.2(as.matrix(t(AUC[1:16])), Rowv=T, Colv=T, dendrogram= "both", margins = c(12,12), breaks = colors,
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # swap out cols for different treatment variables

cols <- c(rep("#F8766D", 14), rep("#7CAE00", 15), rep("#00BFC4", 15), rep("#C77CFF", 15))
heatmap.2(as.matrix(t(AUC[1:16])), Rowv=T, dendrogram= "both", margins = c(12,12), breaks = colors,
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette, ColSideColors = cols) # swap out cols for different treatment variables
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Mock", "MN184", "L1C.5", "Jxwn06"), # category labels
       col = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),  # color key
       lty= 1, lwd = 10) 

cols <- c(rep("grey80", 5), rep("grey40", 5), rep("black", 4), rep("grey80", 5), rep("grey40", 5), rep("black", 5), 
          rep("grey80", 5), rep("grey40", 5), rep("black", 5), rep("grey80", 5), rep("grey40", 5), rep("black", 5))
heatmap.2(as.matrix(t(AUC[1:16])), Rowv=T, dendrogram= "both", margins = c(12,12), breaks = colors,
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette, ColSideColors = cols) # swap out cols for different treatment variables
legend("topright",      # location of the legend on the heatmap plot
       legend = c("2dpi", "6dpi", "10dpi"), # category labels
       col = c("grey80", "grey40", "black"),  # color key
       lty= 1, lwd = 10) 

## Plot individual gene sets & run stats ----

AUC <- AUC %>% 
  pivot_longer(c(colnames(AUC[, 1:(ncol(AUC) - 2)])), names_to = "geneSet", values_to = "AUC")
AUC$group <- factor(AUC$group, levels = c('Mock_2dpi', 'MN184_2dpi','L1C5_2dpi', 'Jxwn06_2dpi', 
                                          'Mock_6dpi', 'MN184_6dpi', 'L1C5_6dpi', 'Jxwn06_6dpi', 
                                          'Mock_10dpi','MN184_10dpi','L1C5_10dpi','Jxwn06_10dpi'))
AUC$dpi <- sub('.*\\_', '', AUC$group)
AUC$DaysPostInfection <- gsub('.{3}$', '', AUC$dpi)
AUC$DaysPostInfection <- as.numeric(AUC$DaysPostInfection)
AUC$Treatment <- sub("_[^_]+$", "", AUC$group)

### c0_DN_C ----
sub <- subset(AUC, geneSet == 'c0_DN_C')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c0_DN_C") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c1_DN_Q ----
sub <- subset(AUC, geneSet == 'c1_DN_Q')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c1_DN_Q") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c2_DP_C_1 ----
sub <- subset(AUC, geneSet == 'c2_DP_C_1')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c2_DP_C_1") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c3_DP_C_2 ----
sub <- subset(AUC, geneSet == 'c3_DP_C_2')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c3_DP_C_2") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c4_DP_Q ----
sub <- subset(AUC, geneSet == 'c4_DP_Q')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c4_DP_Q") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c5_CD2pos_gd_T ----
sub <- subset(AUC, geneSet == 'c5_CD2pos_gd_T')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c5_CD2pos_gd_T") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c6_CD2neg_gd_T ----
sub <- subset(AUC, geneSet == 'c6_CD2neg_gd_T')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c6_CD2neg_gd_T") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c7_T_entry ----
sub <- subset(AUC, geneSet == 'c7_T_entry')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c7_T_entry") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c8_Treg1 ----
sub <- subset(AUC, geneSet == 'c8_Treg1')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c8_Treg1") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c9_Treg2 ----
sub <- subset(AUC, geneSet == 'c9_Treg2')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c9_Treg2") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c10_CD8SP ----
sub <- subset(AUC, geneSet == 'c10_CD8SP')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c10_CD8SP") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c11_CD4SP ----
sub <- subset(AUC, geneSet == 'c11_CD4SP')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c11_CD4SP") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c12_ISG_CD8_T ----
sub <- subset(AUC, geneSet == 'c12_ISG_CD8_T')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c12_ISG_CD8_T") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c13_cytotoxic_CD8_T ----
sub <- subset(AUC, geneSet == 'c13_cytotoxic_CD8_T')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c13_cytotoxic_CD8_T") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c14_CD8aa ----
sub <- subset(AUC, geneSet == 'c14_CD8aa')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c14_CD8aa") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

### c15_B_cells ----
sub <- subset(AUC, geneSet == 'c15_B_cells')
ggplot(sub, aes(x=DaysPostInfection, y = AUC, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(sub[(sub$DaysPostInfection == 2 | sub$DaysPostInfection == 6 | sub$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, AUC, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  ggtitle("c15_B_cells") +
  xlab('Days Post-Infection') +
  ylab('Area under the curve (AUC) score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

df <- split(sub, sub$DaysPostInfection)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

df <- split(sub, sub$Treatment)
wilcox <- lapply(df, function(g) pairwise.wilcox.test(g$AUC, g$group, p.adjust.method="fdr"))
wilcox

report(sessionInfo())
#Analyses were conducted using the R Statistical language (version 4.3.2; R Core Team, 2023) on Ubuntu 22.04.3 LTS, using the packages AUCell (version 1.24.0; Aibar S et al., 2017), lubridate (version 1.9.3; Grolemund G, Wickham H, 2011), report (version 0.5.9; Makowski D et
#                                                                                                                                                                                                                                                      al., 2023), tibble (version 3.2.1; Müller K, Wickham H, 2023), writexl (version 1.5.0; Ooms J, 2024), gplots (version 3.1.3.1; Warnes G et al., 2024), ggplot2 (version 3.5.1; Wickham H, 2016), forcats (version 1.0.0; Wickham H, 2023), stringr (version 1.5.1; Wickham H, 2023),
#tidyverse (version 2.0.0; Wickham H et al., 2019), readxl (version 1.4.3; Wickham H, Bryan J, 2023), dplyr (version 1.1.4; Wickham H et al., 2023), purrr (version 1.0.2; Wickham H, Henry L, 2023), readr (version 2.1.4; Wickham H et al., 2023), scales (version 1.3.0; Wickham H
#                                                                                                                                                                                                                                                            et al., 2023) and tidyr (version 1.3.0; Wickham H et al., 2023).
#
#References
#----------
#  - Aibar S, Bravo Gonzalez-Blas C, Moerman T, Huynh-Thu V, Imrichova H, Hulselmans G, Rambow F, Marine J, Geurts P, Aerts J, van den Oord J, Kalender Atak Z, Wouters J, Aerts S (2017). “SCENIC: Single-Cell Regulatory Network Inference And Clustering.” _Nature Methods_, *14*,
#1083-1086. doi:10.1038/nmeth.4463 <https://doi.org/10.1038/nmeth.4463>. Aibar S, Aerts S (2016). “AUCell.”
#- Grolemund G, Wickham H (2011). “Dates and Times Made Easy with lubridate.” _Journal of Statistical Software_, *40*(3), 1-25. <https://www.jstatsoft.org/v40/i03/>.
#- Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” _CRAN_. <https://easystats.github.io/report/>.
#- Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package version 3.2.1, https://github.com/tidyverse/tibble, <https://tibble.tidyverse.org/>.
#- Ooms J (2024). _writexl: Export Data Frames to Excel 'xlsx' Format_. R package version 1.5.0https://docs.ropensci.org/writexl/ (website) https://github.com/ropensci/writexl (devel) https://libxlsxwriter.github.io (upstream),
#<https://docs.ropensci.org/writexl/%20(website)https://github.com/ropensci/writexl%20(devel)https://libxlsxwriter.github.io%20(upstream)>.
#- R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#- Warnes G, Bolker B, Bonebakker L, Gentleman R, Huber W, Liaw A, Lumley T, Maechler M, Magnusson A, Moeller S, Schwartz M, Venables B (2024). _gplots: Various R Programming Tools for Plotting Data_. R package version 3.1.3.1, <https://CRAN.R-project.org/package=gplots>.
#- Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
#- Wickham H (2023). _forcats: Tools for Working with Categorical Variables (Factors)_. R package version 1.0.0, https://github.com/tidyverse/forcats, <https://forcats.tidyverse.org/>.
#- Wickham H (2023). _stringr: Simple, Consistent Wrappers for Common String Operations_. R package version 1.5.1, https://github.com/tidyverse/stringr, <https://stringr.tidyverse.org>.
#- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the
#tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
#- Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package version 1.4.3, https://github.com/tidyverse/readxl, <https://readxl.tidyverse.org>.
#- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.4, https://github.com/tidyverse/dplyr, <https://dplyr.tidyverse.org>.
#- Wickham H, Henry L (2023). _purrr: Functional Programming Tools_. R package version 1.0.2, https://github.com/tidyverse/purrr, <https://purrr.tidyverse.org/>.
#- Wickham H, Hester J, Bryan J (2023). _readr: Read Rectangular Text Data_. R package version 2.1.4, https://github.com/tidyverse/readr, <https://readr.tidyverse.org>.
#- Wickham H, Pedersen T, Seidel D (2023). _scales: Scale Functions for Visualization_. R package version 1.3.0, https://github.com/r-lib/scales, <https://scales.r-lib.org>.
#- Wickham H, Vaughan D, Girlich M (2023). _tidyr: Tidy Messy Data_. R package version 1.3.0, https://github.com/tidyverse/tidyr, <https://tidyr.tidyverse.org>.