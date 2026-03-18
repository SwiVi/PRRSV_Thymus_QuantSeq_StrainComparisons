lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(readxl)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(gplots)
library(EnhancedVolcano)
library(writexl)
library(report)

# Read in refined gene count table and target file ----
counts <- read.csv('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_GeneCounts_Refined.csv', check.names = FALSE, row.names = 1)
target <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_TargetFile_Refined.xlsx')
identical(colnames(counts), target$name) # make sure the sample IDs match (identical = TRUE)

# Define treatment factor ----
group <- factor(target$group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi',
                                         'MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi',
                                         'L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi',
                                         'Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))
group

# Create DESeq2 matrix: ----
dds<-DESeqDataSetFromMatrix(countData=counts, colData=target, design = ~group)
dds

# Filter genes ----
keep<- rowSums(counts(dds))>=10
dds<-dds[keep,]
dds

# Apply DESeq2 ----
d.dds<- DESeq(dds)
d.dds

# Create PCA ----
vsdb<- varianceStabilizingTransformation(d.dds)
vsdb$group <- factor(vsdb$group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi',
                                            'MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi',
                                            'L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi',
                                            'Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))

#plotPCA(vsdb, intgroup=c("group"), ntop = 2000)
#plotPCA(vsdb, intgroup=c("DaysPostInfection"))
#plotPCA(vsdb, intgroup=c("Treatment"))
#plotPCA(vsdb, intgroup=c("group"), ntop = 2000) + stat_ellipse(aes(fill = group), level = 0.95, geom = 'polygon', alpha = 0.1)
#plotPCA(vsdb, intgroup=c("group")) + stat_ellipse(aes(fill = group), level = 0.95)

d <-plotPCA(vsdb, intgroup=c("group"), ntop = 2000)
d <- d$data
cen <- d %>%
  group_by(group) %>%
  summarise(
    centroidPC1 = mean(PC1),
    centroidPC2 = mean(PC2))
d$centroidPC1 <- cen$centroidPC1[match(d$group, cen$group)]
d$centroidPC2 <- cen$centroidPC2[match(d$group, cen$group)]

#for(i in 1:12){
#  print(hue_pal()(i))
#}

## Plot all samples ----
ggplot(data = d,aes(x=PC1, y=PC2)) +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#FF64B0", "#F564E3", "#C77CFF", "#619CFF", "#00B4F0", "#00BFC4", 
                                                                 "#00C08B", "#00BA38", "#7CAE00", "#B79F00", "#DE8C00", "#F8766D")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi', 'MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi', 
                                                 'L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi', 'Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi')),
                 col = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi', 'MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi', 
                                                'L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi', 'Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#FF64B0", "#F564E3", "#C77CFF", "#619CFF", "#00B4F0", "#00BFC4", 
                                                                   "#00C08B", "#00BA38", "#7CAE00", "#B79F00", "#DE8C00", "#F8766D")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  labs(title="All samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot JXwn06 samples ----
sub <- subset(d, group == 'Jxwn06_2dpi' | group == 'Jxwn06_6dpi' | group == 'Jxwn06_10dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#B79F00", "#DE8C00", "#F8766D")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi')),
                 col = factor(group, levels = c('Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#B79F00", "#DE8C00", "#F8766D")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_path(aes(x=centroidPC1, y=centroidPC2), size = 2, col = 'grey40') +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('Jxwn06_2dpi', 'Jxwn06_6dpi', 'Jxwn06_10dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="Jxwn06 samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot L1C.5 samples ----
sub <- subset(d, group == 'L1C5_2dpi' | group == 'L1C5_6dpi' | group == 'L1C5_10dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#00C08B", "#00BA38", "#7CAE00")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi')),
                 col = factor(group, levels = c('L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#00C08B", "#00BA38", "#7CAE00")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_path(aes(x=centroidPC1, y=centroidPC2), size = 2, col = 'grey40') +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('L1C5_2dpi', 'L1C5_6dpi', 'L1C5_10dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="L1C5 samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot MN184 samples ----
sub <- subset(d, group == 'MN184_2dpi' | group == 'MN184_6dpi' | group == 'MN184_10dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#619CFF", "#00B4F0", "#00BFC4")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi')),
                 col = factor(group, levels = c('MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#619CFF", "#00B4F0", "#00BFC4")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_path(aes(x=centroidPC1, y=centroidPC2), size = 2, col = 'grey40') +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('MN184_2dpi', 'MN184_6dpi', 'MN184_10dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="MN184 samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot mock samples ----
sub <- subset(d, group == 'Mock_2dpi' | group == 'Mock_6dpi' | group == 'Mock_10dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#FF64B0", "#F564E3", "#C77CFF")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi')),
                 col = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#FF64B0", "#F564E3", "#C77CFF")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_path(aes(x=centroidPC1, y=centroidPC2), size = 2, col = 'grey40') +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('Mock_2dpi', 'Mock_6dpi', 'Mock_10dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="Mock samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot 2dpi samples ----
sub <- subset(d, group == 'Mock_2dpi' | group == 'MN184_2dpi' | group == 'L1C5_2dpi' | group == 'Jxwn06_2dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('Mock_2dpi', 'MN184_2dpi', 'L1C5_2dpi', 'Jxwn06_2dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('Mock_2dpi', 'MN184_2dpi', 'L1C5_2dpi', 'Jxwn06_2dpi')),
                 col = factor(group, levels = c('Mock_2dpi', 'MN184_2dpi', 'L1C5_2dpi', 'Jxwn06_2dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('Mock_2dpi', 'MN184_2dpi', 'L1C5_2dpi', 'Jxwn06_2dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('Mock_2dpi', 'MN184_2dpi', 'L1C5_2dpi', 'Jxwn06_2dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="2dpi samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot 6dpi samples ----
sub <- subset(d, group == 'Mock_6dpi' | group == 'MN184_6dpi' | group == 'L1C5_6dpi' | group == 'Jxwn06_6dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('Mock_6dpi', 'MN184_6dpi', 'L1C5_6dpi', 'Jxwn06_6dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('Mock_6dpi', 'MN184_6dpi', 'L1C5_6dpi', 'Jxwn06_6dpi')),
                 col = factor(group, levels = c('Mock_6dpi', 'MN184_6dpi', 'L1C5_6dpi', 'Jxwn06_6dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('Mock_6dpi', 'MN184_6dpi', 'L1C5_6dpi', 'Jxwn06_6dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('Mock_6dpi', 'MN184_6dpi', 'L1C5_6dpi', 'Jxwn06_6dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="6dpi samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

## Plot 10dpi samples ----
sub <- subset(d, group == 'Mock_10dpi' | group == 'MN184_10dpi' | group == 'L1C5_10dpi' | group == 'Jxwn06_10dpi')
ggplot(data = sub,aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(x=PC1, y=PC2, group=group, 
                   fill = factor(group, levels = c('Mock_10dpi', 'MN184_10dpi', 'L1C5_10dpi', 'Jxwn06_10dpi'))), 
               inherit.aes = FALSE, alpha=0.1, color = 'grey40', lwd = 0.1, geom = "polygon") +
  scale_fill_manual("Strain + Days Post-Inoculation", values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  geom_point(aes(x=PC1, y=PC2, group=group, 
                 fill = factor(group, levels = c('Mock_10dpi', 'MN184_10dpi', 'L1C5_10dpi', 'Jxwn06_10dpi')),
                 col = factor(group, levels = c('Mock_10dpi', 'MN184_10dpi', 'L1C5_10dpi', 'Jxwn06_10dpi'))),
             size=4, alpha = 0.8, shape = 19) +
  geom_segment(aes(xend=centroidPC1, yend=centroidPC2, col = factor(group, levels = c('Mock_10dpi', 'MN184_10dpi', 'L1C5_10dpi', 'Jxwn06_10dpi'))),
               linejoin = 'round', size=1, alpha = 0.8) + 
  scale_colour_manual("Strain + Days Post-Inoculation", values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  theme_bw() +
  scale_y_continuous(limits = c(-53, 43), breaks = seq(-50, 25, by = 25)) +
  scale_x_continuous(limits = c(-50, 100), breaks = seq(-50, 100, by = 25)) +
  geom_point(aes(x=centroidPC1, y=centroidPC2,
                 col = factor(group, levels = c('Mock_10dpi', 'MN184_10dpi', 'L1C5_10dpi', 'Jxwn06_10dpi'))),
             size = 8, alpha = 0.8) +
  labs(title="10dpi samples", x ="PC1 (45% variance)", y = "PC2 (8% variance)")

# Make heatmap ----
topVarGenes <- head(order(rowVars(assay(vsdb)),decreasing=TRUE),2000)

vsdb_table<- as.data.frame(assay(vsdb))
a<-vsdb_table[rownames(vsdb_table)[topVarGenes],]
my_palette <- colorRampPalette(c("navy", "grey95", "red2"))(n = 399)
colors = seq(-3, 3,length=400)
#cols <- c(rep("#B79F00", 5), rep("#DE8C00", 5), rep("#F8766D", 4), rep("#00C08B", 5), rep("#00BA38", 5), rep("#7CAE00", 5),
#          rep("#619CFF", 5), rep("#00B4F0", 5), rep("#00BFC4", 5), rep("#FF64B0", 5), rep("#F564E3", 5), rep("#C77CFF", 5))
cols <- c(rep("#F8766D", 14), rep("#7CAE00", 15), rep("#00BFC4", 15), rep("#C77CFF", 15))
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both", margins = c(12,12),
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette, ColSideColors = cols, breaks = colors) # swap out cols for different treatment variables
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Mock", "MN184", "LC1.5", "Jxwn06"), # category labels
       col = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),  # color key
       lty= 1, lwd = 10) 

cols <- c(rep("grey80", 5), rep("grey40", 5), rep("black", 4), rep("grey80", 5), rep("grey40", 5), rep("black", 5), 
          rep("grey80", 5), rep("grey40", 5), rep("black", 5), rep("grey80", 5), rep("grey40", 5), rep("black", 5))
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both", margins = c(12,12),
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette, ColSideColors = cols, breaks = colors) # swap out cols for different treatment variables
legend("topright",      # location of the legend on the heatmap plot
       legend = c("2dpi", "6dpi", "10dpi"), # category labels
       col = c("grey80", "grey40", "black"),  # color key
       lty= 1, lwd = 10) 

# DGE results ----

## L1C.5 vs mock ----
# DGE results

LCvMock_2<-results(d.dds, contrast=c("group","L1C5_2dpi", "Mock_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(LCvMock_2)

LCvMock_6<-results(d.dds, contrast=c("group","L1C5_6dpi", "Mock_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(LCvMock_6)

LCvMock_10<-results(d.dds, contrast=c("group","L1C5_10dpi", "Mock_10dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(LCvMock_10)

# volcano plots for each comp:
EnhancedVolcano(toptable = LCvMock_2,
                title = "Mock (left) v L1C.5 (right) at 2dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = LCvMock_6,
                title = "Mock (left) v L1C.5 (right) at 6dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = LCvMock_10,
                title = "Mock (left) v L1C.5 (right) at 10dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## JXwno6 vs mock ----
# DGE results

JXvMock_2<-results(d.dds, contrast=c("group","Jxwn06_2dpi", "Mock_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JXvMock_2)

JXvMock_6<-results(d.dds, contrast=c("group","Jxwn06_6dpi", "Mock_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JXvMock_6)

JXvMock_10<-results(d.dds, contrast=c("group","Jxwn06_10dpi", "Mock_10dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JXvMock_10)

# volcano plots for each comp:
EnhancedVolcano(toptable = JXvMock_2,
                title = "Mock (left) v Jxwn06 (right) at 2dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = JXvMock_6,
                title = "Mock (left) v Jxwn06 (right) at 6dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = JXvMock_10,
                title = "Mock (left) v Jxwn06 (right) at 10dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## MN184 vs mock ----
# DGE results

MNvMock_2<-results(d.dds, contrast=c("group","MN184_2dpi", "Mock_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(MNvMock_2)

MNvMock_6<-results(d.dds, contrast=c("group","MN184_6dpi", "Mock_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(MNvMock_6)

MNvMock_10<-results(d.dds, contrast=c("group","MN184_10dpi", "Mock_10dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(MNvMock_10)

# volcano plots for each comp:
EnhancedVolcano(toptable = MNvMock_2,
                title = "Mock (left) v MN184 (right) at 2dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = MNvMock_6,
                title = "Mock (left) v MN184 (right) at 6dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = MNvMock_10,
                title = "Mock (left) v MN184 (right) at 10dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## L1C.5 vs MN184 ----
# DGE results

LCvMN_2<-results(d.dds, contrast=c("group","L1C5_2dpi", "MN184_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(LCvMN_2)

LCvMN_6<-results(d.dds, contrast=c("group","L1C5_6dpi", "MN184_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(LCvMN_6)

LCvMN_10<-results(d.dds, contrast=c("group","L1C5_10dpi", "MN184_10dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(LCvMN_10)

# volcano plots for each comp:
EnhancedVolcano(toptable = LCvMN_2,
                title = "MN184 (left) v L1C5 (right) at 2dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = LCvMN_6,
                title = "MN184 (left) v L1C5 (right) at 6dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = LCvMN_10,
                title = "MN184 (left) v L1C5 (right) at 10dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## JXwn06 vs MN184 ----
# DGE results

JxvMN_2<-results(d.dds, contrast=c("group","Jxwn06_2dpi", "MN184_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JxvMN_2)

JxvMN_6<-results(d.dds, contrast=c("group","Jxwn06_6dpi", "MN184_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JxvMN_6)

JxvMN_10<-results(d.dds, contrast=c("group","Jxwn06_10dpi", "MN184_10dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JxvMN_10)

# volcano plots for each comp:
EnhancedVolcano(toptable = JxvMN_2,
                title = "MN184 (left) v Jxwn06 (right) at 2dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = JxvMN_6,
                title = "MN184 (left) v Jxwn06 (right) at 6dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = JxvMN_10,
                title = "MN184 (left) v Jxwn06 (right) at 10dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## JXwn06 vs L1C.5 ----
# DGE results

JxvLC_2<-results(d.dds, contrast=c("group","Jxwn06_2dpi", "L1C5_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JxvLC_2)

JxvLC_6<-results(d.dds, contrast=c("group","Jxwn06_6dpi", "L1C5_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JxvLC_6)

JxvLC_10<-results(d.dds, contrast=c("group","Jxwn06_10dpi", "L1C5_10dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(JxvLC_10)

# volcano plots for each comp:
EnhancedVolcano(toptable = JxvLC_2,
                title = "L1C5 (left) v Jxwn06 (right) at 2dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = JxvLC_6,
                title = "L1C5 (left) v Jxwn06 (right) at 6dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = JxvLC_10,
                title = "L1C5 (left) v Jxwn06 (right) at 10dpi",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## Mock timepoints ----
# DGE results

Mock2v6<-results(d.dds, contrast=c("group","Mock_6dpi", "Mock_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(Mock2v6)

Mock2v10<-results(d.dds, contrast=c("group","Mock_10dpi", "Mock_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(Mock2v10)

Mock6v10<-results(d.dds, contrast=c("group","Mock_10dpi", "Mock_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(Mock6v10)

# volcano plots for each comp:
EnhancedVolcano(toptable = Mock2v6,
                title = "Mock 2dpi (left) vs 6dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = Mock2v10,
                title = "Mock 2dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = Mock6v10,
                title = "Mock 6dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## MN184 timepoints ----
# DGE results

MN1842v6<-results(d.dds, contrast=c("group","MN184_6dpi", "MN184_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(MN1842v6)

MN1842v10<-results(d.dds, contrast=c("group","MN184_10dpi", "MN184_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(MN1842v10)

MN1846v10<-results(d.dds, contrast=c("group","MN184_10dpi", "MN184_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(MN1846v10)

# volcano plots for each comp:
EnhancedVolcano(toptable = MN1842v6,
                title = "MN184 2dpi (left) vs 6dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = MN1842v10,
                title = "MN184 2dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = MN1846v10,
                title = "MN184 6dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## L1C.5 timepoints ----
# DGE results

L1C52v6<-results(d.dds, contrast=c("group","L1C5_6dpi", "L1C5_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(L1C52v6)

L1C52v10<-results(d.dds, contrast=c("group","L1C5_10dpi", "L1C5_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(L1C52v10)

L1C56v10<-results(d.dds, contrast=c("group","L1C5_10dpi", "L1C5_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(L1C56v10)

# volcano plots for each comp:
EnhancedVolcano(toptable = L1C52v6,
                title = "L1C5 2dpi (left) vs 6dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = L1C52v10,
                title = "L1C5 2dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = L1C56v10,
                title = "L1C5 6dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## JXwn06 timepoints ----
# DGE results

Jxwn062v6<-results(d.dds, contrast=c("group","Jxwn06_6dpi", "Jxwn06_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(Jxwn062v6)

Jxwn062v10<-results(d.dds, contrast=c("group","Jxwn06_10dpi", "Jxwn06_2dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(Jxwn062v10)

Jxwn066v10<-results(d.dds, contrast=c("group","Jxwn06_10dpi", "Jxwn06_6dpi"),alpha=0.05, lfcThreshold = 0, pAdjustMethod = 'fdr')
summary(Jxwn066v10)

# volcano plots for each comp:
EnhancedVolcano(toptable = Jxwn062v6,
                title = "Jxwn06 2dpi (left) vs 6dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = Jxwn062v10,
                title = "Jxwn06 2dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

EnhancedVolcano(toptable = Jxwn066v10,
                title = "Jxwn06 6dpi (left) vs 10dpi (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 0.5, pointSize = 1,
                ylim = c(0, 11),
                xlim = c(-7.5,7.5)) #+ coord_flip()

## Save DE results ----
dir.create('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/')
ann <- read_excel('/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # read in annotation key that has Ensembl IDs

LCvMock_2$Gene <- rownames(LCvMock_2)
LCvMock_2 <- data.frame(LCvMock_2)
LCvMock_2$EnsemblID <- ann$ENSID[match(LCvMock_2$Gene, ann$FinalList)]

LCvMock_6$Gene <- rownames(LCvMock_6)
LCvMock_6 <- data.frame(LCvMock_6)
LCvMock_6$EnsemblID <- ann$ENSID[match(LCvMock_6$Gene, ann$FinalList)]

LCvMock_10$Gene <- rownames(LCvMock_10)
LCvMock_10 <- data.frame(LCvMock_10)
LCvMock_10$EnsemblID <- ann$ENSID[match(LCvMock_10$Gene, ann$FinalList)]

write_xlsx(data.frame(LCvMock_2), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_LCvsMock_2dpi.xlsx')
write_xlsx(data.frame(LCvMock_6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_LCvsMock_6dpi.xlsx')
write_xlsx(data.frame(LCvMock_10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_LCvsMock_10dpi.xlsx')

JXvMock_2$Gene <- rownames(JXvMock_2)
JXvMock_2 <- data.frame(JXvMock_2)
JXvMock_2$EnsemblID <- ann$ENSID[match(JXvMock_2$Gene, ann$FinalList)]

JXvMock_6$Gene <- rownames(JXvMock_6)
JXvMock_6 <- data.frame(JXvMock_6)
JXvMock_6$EnsemblID <- ann$ENSID[match(JXvMock_6$Gene, ann$FinalList)]

JXvMock_10$Gene <- rownames(JXvMock_10)
JXvMock_10 <- data.frame(JXvMock_10)
JXvMock_10$EnsemblID <- ann$ENSID[match(JXvMock_10$Gene, ann$FinalList)]

write_xlsx(data.frame(JXvMock_2), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsMock_2dpi.xlsx')
write_xlsx(data.frame(JXvMock_6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsMock_6dpi.xlsx')
write_xlsx(data.frame(JXvMock_10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsMock_10dpi.xlsx')

MNvMock_2$Gene <- rownames(MNvMock_2)
MNvMock_2 <- data.frame(MNvMock_2)
MNvMock_2$EnsemblID <- ann$ENSID[match(MNvMock_2$Gene, ann$FinalList)]

MNvMock_6$Gene <- rownames(MNvMock_6)
MNvMock_6 <- data.frame(MNvMock_6)
MNvMock_6$EnsemblID <- ann$ENSID[match(MNvMock_6$Gene, ann$FinalList)]

MNvMock_10$Gene <- rownames(MNvMock_10)
MNvMock_10 <- data.frame(MNvMock_10)
MNvMock_10$EnsemblID <- ann$ENSID[match(MNvMock_10$Gene, ann$FinalList)]

write_xlsx(data.frame(MNvMock_2), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_MNvsMock_2dpi.xlsx')
write_xlsx(data.frame(MNvMock_6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_MNvsMock_6dpi.xlsx')
write_xlsx(data.frame(MNvMock_10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_MNvsMock_10dpi.xlsx')

LCvMN_2$Gene <- rownames(LCvMN_2)
LCvMN_2 <- data.frame(LCvMN_2)
LCvMN_2$EnsemblID <- ann$ENSID[match(LCvMN_2$Gene, ann$FinalList)]

LCvMN_6$Gene <- rownames(LCvMN_6)
LCvMN_6 <- data.frame(LCvMN_6)
LCvMN_6$EnsemblID <- ann$ENSID[match(LCvMN_6$Gene, ann$FinalList)]

LCvMN_10$Gene <- rownames(LCvMN_10)
LCvMN_10 <- data.frame(LCvMN_10)
LCvMN_10$EnsemblID <- ann$ENSID[match(LCvMN_10$Gene, ann$FinalList)]

write_xlsx(data.frame(LCvMN_2), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_LCvsMN_2dpi.xlsx')
write_xlsx(data.frame(LCvMN_6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_LCvsMN_6dpi.xlsx')
write_xlsx(data.frame(LCvMN_10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_LCvsMN_10dpi.xlsx')

JxvMN_2$Gene <- rownames(JxvMN_2)
JxvMN_2 <- data.frame(JxvMN_2)
JxvMN_2$EnsemblID <- ann$ENSID[match(JxvMN_2$Gene, ann$FinalList)]

JxvMN_6$Gene <- rownames(JxvMN_6)
JxvMN_6 <- data.frame(JxvMN_6)
JxvMN_6$EnsemblID <- ann$ENSID[match(JxvMN_6$Gene, ann$FinalList)]

JxvMN_10$Gene <- rownames(JxvMN_10)
JxvMN_10 <- data.frame(JxvMN_10)
JxvMN_10$EnsemblID <- ann$ENSID[match(JxvMN_10$Gene, ann$FinalList)]

write_xlsx(data.frame(JxvMN_2), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsMN_2dpi.xlsx')
write_xlsx(data.frame(JxvMN_6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsMN_6dpi.xlsx')
write_xlsx(data.frame(JxvMN_10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsMN_10dpi.xlsx')

JxvLC_2$Gene <- rownames(JxvLC_2)
JxvLC_2 <- data.frame(JxvLC_2)
JxvLC_2$EnsemblID <- ann$ENSID[match(JxvLC_2$Gene, ann$FinalList)]

JxvLC_6$Gene <- rownames(JxvLC_6)
JxvLC_6 <- data.frame(JxvLC_6)
JxvLC_6$EnsemblID <- ann$ENSID[match(JxvLC_6$Gene, ann$FinalList)]

JxvLC_10$Gene <- rownames(JxvLC_10)
JxvLC_10 <- data.frame(JxvLC_10)
JxvLC_10$EnsemblID <- ann$ENSID[match(JxvLC_10$Gene, ann$FinalList)]

write_xlsx(data.frame(JxvLC_2), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsLC_2dpi.xlsx')
write_xlsx(data.frame(JxvLC_6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsLC_6dpi.xlsx')
write_xlsx(data.frame(JxvLC_10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_JXvsLC_10dpi.xlsx')

Mock2v6$Gene <- rownames(Mock2v6)
Mock2v6 <- data.frame(Mock2v6)
Mock2v6$EnsemblID <- ann$ENSID[match(Mock2v6$Gene, ann$FinalList)]

Mock2v10$Gene <- rownames(Mock2v10)
Mock2v10 <- data.frame(Mock2v10)
Mock2v10$EnsemblID <- ann$ENSID[match(Mock2v10$Gene, ann$FinalList)]

Mock6v10$Gene <- rownames(Mock6v10)
Mock6v10 <- data.frame(Mock6v10)
Mock6v10$EnsemblID <- ann$ENSID[match(Mock6v10$Gene, ann$FinalList)]

write_xlsx(data.frame(Mock2v6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_Mock_2v6dpi.xlsx')
write_xlsx(data.frame(Mock2v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_Mock_2v10dpi.xlsx')
write_xlsx(data.frame(Mock6v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_Mock_6v10dpi.xlsx')

MN1842v6$Gene <- rownames(MN1842v6)
MN1842v6 <- data.frame(MN1842v6)
MN1842v6$EnsemblID <- ann$ENSID[match(MN1842v6$Gene, ann$FinalList)]

MN1842v10$Gene <- rownames(MN1842v10)
MN1842v10 <- data.frame(MN1842v10)
MN1842v10$EnsemblID <- ann$ENSID[match(MN1842v10$Gene, ann$FinalList)]

MN1846v10$Gene <- rownames(MN1846v10)
MN1846v10 <- data.frame(MN1846v10)
MN1846v10$EnsemblID <- ann$ENSID[match(MN1846v10$Gene, ann$FinalList)]

write_xlsx(data.frame(MN1842v6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_MN184_2v6dpi.xlsx')
write_xlsx(data.frame(MN1842v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_MN184_2v10dpi.xlsx')
write_xlsx(data.frame(MN1846v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_MN184_6v10dpi.xlsx')

L1C52v6$Gene <- rownames(L1C52v6)
L1C52v6 <- data.frame(L1C52v6)
L1C52v6$EnsemblID <- ann$ENSID[match(L1C52v6$Gene, ann$FinalList)]

L1C52v10$Gene <- rownames(L1C52v10)
L1C52v10 <- data.frame(L1C52v10)
L1C52v10$EnsemblID <- ann$ENSID[match(L1C52v10$Gene, ann$FinalList)]

L1C56v10$Gene <- rownames(L1C56v10)
L1C56v10 <- data.frame(L1C56v10)
L1C56v10$EnsemblID <- ann$ENSID[match(L1C56v10$Gene, ann$FinalList)]

write_xlsx(data.frame(L1C52v6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_L1C5_2v6dpi.xlsx')
write_xlsx(data.frame(L1C52v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_L1C5_2v10dpi.xlsx')
write_xlsx(data.frame(L1C56v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_L1C5_6v10dpi.xlsx')


Jxwn062v6$Gene <- rownames(Jxwn062v6)
Jxwn062v6 <- data.frame(Jxwn062v6)
Jxwn062v6$EnsemblID <- ann$ENSID[match(Jxwn062v6$Gene, ann$FinalList)]

Jxwn062v10$Gene <- rownames(Jxwn062v10)
Jxwn062v10 <- data.frame(Jxwn062v10)
Jxwn062v10$EnsemblID <- ann$ENSID[match(Jxwn062v10$Gene, ann$FinalList)]

Jxwn066v10$Gene <- rownames(Jxwn066v10)
Jxwn066v10 <- data.frame(Jxwn066v10)
Jxwn066v10$EnsemblID <- ann$ENSID[match(Jxwn066v10$Gene, ann$FinalList)]

write_xlsx(data.frame(Jxwn062v6), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_Jxwn06_2v6dpi.xlsx')
write_xlsx(data.frame(Jxwn062v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_Jxwn06_2v10dpi.xlsx')
write_xlsx(data.frame(Jxwn066v10), '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/UnfilteredDGEResults_Jxwn06_6v10dpi.xlsx')

## Save DESeq objects ----
saveRDS(vsdb, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/vsdbObject.rds')
saveRDS(d.dds, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/DGEResults/dddsObject.rds')

report(sessionInfo())
#Analyses were conducted using the R Statistical language (version 4.3.2; R Core Team, 2023) on Ubuntu 22.04.3 LTS, using the packages MatrixGenerics (version 1.14.0; Ahlmann-Eltze C et al., 2023), GenomeInfoDb (version 1.38.8; Arora S et al., 2024), matrixStats (version 1.3.0;
#                                                                                                                                                                                                                                                                       Bengtsson H, 2024), EnhancedVolcano (version 1.20.0; Blighe K et al., 2023), lubridate (version 1.9.3; Grolemund G, Wickham H, 2011), Biobase (version 2.62.0; Huber W et al., 2015), BiocGenerics (version 0.48.1; Huber et al., 2015), GenomicRanges (version 1.54.1; Lawrence M et
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               al., 2013), IRanges (version 2.36.0; Lawrence M et al., 2013), DESeq2 (version 1.42.1; Love MI et al., 2014), report (version 0.5.9; Makowski D et al., 2023), SummarizedExperiment (version 1.32.0; Morgan M et al., 2023), tibble (version 3.2.1; Müller K, Wickham H, 2023), writexl
#(version 1.5.0; Ooms J, 2024), S4Vectors (version 0.40.2; Pagès H et al., 2023), ggrepel (version 0.9.5; Slowikowski K, 2024), gplots (version 3.1.3.1; Warnes G et al., 2024), ggplot2 (version 3.5.1; Wickham H, 2016), forcats (version 1.0.0; Wickham H, 2023), stringr (version
#                                                                                                                                                                                                                                                                             1.5.1; Wickham H, 2023), tidyverse (version 2.0.0; Wickham H et al., 2019), readxl (version 1.4.3; Wickham H, Bryan J, 2023), dplyr (version 1.1.4; Wickham H et al., 2023), purrr (version 1.0.2; Wickham H, Henry L, 2023), readr (version 2.1.4; Wickham H et al., 2023) and tidyr
#(version 1.3.0; Wickham H et al., 2023).

#References
#----------
#- Ahlmann-Eltze C, Hickey P, Pagès H (2023). _MatrixGenerics: S4 Generic Summary Statistic Functions that Operate on Matrix-Like Objects_. doi:10.18129/B9.bioc.MatrixGenerics <https://doi.org/10.18129/B9.bioc.MatrixGenerics>, R package version 1.14.0,
#<https://bioconductor.org/packages/MatrixGenerics>.
#- Arora S, Morgan M, Carlson M, Pagès H (2024). _GenomeInfoDb: Utilities for manipulating chromosome names, including modifying them to follow a particular naming style_. R package version 1.38.8, <https://bioconductor.org/packages/GenomeInfoDb>.
#- Bengtsson H (2024). _matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors)_. R package version 1.3.0, <https://CRAN.R-project.org/package=matrixStats>.
#- Blighe K, Rana S, Lewis M (2023). _EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling_. doi:10.18129/B9.bioc.EnhancedVolcano <https://doi.org/10.18129/B9.bioc.EnhancedVolcano>, R package version 1.20.0,
#<https://bioconductor.org/packages/EnhancedVolcano>.
#- Grolemund G, Wickham H (2011). “Dates and Times Made Easy with lubridate.” _Journal of Statistical Software_, *40*(3), 1-25. <https://www.jstatsoft.org/v40/i03/>.
#- Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Ole's AK, Pag`es H, Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M
#(2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” _Nature Methods_, *12*(2), 115-121. <http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
#  - Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M., Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M, MacDonald, J., Obenchain, V., Ole's, K. A, Pag`es, H., Reyes, A.,
#Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan, M. (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” _Nature Methods_, *12*(2), 115-121. <http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
#- Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” _PLoS Computational Biology_, *9*. doi:10.1371/journal.pcbi.1003118 <https://doi.org/10.1371/journal.pcbi.1003118>,
#<http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
#- Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” _PLoS Computational Biology_, *9*. doi:10.1371/journal.pcbi.1003118 <https://doi.org/10.1371/journal.pcbi.1003118>,
#<http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
#- Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” _Genome Biology_, *15*, 550. doi:10.1186/s13059-014-0550-8 <https://doi.org/10.1186/s13059-014-0550-8>.
#- Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” _CRAN_. <https://easystats.github.io/report/>.
#- Morgan M, Obenchain V, Hester J, Pagès H (2023). _SummarizedExperiment: SummarizedExperiment container_. doi:10.18129/B9.bioc.SummarizedExperiment <https://doi.org/10.18129/B9.bioc.SummarizedExperiment>, R package version 1.32.0,
#<https://bioconductor.org/packages/SummarizedExperiment>.
#- Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package version 3.2.1, https://github.com/tidyverse/tibble, <https://tibble.tidyverse.org/>.
#- Ooms J (2024). _writexl: Export Data Frames to Excel 'xlsx' Format_. R package version 1.5.0https://docs.ropensci.org/writexl/ (website) https://github.com/ropensci/writexl (devel) https://libxlsxwriter.github.io (upstream),
#<https://docs.ropensci.org/writexl/%20(website)https://github.com/ropensci/writexl%20(devel)https://libxlsxwriter.github.io%20(upstream)>.
#- Pagès H, Lawrence M, Aboyoun P (2023). _S4Vectors: Foundation of vector-like and list-like containers in Bioconductor_. R package version 0.40.2, <https://bioconductor.org/packages/S4Vectors>.
#- R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#- Slowikowski K (2024). _ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'_. R package version 0.9.5, <https://CRAN.R-project.org/package=ggrepel>.
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
#- Wickham H, Vaughan D, Girlich M (2023). _tidyr: Tidy Messy Data_. R package version 1.3.0, https://github.com/tidyverse/tidyr, <https://tidyr.tidyverse.org>.