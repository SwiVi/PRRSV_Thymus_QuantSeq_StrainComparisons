lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(AUCell)
library(hypeR)
library(data.table)
library(writexl)
library(readxl)
library(tidyverse)
library(gplots)
library(ggplot2)
library(scales)
library(report)

## Make hallmark gene sets ----
#msigdb_info()
hallmark <- msigdb_gsets(species="Homo sapiens", collection="H", clean = TRUE) # extract hallmark gene sets from MSigDB (these are human gene sets)
#print(hallmark$genesets)

setlist = rbindlist(
  lapply(hallmark$genesets, function(x) data.table(x)),
  fill = TRUE, idcol = TRUE
)
colnames(setlist) <- c('hallmark', 'gene')

write_xlsx(setlist, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/GSEA_HallmarkGeneSets_human.xlsx')
#Accessed https://biit.cs.ut.ee/gprofiler_archive3/e97_eg44_p13/orth to identify Ensembl IDs of pig ortholog genes from v97 Ensembl annotation by copying in human gene symbols from the .xlsx file above. That list is in the .csv file below.
#Manually added additional column for hallmark process to the file, then exported as .xlsx
orthos <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/20260227_gProfiler_v97archive_GSEA_HallmarkgeneSets_humanToPigOrthologs.xlsx')

anngenes <- read_xlsx('/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # load in gene naming file that was used for custom genome annotation
orthos$gene = anngenes$FinalList[match(orthos$ortholog_ensg, anngenes$ENSID)] # identify gene names used in annotation file

geneSets <- orthos %>% pivot_wider(names_from = hallmark_process, values_from = gene)
geneSets <- geneSets[,c((ncol(geneSets)-50+1):ncol(geneSets))]
write_xlsx(geneSets, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/GSEA_HallmarkGeneSets_HumanToPigOrthos.xlsx')

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
write_xlsx(AUC, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/GSEA_HallmarkGeneSets_AUCScores.xlsx')

## Make heatmaps (all samples hierarchical) ----
my_palette <- colorRampPalette(c("navy", "grey95", "red2"))(n = 399)
colors = seq(-4, 4,length=400)
heatmap.2(as.matrix(t(AUC[1:50])), Rowv=T, Colv=T, dendrogram= "both", margins = c(12,12), breaks = colors,
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # swap out cols for different treatment variables

cols <- c(rep("#F8766D", 14), rep("#7CAE00", 15), rep("#00BFC4", 15), rep("#C77CFF", 15))
heatmap.2(as.matrix(t(AUC[1:50])), Rowv=T, dendrogram= "both", margins = c(12,12), breaks = colors,
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette, ColSideColors = cols) # swap out cols for different treatment variables
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Mock", "MN184", "L1C.5", "Jxwn06"), # category labels
       col = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),  # color key
       lty= 1, lwd = 10) 

cols <- c(rep("grey80", 5), rep("grey40", 5), rep("black", 4), rep("grey80", 5), rep("grey40", 5), rep("black", 5), 
          rep("grey80", 5), rep("grey40", 5), rep("black", 5), rep("grey80", 5), rep("grey40", 5), rep("black", 5))
heatmap.2(as.matrix(t(AUC[1:50])), Rowv=T, dendrogram= "both", margins = c(12,12), breaks = colors,
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

### Adipogenesis ----
sub <- subset(AUC, geneSet == 'Adipogenesis')
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
  ggtitle("Adipogenesis") +
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

### Allograft Rejection ----
sub <- subset(AUC, geneSet == 'Allograft.Rejection')
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
  ggtitle("Allograft Rejection") +
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

### Androgen Response ----
sub <- subset(AUC, geneSet == 'Androgen.Response')
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
  ggtitle("Androgen Response") +
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

### Angiogenesis ----
sub <- subset(AUC, geneSet == 'Angiogenesis')
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
  ggtitle("Angiogenesis") +
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

### Apical Junction ----
sub <- subset(AUC, geneSet == 'Apical.Junction')
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
  ggtitle("Apical Junction") +
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

### Apical Surface ----
sub <- subset(AUC, geneSet == 'Apical.Surface')
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
  ggtitle("Apical Surface") +
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

### Apoptosis ----
sub <- subset(AUC, geneSet == 'Apoptosis')
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
  ggtitle("Apoptosis") +
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

### Bile Acid Metabolism ----
sub <- subset(AUC, geneSet == 'Bile.Acid.Metabolism')
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
  ggtitle("Bile Acid Metabolism") +
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

### Cholesterol Homeostasis ----
sub <- subset(AUC, geneSet == 'Cholesterol.Homeostasis')
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
  ggtitle("Cholesterol Homeostasis") +
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

### Coagulation ----
sub <- subset(AUC, geneSet == 'Coagulation')
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
  ggtitle("Coagulation") +
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

### Complement ----
sub <- subset(AUC, geneSet == 'Complement')
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
  ggtitle("Complement") +
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

### DNA Repair ----
sub <- subset(AUC, geneSet == 'Dna.Repair')
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
  ggtitle("DNA Repair") +
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

### E2F Targets ----
sub <- subset(AUC, geneSet == 'E2f.Targets')
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
  ggtitle("E2F Targets") +
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

### Epithelial Mesenchymal Transition ----
sub <- subset(AUC, geneSet == 'Epithelial.Mesenchymal.Transition')
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
  ggtitle("Epithelial Mesenchymal Transition") +
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

### Estrogen Response Early ----
sub <- subset(AUC, geneSet == 'Estrogen.Response.Early')
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
  ggtitle("Estrogen Response Early") +
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

### Estrogen Response Late ----
sub <- subset(AUC, geneSet == 'Estrogen.Response.Late')
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
  ggtitle("Estrogen Response Late") +
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

### Fatty Acid Metabolism ----
sub <- subset(AUC, geneSet == 'Fatty.Acid.Metabolism')
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
  ggtitle("Fatty Acid Metabolism") +
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

### G2M Checkpoint ----
sub <- subset(AUC, geneSet == 'G2m.Checkpoint')
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
  ggtitle("G2M Checkpoint") +
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

### Glycolysis ----
sub <- subset(AUC, geneSet == 'Glycolysis')
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
  ggtitle("Glycolysis") +
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

### Hedgehog Signaling ----
sub <- subset(AUC, geneSet == 'Hedgehog.Signaling')
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
  ggtitle("Hedgehog Signaling") +
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

### Heme Metabolism ----
sub <- subset(AUC, geneSet == 'Heme.Metabolism')
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
  ggtitle("Heme Metabolism") +
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

### Hypoxia ----
sub <- subset(AUC, geneSet == 'Hypoxia')
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
  ggtitle("Hypoxia") +
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

### IL-2/STAT5 Signaling ----
sub <- subset(AUC, geneSet == 'Il2.Stat5.Signaling')
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
  ggtitle("IL-2/STAT5 Signaling") +
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

### IL-6/JAK/STAT3 Signaling ----
sub <- subset(AUC, geneSet == 'Il6.Jak.Stat3.Signaling')
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
  ggtitle("IL-6/JAK/STAT3 Signaling") +
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

### Inflammatory Response ----
sub <- subset(AUC, geneSet == 'Inflammatory.Response')
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
  ggtitle("Inflammatory Response") +
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

### IFN-a Response ----
sub <- subset(AUC, geneSet == 'Interferon.Alpha.Response')
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
  ggtitle("IFN-a Response") +
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

### IFN-g Response ----
sub <- subset(AUC, geneSet == 'Interferon.Gamma.Response')
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
  ggtitle("IFN-g Response") +
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

### Kras Signaling (Down) ----
sub <- subset(AUC, geneSet == 'Kras.Signaling.Dn')
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
  ggtitle("Kras Signaling (Down)") +
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

### Kras Signaling (Up) ----
sub <- subset(AUC, geneSet == 'Kras.Signaling.Up')
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
  ggtitle("Kras Signaling (Up)") +
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

### Mitotic Spindle ----
sub <- subset(AUC, geneSet == 'Mitotic.Spindle')
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
  ggtitle("Mitotic Spindle") +
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

### MTORC1 Signaling ----
sub <- subset(AUC, geneSet == 'Mtorc1.Signaling')
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
  ggtitle("MTORC1 Signaling") +
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

### MYC Targets (v1) ----
sub <- subset(AUC, geneSet == 'Myc.Targets.V1')
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
  ggtitle("MYC Targets (v1)") +
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

### MYC Targets (v2) ----
sub <- subset(AUC, geneSet == 'Myc.Targets.V2')
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
  ggtitle("MYC Targets (v2)") +
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

### Myogenesis ----
sub <- subset(AUC, geneSet == 'Myogenesis')
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
  ggtitle("Myogenesis") +
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

### Notch Signaling ----
sub <- subset(AUC, geneSet == 'Notch.Signaling')
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
  ggtitle("Notch Signaling") +
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

### Oxidative Phosphorylation ----
sub <- subset(AUC, geneSet == 'Oxidative.Phosphorylation')
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
  ggtitle("Oxidative Phosphorylation") +
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

### p53 Pathway ----
sub <- subset(AUC, geneSet == 'P53.Pathway')
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
  ggtitle("p53 Pathway") +
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

### Pancreas Beta Cells ----
sub <- subset(AUC, geneSet == 'Pancreas.Beta.Cells')
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
  ggtitle("Pancreas Beta Cells") +
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

### Peroxisome ----
sub <- subset(AUC, geneSet == 'Peroxisome')
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
  ggtitle("Peroxisome") +
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

### PI3K/AKT/MTOR Signaling ----
sub <- subset(AUC, geneSet == 'Pi3k.Akt.Mtor.Signaling')
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
  ggtitle("PI3K/AKT/MTOR Signaling") +
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

### Protein Secretion ----
sub <- subset(AUC, geneSet == 'Protein.Secretion')
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
  ggtitle("Protein Secretion") +
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

### Reactive Oxygen Species Pathway ----
sub <- subset(AUC, geneSet == 'Reactive.Oxygen.Species.Pathway')
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
  ggtitle("Reactive Oxygen Species Pathway") +
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

### Spermatogenesis ----
sub <- subset(AUC, geneSet == 'Spermatogenesis')
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
  ggtitle("Spermatogenesis") +
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

### TGF-b Signaling ----
sub <- subset(AUC, geneSet == 'Tgf.Beta.Signaling')
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
  ggtitle("TGF-b Signaling") +
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

### TNF-a Signaling Via NFk-b ----
sub <- subset(AUC, geneSet == 'Tnfa.Signaling.Via.Nfkb')
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
  ggtitle("TNF-a Signaling Via NFk-b") +
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

### Unfolded Protein Response ----
sub <- subset(AUC, geneSet == 'Unfolded.Protein.Response')
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
  ggtitle("Unfolded Protein Response") +
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

### UV Response (Down) ----
sub <- subset(AUC, geneSet == 'Uv.Response.Dn')
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
  ggtitle("UV Response (Down)") +
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

### UV Response (Up) ----
sub <- subset(AUC, geneSet == 'Uv.Response.Up')
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
  ggtitle("UV Response (Up)") +
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

### Wnt/b-catenin Signaling ----
sub <- subset(AUC, geneSet == 'Wnt.Beta.Catenin.Signaling')
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
  ggtitle("Wnt/b-catenin Signaling") +
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

### Xenobiotic Metabolism ----
sub <- subset(AUC, geneSet == 'Xenobiotic.Metabolism')
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
  ggtitle("Xenobiotic Metabolism") +
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
#Analyses were conducted using the R Statistical language (version 4.3.2; R Core Team, 2023) on Ubuntu 22.04.3 LTS, using the packages AUCell (version 1.24.0; Aibar S et al., 2017), data.table (version 1.14.10; Barrett T et al., 2023), hypeR (version
#                                                                                                                                                                                                                                                  2.8.2; Federico A, Monti S, 2020), lubridate (version 1.9.3; Grolemund G, Wickham H, 2011), report (version 0.5.9; Makowski D et al., 2023), tibble (version 3.2.1; Müller K, Wickham H, 2023), writexl (version 1.5.0; Ooms J, 2024), gplots (version 3.1.3.1;
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Warnes G et al., 2024), ggplot2 (version 3.5.1; Wickham H, 2016), forcats (version 1.0.0; Wickham H, 2023), stringr (version 1.5.1; Wickham H, 2023), tidyverse (version 2.0.0; Wickham H et al., 2019), readxl (version 1.4.3; Wickham H, Bryan J, 2023),
#dplyr (version 1.1.4; Wickham H et al., 2023), purrr (version 1.0.2; Wickham H, Henry L, 2023), readr (version 2.1.4; Wickham H et al., 2023), scales (version 1.3.0; Wickham H et al., 2023) and tidyr (version 1.3.0; Wickham H et al., 2023).

#References
#----------
#  - Aibar S, Bravo Gonzalez-Blas C, Moerman T, Huynh-Thu V, Imrichova H, Hulselmans G, Rambow F, Marine J, Geurts P, Aerts J, van den Oord J, Kalender Atak Z, Wouters J, Aerts S (2017). “SCENIC: Single-Cell Regulatory Network Inference And Clustering.”
#_Nature Methods_, *14*, 1083-1086. doi:10.1038/nmeth.4463 <https://doi.org/10.1038/nmeth.4463>. Aibar S, Aerts S (2016). “AUCell.”
#- Barrett T, Dowle M, Srinivasan A (2023). _data.table: Extension of `data.frame`_. R package version 1.14.10, https://Rdatatable.gitlab.io/data.table, https://github.com/Rdatatable/data.table, <https://r-datatable.com>.
#- Federico A, Monti S (2020). “hypeR: an R package for geneset enrichment workflows.” _Bioinformatics_, *36*(4), 1307-1308. <https://doi.org/10.1093/bioinformatics/btz700>.
#- Grolemund G, Wickham H (2011). “Dates and Times Made Easy with lubridate.” _Journal of Statistical Software_, *40*(3), 1-25. <https://www.jstatsoft.org/v40/i03/>.
#- Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” _CRAN_. <https://easystats.github.io/report/>.
#- Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package version 3.2.1, https://github.com/tidyverse/tibble, <https://tibble.tidyverse.org/>.
#- Ooms J (2024). _writexl: Export Data Frames to Excel 'xlsx' Format_. R package version 1.5.0https://docs.ropensci.org/writexl/ (website) https://github.com/ropensci/writexl (devel) https://libxlsxwriter.github.io (upstream),
#<https://docs.ropensci.org/writexl/%20(website)https://github.com/ropensci/writexl%20(devel)https://libxlsxwriter.github.io%20(upstream)>.
#- R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#- Warnes G, Bolker B, Bonebakker L, Gentleman R, Huber W, Liaw A, Lumley T, Maechler M, Magnusson A, Moeller S, Schwartz M, Venables B (2024). _gplots: Various R Programming Tools for Plotting Data_. R package version 3.1.3.1,
#<https://CRAN.R-project.org/package=gplots>.
#- Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
#- Wickham H (2023). _forcats: Tools for Working with Categorical Variables (Factors)_. R package version 1.0.0, https://github.com/tidyverse/forcats, <https://forcats.tidyverse.org/>.
#- Wickham H (2023). _stringr: Simple, Consistent Wrappers for Common String Operations_. R package version 1.5.1, https://github.com/tidyverse/stringr, <https://stringr.tidyverse.org>.
#- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019).
#“Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
#- Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package version 1.4.3, https://github.com/tidyverse/readxl, <https://readxl.tidyverse.org>.
#- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.4, https://github.com/tidyverse/dplyr, <https://dplyr.tidyverse.org>.
#- Wickham H, Henry L (2023). _purrr: Functional Programming Tools_. R package version 1.0.2, https://github.com/tidyverse/purrr, <https://purrr.tidyverse.org/>.
#- Wickham H, Hester J, Bryan J (2023). _readr: Read Rectangular Text Data_. R package version 2.1.4, https://github.com/tidyverse/readr, <https://readr.tidyverse.org>.
#- Wickham H, Pedersen T, Seidel D (2023). _scales: Scale Functions for Visualization_. R package version 1.3.0, https://github.com/r-lib/scales, <https://scales.r-lib.org>.
#- Wickham H, Vaughan D, Girlich M (2023). _tidyr: Tidy Messy Data_. R package version 1.3.0, https://github.com/tidyverse/tidyr, <https://tidyr.tidyverse.org>.