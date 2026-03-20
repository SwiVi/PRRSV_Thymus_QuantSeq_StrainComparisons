library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(data.table)
library(tidyr)
library(rstatix)
library(writexl)
library(ggpmisc)

# Atrophy score ----

## Stats ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.numeric(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)

## Kruskal-Wallis rank sum test for multi-group, non-parametric, unpaired alternative to the ANOVA
df$TreatmentxDaysPostInfection <- paste(df$Treatment, df$DaysPostInfection, sep = '_dpi')
kruskal <- df %>% kruskal_test(CorticalAtrophyScore ~ TreatmentxDaysPostInfection)
kruskal # if the p-value is significant, can move on to running post-hoc multiple comparison tests

##Post-hoc Mann-Whitney U test for multi-group, non-parametric, unpaired multiple comparisons
## For unpaired comparison across all groups
## Only select relevant comparisons of interest to look at. Do not look at any comparisons made within a treatment group across time, as these should instead be run with paired comparisons.
## No p value correction used. We pre-determined to only look at time-matched mock vs treatment group comparisons

df2 <- split(df, df$DaysPostInfection)
wilcox <- lapply(df2, function(g) pairwise.wilcox.test(g$CorticalAtrophyScore, g$Treatment, p.adjust.method="fdr"))
wilcox

df2 <- split(df, df$Treatment)
wilcox <- lapply(df2, function(g) pairwise.wilcox.test(g$CorticalAtrophyScore, g$DaysPostInfection, p.adjust.method="fdr"))
wilcox

## Plotting ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.numeric(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)
ggplot(df, aes(y= CorticalAtrophyScore, x=Treatment)) +  
  #geom_boxplot(outlier.shape = NA, alpha = 0.3, coef= 0, aes(fill = Treatment, color = Treatment))+ 
  geom_jitter(width = 0.3, height = 0, alpha = 0.3, size = 2, shape = 21, aes(fill = Treatment, color = Treatment)) + 
  stat_summary(fun = mean, 
               geom = "point", size = 10, shape=45,
               aes(fill = Treatment, color = Treatment)) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", size = 0.8, width = 0.4,
               aes(fill = Treatment, color = Treatment)) +
  scale_y_continuous(limits = c(0, 3.75), breaks = seq(0, 3, by = 1)) +
  scale_x_discrete(limits=rev) +
  facet_wrap(~DaysPostInfection, scale = 'fixed', 
             labeller = as_labeller(c(`2` = "2dpi", `6` = "6dpi", `10` = "10dpi"))) + 
  theme_bw() +
  ylab('Atrophy Score') +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'),
        strip.text = element_text(size = 15, face = "bold")) 

ggplot(df, aes(x=DaysPostInfection, y = CorticalAtrophyScore, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(df[(df$DaysPostInfection == 2 | df$DaysPostInfection == 6 | df$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, CorticalAtrophyScore),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(df[(df$DaysPostInfection == 2 | df$DaysPostInfection == 6 | df$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, CorticalAtrophyScore, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  xlab('Days Post-Infection') +
  ylab('Atrophy Score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

# PRRSV IHC ----

## Stats ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.numeric(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)

## Kruskal-Wallis rank sum test for multi-group, non-parametric, unpaired alternative to the ANOVA
df$TreatmentxDaysPostInfection <- paste(df$Treatment, df$DaysPostInfection, sep = '_dpi')
kruskal <- df %>% kruskal_test(PRRSVIHCScore ~ TreatmentxDaysPostInfection)
kruskal # if the p-value is significant, can move on to running post-hoc multiple comparison tests

##Post-hoc Mann-Whitney U test for multi-group, non-parametric, unpaired multiple comparisons
## For unpaired comparison across all groups
## Only select relevant comparisons of interest to look at. Do not look at any comparisons made within a treatment group across time, as these should instead be run with paired comparisons.
## No p value correction used. We pre-determined to only look at time-matched mock vs treatment group comparisons

df2 <- split(df, df$DaysPostInfection)
wilcox <- lapply(df2, function(g) pairwise.wilcox.test(g$PRRSVIHCScore, g$Treatment, p.adjust.method="fdr"))
wilcox

df2 <- split(df, df$Treatment)
wilcox <- lapply(df2, function(g) pairwise.wilcox.test(g$PRRSVIHCScore, g$DaysPostInfection, p.adjust.method="fdr"))
wilcox

## Plotting ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.numeric(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)
ggplot(df, aes(y= PRRSVIHCScore, x=Treatment)) +  
  #geom_boxplot(outlier.shape = NA, alpha = 0.3, coef= 0, aes(fill = Treatment, color = Treatment))+ 
  geom_jitter(width = 0.3, height = 0, alpha = 0.3, size = 2, shape = 21, aes(fill = Treatment, color = Treatment)) + 
  stat_summary(fun = mean, 
               geom = "point", size = 10, shape=45,
               aes(fill = Treatment, color = Treatment)) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", size = 0.8, width = 0.4,
               aes(fill = Treatment, color = Treatment)) +
  scale_y_continuous(limits = c(0, 4.75), breaks = seq(0, 4, by = 1)) +
  scale_x_discrete(limits=rev) +
  facet_wrap(~DaysPostInfection, scale = 'fixed', 
             labeller = as_labeller(c(`2` = "2dpi", `6` = "6dpi", `10` = "10dpi"))) + 
  theme_bw() +
  ylab('PRRSV IHC Score') +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'),
        strip.text = element_text(size = 15, face = "bold")) 

ggplot(df, aes(x=DaysPostInfection, y = PRRSVIHCScore, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(df[(df$DaysPostInfection == 2 | df$DaysPostInfection == 6 | df$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, PRRSVIHCScore),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(df[(df$DaysPostInfection == 2 | df$DaysPostInfection == 6 | df$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, PRRSVIHCScore, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  xlab('Days Post-Infection') +
  ylab('PRRSV IHC Score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'))

# TUNEL ----

## Stats ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.numeric(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)

## Kruskal-Wallis rank sum test for multi-group, non-parametric, unpaired alternative to the ANOVA
df$TreatmentxDaysPostInfection <- paste(df$Treatment, df$DaysPostInfection, sep = '_dpi')
kruskal <- df %>% kruskal_test(TUNELScore ~ TreatmentxDaysPostInfection)
kruskal # if the p-value is significant, can move on to running post-hoc multiple comparison tests

##Post-hoc Mann-Whitney U test for multi-group, non-parametric, unpaired multiple comparisons
## For unpaired comparison across all groups
## Only select relevant comparisons of interest to look at. Do not look at any comparisons made within a treatment group across time, as these should instead be run with paired comparisons.
## No p value correction used. We pre-determined to only look at time-matched mock vs treatment group comparisons

df2 <- split(df, df$DaysPostInfection)
wilcox <- lapply(df2, function(g) pairwise.wilcox.test(g$TUNELScore, g$Treatment, p.adjust.method="fdr"))
wilcox

df2 <- split(df, df$Treatment)
wilcox <- lapply(df2, function(g) pairwise.wilcox.test(g$TUNELScore, g$DaysPostInfection, p.adjust.method="fdr"))
wilcox

## Plotting ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.numeric(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)
ggplot(df, aes(y= TUNELScore, x=Treatment)) +  
  #geom_boxplot(outlier.shape = NA, alpha = 0.3, coef= 0, aes(fill = Treatment, color = Treatment))+ 
  geom_jitter(width = 0.3, height = 0, alpha = 0.3, size = 2, shape = 21, aes(fill = Treatment, color = Treatment)) + 
  stat_summary(fun = mean, 
               geom = "point", size = 10, shape=45,
               aes(fill = Treatment, color = Treatment)) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", size = 0.8, width = 0.4,
               aes(fill = Treatment, color = Treatment)) +
  scale_y_continuous(limits = c(0, 4.75), breaks = seq(0, 4, by = 1)) +
  scale_x_discrete(limits=rev) +
  facet_wrap(~DaysPostInfection, scale = 'fixed', 
             labeller = as_labeller(c(`2` = "2dpi", `6` = "6dpi", `10` = "10dpi"))) + 
  theme_bw() +
  ylab('TUNEL Score') +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold'),
        strip.text = element_text(size = 15, face = "bold")) 

ggplot(df, aes(x=DaysPostInfection, y = TUNELScore, color = Treatment, group = Treatment)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.02), alpha = .4) + 
  #geom_line(alpha = 0.2, aes(group = AnimalID)) +
  scale_x_continuous(limits = c(-.5, 10.5), breaks = seq(0, 10, 1), minor_breaks = NULL) +
  theme_bw() +
  stat_summary(df[(df$DaysPostInfection == 2 | df$DaysPostInfection == 6 | df$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, TUNELScore),
               fun = mean, 
               geom = "line", size = 1) +
  stat_summary(df[(df$DaysPostInfection == 2 | df$DaysPostInfection == 6 | df$DaysPostInfection == 10),],
               mapping = aes(DaysPostInfection, TUNELScore, fill = Treatment),
               fun = mean, shape = 21,
               geom = "point", size = 3, color = 'black') +
  xlab('Days Post-Infection') +
  ylab('TUNEL Score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold')) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1))

# Correlation plot ----
df <- read_xlsx('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/SV671_ThymusPathologyScoring.xlsx')
colnames(df)
df$DaysPostInfection <- gsub('.{3}$', '', df$DaysPostInfection)
df$DaysPostInfection <- as.character(df$DaysPostInfection)
df$AnimalID <- as.character(df$AnimalID)

ggplot(df, aes(y = CorticalAtrophyScore, x = TUNELScore)) + 
  geom_point(position = position_jitter(width = 0.04, height = 0.04), alpha = .4, size = 1, fill = 'grey30') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  ylab('Atrophy Score') +
  xlab('TUNEL Score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold')) +
  stat_poly_eq(method = 'lm', use_label("eq")) +
  stat_poly_eq(method = 'lm', use_label("R2"), label.y = 0.9) +
  stat_poly_eq(method = 'lm', use_label("p"), label.y = 0.82) +
  stat_poly_line(method = 'lm', color = 'red3', se = TRUE) 

ggplot(df, aes(y = CorticalAtrophyScore, x = PRRSVIHCScore)) + 
  geom_point(position = position_jitter(width = 0.04, height = 0.04), alpha = .4, size = 1, fill = 'grey30') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  ylab('Atrophy Score') +
  xlab('PRRSV IHC Score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold')) +
  stat_poly_eq(method = 'lm', use_label("eq")) +
  stat_poly_eq(method = 'lm', use_label("R2"), label.y = 0.9) +
  stat_poly_eq(method = 'lm', use_label("p"), label.y = 0.82) +
  stat_poly_line(method = 'lm', color = 'red3', se = TRUE) 

ggplot(df, aes(x = TUNELScore, y = PRRSVIHCScore)) + 
  geom_point(position = position_jitter(width = 0.04, height = 0.04), alpha = .4, size = 1, fill = 'grey30') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  xlab('TUNEL Score') +
  ylab('PRRSV IHC Score') +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 20, face = 'bold'),
        title = element_text(size = 20, face = 'bold')) +
  stat_poly_eq(method = 'lm', use_label("eq")) +
  stat_poly_eq(method = 'lm', use_label("R2"), label.y = 0.9) +
  stat_poly_eq(method = 'lm', use_label("p"), label.y = 0.82) +
  stat_poly_line(method = 'lm', color = 'red3', se = TRUE) 

