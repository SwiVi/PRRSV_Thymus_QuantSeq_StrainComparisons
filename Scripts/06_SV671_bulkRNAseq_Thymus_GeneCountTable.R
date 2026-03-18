lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(readxl)
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(writexl)
library(report)

### Upload count data and target file ----

# Load target file into R ----
target <-read_excel('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_TargetFile.xlsx')
#View(target) # view full target file
head(target)

# Merge count data ----
# This code snippet is from post-doc Tyron Chang
setwd("/project/nadc_prrsv/HTseqCount_SV671_Thymus_TC")
dir <- "/project/nadc_prrsv/HTseqCount_SV671_Thymus_TC/HT_seq_output"
files <- list.files(path=dir,pattern="*_trimmed.fastq.gz_Aligned.out.txt",full.names = TRUE)
length(files)
count_list <- list()## create a list to store these files
total_matches <- 0
for (file in files) {
  df <- read.table(file,header=FALSE,row.names = 1,sep="\t")
  # Ensure df remains a data frame.Otherwise, without", drop = FALSE" it will convert df as a vector.
  df <- df[!grepl("^__", rownames(df)), , drop = FALSE] 
  file_name <- basename(file)
  # Update the regex pattern to be more flexible with possible characters around the pattern
  pattern <- "(\\d+-Mock-\\d+dpi-Thymus_S\\d+)|(\\d+-L1C-5-\\d+dpi-Thymus_S\\d+)|(\\d+-Jxwn06ic-\\d+dpi-Thymus_S\\d+)|(\\d+-LowPath-\\d+dpi-Thymus_S\\d+)"
  sample_name <- str_extract(file_name, pattern)
  
  # calculate the total numbers that match
  total_matches<- total_matches + sum(!is.na(sample_name))
  
  colnames(df)[1] <- sample_name
  df$gene_ID <- rownames(df)
  count_list[[file]] <- df
  
}
print(paste0("Total numbers of samples match for regex is:", " ", total_matches))
head(count_list[[1]])
merged_count <- Reduce(function(x,y) merge(x,y,by="gene_ID",all=TRUE),count_list)

# Convert to custom annotation file gene nomenclature ----
## Convert as many Ensembl IDs as possible to gene symbols from custom annotation file:
ann <- read_excel('/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx')
merged_count <- merged_count %>% mutate(newGeneID = ifelse(startsWith(gene_ID, "ENSSSCG"), ann$FinalList[match(merged_count$gene_ID, ann$ENSID)], merged_count$gene_ID))

# do some manual work to get all rownames to match custom annotation file
setdiff(merged_count$newGeneID, ann$FinalList)
merged_count$newGeneID[merged_count$newGeneID == "AQN_1"] <- "AQN-1"
merged_count$newGeneID[merged_count$newGeneID == "CFAP298_TCP10L"] <- "CFAP298-TCP10L"
merged_count$newGeneID[merged_count$newGeneID == "CRSP_2"] <- "CRSP-2"
merged_count$newGeneID[merged_count$newGeneID == "FFAR2_L"] <- "FFAR2-L"
merged_count$newGeneID[merged_count$newGeneID == "GP91_PHOX_36246"] <- "GP91-PHOX_36246"
merged_count$newGeneID[merged_count$newGeneID == "H1_2"] <- "H1-2"
merged_count$newGeneID[merged_count$newGeneID == "MIR4331_2"] <- "MIR4331-2"
merged_count$newGeneID[merged_count$newGeneID == "MT_2B"] <- "MT-2B"
merged_count$newGeneID[merged_count$newGeneID == "NKX6_3"] <- "NKX6-3"
merged_count$newGeneID[merged_count$newGeneID == "OLF42_3"] <- "OLF42-3"
merged_count$newGeneID[merged_count$newGeneID == "PIGE_108A11.3"] <- "PIGE-108A11.3"
merged_count$newGeneID[merged_count$newGeneID == "RPL36A_HNRNPH2"] <- "RPL36A-HNRNPH2"
merged_count$newGeneID[merged_count$newGeneID == "SERPINA3_2_30371"] <- "SERPINA3-2_30371"
merged_count$newGeneID[merged_count$newGeneID == "SLA_2"] <- "SLA-2"
merged_count$newGeneID[merged_count$newGeneID == "SLA_8"] <- "SLA-8"
merged_count$newGeneID[merged_count$newGeneID == "SLA_DQA1"] <- "SLA-DQA1"
merged_count$newGeneID[merged_count$newGeneID == "SLA_DRB1"] <- "SLA-DRB1"
setdiff(merged_count$newGeneID, ann$FinalList)
rownames(merged_count) <- merged_count$newGeneID #set gene_ID as rownames
merged_count$gene_ID <- NULL #drop gene_ID columns
merged_count$newGeneID <- NULL #drop gene_ID columns
head(merged_count)

# Assess variability of technical replicates (biological replicates sequenced across two lanes): ----
counts <- merged_count
sums <- data.frame(colSums(counts))
barplot(sums$colSums.counts.) # check that every two samples have similar reads detected (indicating biological replicate sequenced across two lanes)
cor <- cor(counts)
cor <- melt(cor)
head(cor)
ggplot(data = cor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = .9, limit = c(.8, 1), oob=scales::squish) # should see extremely high correlations for every two samples (indicating biological replicate sequenced across two lanes)
# All looks well across technical replicates

# Merge technical replicates: ----
colnames(counts) <- str_extract(colnames(counts), "[^_]+")
target$name <- paste(paste(paste(target$AnimalID, target$Treatment, sep = '-'), target$DaysPostInfection, sep = '-'), 'Thymus', sep = '-')
target$name <- gsub("Jxwn06", "Jxwn06ic", target$name)
target$name <- gsub("MN184", "LowPath", target$name)
target$name <- gsub("L1C5", "L1C-5", target$name)
identical(colnames(counts), target$name) # make sure the sample IDs match (identical = TRUE) before proceeding!

counts <- t(rowsum(t(counts), group = colnames(counts), na.rm = T)) # sum counts for columns with identical names

sums <- data.frame(colSums(counts))
barplot(sums$colSums.counts.) # see total library size per sample

# Further refine target file: ----
target$group <- paste(target$Treatment, target$DaysPostInfection, sep = '_')
#colnames(target)
target <- target[,c(1:3,6, 7)] # subset to only experimental variables of animal ID, treatment, and dpi
target <- distinct(target) # remove duplicate rows
target <- as.data.frame(target)
head(target)
target <- target[match(colnames(counts), target$name),]

identical(colnames(counts), target$name) # make sure the sample IDs match (identical = TRUE)

# Remove known non-thymus sample: ----
## This was validated via assessment of matching H&E tissue
target <- target[target$AnimalID != '854-Jxwn06ic-10dpi-Thymus', ]
counts <- subset(counts, select = -c(`854-Jxwn06ic-10dpi-Thymus`))

target <- target[match(colnames(counts), target$name),]
identical(colnames(counts), target$name) # make sure the sample IDs match (identical = TRUE)

# Save refined target file and gene counts files:
write_xlsx(target, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_TargetFile_Refined.xlsx')
write.csv(counts, '/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_GeneCounts_Refined.csv', row.names = TRUE)
#counts <- read.csv('/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/MetaData/SV671_bulkRNAseq_Thymus_GeneCounts_Refined.csv', check.names = FALSE, row.names = 1)

report(sessionInfo())
#Analyses were conducted using the R Statistical language (version 4.3.2; R Core Team, 2023) on Ubuntu 22.04.3 LTS, using the packages report (version 0.5.9; Makowski D et al., 2023), writexl (version 1.5.0; Ooms J, 2024), reshape2 (version 1.4.4; Wickham H, 2007), ggplot2 (version
#                                                                                                                                                                                                                                                                                 3.5.1; Wickham H, 2016), stringr (version 1.5.1; Wickham H, 2023), readxl (version 1.4.3; Wickham H, Bryan J, 2023), dplyr (version 1.1.4; Wickham H et al., 2023) and scales (version 1.3.0; Wickham H et al., 2023).

#References
#----------
#- Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” _CRAN_. <https://easystats.github.io/report/>.
#- Ooms J (2024). _writexl: Export Data Frames to Excel 'xlsx' Format_. R package version 1.5.0https://docs.ropensci.org/writexl/ (website) https://github.com/ropensci/writexl (devel) https://libxlsxwriter.github.io (upstream),
#<https://docs.ropensci.org/writexl/%20(website)https://github.com/ropensci/writexl%20(devel)https://libxlsxwriter.github.io%20(upstream)>.
#- R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#- Wickham H (2007). “Reshaping Data with the reshape Package.” _Journal of Statistical Software_, *21*(12), 1-20. <http://www.jstatsoft.org/v21/i12/>.
#- Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
#- Wickham H (2023). _stringr: Simple, Consistent Wrappers for Common String Operations_. R package version 1.5.1, https://github.com/tidyverse/stringr, <https://stringr.tidyverse.org>.
#- Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package version 1.4.3, https://github.com/tidyverse/readxl, <https://readxl.tidyverse.org>.
#- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.4, https://github.com/tidyverse/dplyr, <https://dplyr.tidyverse.org>.
#- Wickham H, Pedersen T, Seidel D (2023). _scales: Scale Functions for Visualization_. R package version 1.3.0, https://github.com/r-lib/scales, <https://scales.r-lib.org>.