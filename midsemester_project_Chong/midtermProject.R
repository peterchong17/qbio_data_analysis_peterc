library(TCGAbiolinks)
library(SummarizedExperiment)
library("DESeq2")
library(survival)
library(survminer)
getwd()
setwd("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data")
getwd()

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)
counts = assays(sum_exp)$"HTSeq - Counts"

colData(sum_exp)$gender # category that I will be using to compare expression between
bool_gender_na = is.na(colData(sum_exp)$gender) # create a vector that isolates which rows have NA

# cleaned colData without the two patients that were missing gender information
gender_no_NAs = colData(sum_exp)$gender[!bool_gender_na]

# create a mask only containing T for the gene we are searching for, CFTR
name_mask = (rowData(sum_exp)$external_gene_name == "CFTR")

# creates a subset of counts that only contains samples of CFTR
gene_counts = counts[name_mask, !bool_gender_na]

# makes comparative boxplot
box=boxplot(gene_counts ~ gender_no_NAs, xlab = "Gender", ylab = "Counts of CFTR")
box$stats # summarizes the stats for the two boxplots to find the medians. For female, it is 5376 and for male it is 5742

# makes bar graph
barplot(table(gender_no_NAs), 
        col = rainbow(2), 
        xlab = "Gender", ylab = "Counts of CFTR", main = "CFTR Counts Between Gender")

#testing
sum(gender_no_NAs=="male") # number of male samples (272)
sum(gender_no_NAs=="female") # number of female samples (247)

head(sum_exp)
#SURVIVAL PLOT
# if there is no days_to_death information, then it uses days_to_last_follow_up as an estimate
sum_exp$days_to_death = ifelse(is.na(sum_exp$days_to_death), sum_exp$days_to_last_follow_up, sum_exp$days_to_death)
# removes rows that do not have info
cleaned_sum_exp = sum_exp[!is.na(sum_exp$days_to_death),]
cleaned_sum_exp$days_to_death

sum(is.na(sum_exp$days_to_death))
sum_exp$death_event = as.integer(sum_exp$vital_status == "Dead")

surv_object <- Surv(time = cleaned_clinic$days_to_death, 
                    event = cleaned_clinic$death_event)



######
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
#BiocManager::install("DESeq2")
library(TCGAbiolinks)
library(SummarizedExperiment)
library("DESeq2")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)

head(colData(sum_exp))
is.na(colData(sum_exp))
patient_data = colData(sum_exp)
counts = assays(sum_exp)$"HTSeq - Counts"

#clean counts of the gender na's
patient_data_na_mask = is.na(colData(sum_exp)$gender)
counts = counts[ , !(patient_data_na_mask) ]
patient_data = patient_data[ !(patient_data_na_mask) , ]

if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
counts_row_sums = rowSums(counts)
low_counts_mask = ifelse(counts_row_sums < 10, FALSE, TRUE)
sum(low_counts_mask) 
counts = counts[ low_counts_mask, ]
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = patient_data, 
                             design = ~gender)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the female vs. male comparison
results = results(dds_obj, format = "DataFrame", contrast = c("gender", "female", "male"))

log2FoldChange_threshold =  ifelse(results$log2FoldChange_threshold < 1, "female", "male") 
padj_threshold = ifelse(results$padj_threshold < 0.05, "female", "male") 
results$log2FoldChange_threshold > results$log2FoldChange 
# Tried indexing to see only the values where log2FoldChange_threshold is greater than log2FoldChange
results$padj_threshold > results$padj 
results = results [ results$padj_threshold > results$padj  , results$log2FoldChange_threshold > results$log2FoldChange ]



fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

plot(x = log2FoldChange_threshold,
     y = -log10(padj_threshold),
     xlab = "-log10(p value)", # be sure the specify that it's female over male!
     ylab = "Log2 Fold Change",
     pch = 20) # smaller solid circles

library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange_threshold, y = -log10(padj_threshold))) + 
  geom_point(aes(color = ifelse(log2FoldChange_threshold < -1 & padj_threshold < 0.05, "lower in female",
                                ifelse(log2FoldChange_threshold > 1 & padj_threshold < 0.05, "higher in female", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Female/Male)",
       y = "-log10 Adjusted p-value")


volcano_plot







