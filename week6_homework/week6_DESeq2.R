#1.1
#BiocManager::install("DESeq2")
library(SummarizedExperiment)
library(DESeq2)
library(TCGAbiolinks)

getwd()
setwd("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data")
getwd()

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)

#1.2
patient_data=colData(sum_exp)
counts = assays(sum_exp)$"HTSeq - Counts"
counts = counts[ ,!is.na(patient_data$age_at_index)]

patient_data=patient_data[!is.na(patient_data$age_at_index),]
patient_data$age_category=ifelse(patient_data$age_at_index<50, "young", "old")
patient_data$age_category = factor(patient_data$age_category, levels = c("young", "old"))      

#1.3
counts_row_sums=rowSums(counts)
low_counts_mask=ifelse(counts_row_sums<10, FALSE, TRUE)
low_counts_mask
sum(!low_counts_mask) #count how many <10 there are
counts=counts[low_counts_mask,]

#2.1
dds=DESeqDataSetFromMatrix(countData=counts, colData=patient_data, design=~age_category)
dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

#2.2
head(results)
str(results)
my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))

order_indices = order(my_df$y)
# we expect c(5, 1, 3, 2, 4) because:
# 1 is index 5
# 2 is index 1
# 3 is index 3
# 4 is index 2
# 5 is index 4
order_indices  # note the order!

my_df = my_df[order_indices, ]
my_df

#2.3
row_order=order(results)
head(results, 20)
#I chose the gene ENSG00000001626. It is expressed more in older patients as it
#has a negative log2FoldChange. It is known as CFTR, or CF Transmembrane
#Conductance Regulator. This gene codes for a protein that functions as a
#chloride channel for controlling ion and water secretion and absorption in
#epithelial tissues.

log2FoldChange_threshold = ifelse(results$log2FoldChange<1, "young", "old") 
padj_threshold = ifelse(results$padj<0.05, "young", "old") 
#I think I messed up my operators somehow here
mask1=results$log2FoldChange>results$log2FoldChange_threshold || results$log2FoldChange<(-1*results$log2FoldChange_threshold)
mask2=results$padj < results$padj_threshold
results = results [mask1 &&mask2,]

#2.5
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

plot(x = results$log2FoldChange,
     y = -log10(results$padj),
     xlab = "Log2 Fold Change)", # be sure the specify that it's young over old!
     ylab = "-log10(p value)",
     pch = 20) # smaller solid circles

abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

#2.6
write.csv(x = results,
          file = "/Users/peterchong/Documents/qbio_data_analysis_peterc/week6_homework/results.csv",
          row.names = FALSE)





