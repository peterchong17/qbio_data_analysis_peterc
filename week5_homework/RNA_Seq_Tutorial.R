#1.1
BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
library(SummarizedExperiment)
getwd()
setwd("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data")
getwd()

#2.1
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)

#2.2
counts = assays(sum_exp)$"HTSeq - Counts"
assays(sum_exp)$"HTSeq - Counts" [1:5, 1:5]
head(rowData(sum_exp))
colData(sum_exp)[1:5, 25:29]
metadata(sum_exp)

#2.3
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")

#2.4
str(colData(sum_exp))
head(colData(sum_exp))

#2.5
colnames(colData(sum_exp))
#age_at_diagnosis

#2.6
colData(sum_exp)$age_at_diagnosis[1:10]

#2.7
colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis/365

#2.8
colData(sum_exp)$age_category=ifelse(colData(sum_exp)$age_at_diagnosis<50,"Young", "Old")
colData(sum_exp)$age_category=ifelse(is.na(colData(sum_exp)$age_at_diagnosis), NA, colData(sum_exp)$age_category)
colData(sum_exp)$age_category

#2.9
head(rowData(sum_exp))
dim(rowData(sum_exp))

#2.10
"MSH2" %in% rowData(sum_exp)$external_gene_name
"MSH6" %in% rowData(sum_exp)$external_gene_name

#2.11
assays(sum_exp)$"HTSeq - Counts"[20:25,]
assays(sum_exp)$"HTSeq - Counts"[,30-35]

#2.12
mask1=rowData(sum_exp)$external_gene_name=="MSH2"
sum(mask1)
ensembl_geneA=rowData(sum_exp)[mask1,]$ensembl_gene_id

mask2=rowData(sum_exp)$external_gene_name=="MSH6"
sum(mask2)
ensembl_geneB=rowData(sum_exp)[mask2,]$ensembl_gene_id

#2.13
#Ensembl gene ID is a column in assays(sum_exp)"HTSeq - Counts"

#2.14
min(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA , ])
max(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA , ])
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB , ])

assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA , ]

#2.15 couldn't get this to work :'(
plot(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ],
     xlab = "Gene A",
     ylab = "Gene B"
)

#2.16
bool_age_na=is.na(colData(sum_exp)$age_category)
num_na=sum(bool_age_na)
num_na
# Removing the NAs would create data with unequal lengths in the counts data.

#2.17
age_cat_no_NAs=colData(sum_exp)[!bool_age_na,]
sum(is.na(age_cat_no_NAs$age_category))

#2.18
length(age_cat_no_NAs)
dim(colData(sum_exp))
dim(colData(sum_exp))[1]
dim(colData(sum_exp))[2]

#2.19
dim(assays(sum_exp)$"HTSeq - Counts")
dim(age_cat_no_NAs)
#no the number of patient age categories are not the same, because

#2.20
identical(colnames(assays(sum_exp)$"HTSeq - Counts"), rownames(colData(sum_exp)))
assays(sum_exp)$"HTSeq - Counts"[,!bool_age_na]
assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,]
gene_counts=assays(sum_exp)$"HTSeq - Counts"[,!bool_age_na]

#2.21
length(age_cat_no_NAs)==length(gene_counts)
# I'm not entirely sure what the variable gene_counts is supposed to be.
# I had some difficulty with this section, I think because I'm not sure how to access the number of counts from which dataframe and from which column.

#2.22
# Because I didn't know how to make the gene_counts variable, I can't figure out how to make the boxplot.

#3
#1. You access it through assays(sum_exp)$"HTSeq - Counts". The columns are the patient IDs and the rows are the Ensembl Gene ID.
#2. You can access it through rowData(sum_exp). It stores the gene id, its common gene name, as well as its gene id. The rows of counts holds the gene id.
#3. You can access it through colData(sum_exp). It stores various clinical data, with the rows being each patient's id. Counts holds patient ids in its columns.
  
  

