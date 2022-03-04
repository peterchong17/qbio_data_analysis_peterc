# Exercise 1.1
BiocManager::install("maftools")
library(TCGAbiolinks)

getwd()
setwd("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data")
getwd()

# Exercise 1.2
clinic <- data.table::fread("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data/coad_clinical_data.csv",
                            data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

# Exercise 1.3
length(colnames(clinic))
# colnames(clinic)is 78 names long.
length(clinic[, "Tumor_Sample_Barcode"])
# The newly renamed column into Tumor Sample Barcode is 524 strings long.
sum(colnames(clinic) == "bcr_patient_barcode")
# Before the change, there was only one column name that matched "bcr_patient_barcode." Now, there are 0 as it got renamed to "Tumor_Sample_Barcode."

# Exercise 1.4
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
# Exercise 1.5
getwd()
setwd("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data/GDCdata/")
list.files()
maf_dataframe = data.table::fread("TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", data.table=F)
clinic = data.table::fread("/Users/peterchong/Documents/qbio_data_analysis_peterc/analysis_data/coad_clinical_data.csv",
                         data.table = F)
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] = "Tumor_Sample_Barcode"
maf_object = read.maf(maf = maf_dataframe, clinicalData = clinic, isTCGA = T)

# Exercise 2.1
maf_object
str(maf_object)
maf_object@data
str(maf_object@data)
maf_object@clinical.data
str(maf_object@clinical.data)

colnames(maf_object@data)
colnames(maf_object@clinical.data)
# the shared column is Tumor_Sample_Barcode
# this does make sense as clinical data requires recording the barcode with with each patient is labelled with. Similarly, data about the mutation also needs to reference in
      # in which patient the mutation occurred in.

# Exercise 3.1
oncoplot(maf = maf_object,
         top = 19) 

ggsave("/Users/peterchong/Documents/qbio_data_analysis_peterc/week7_maf/oncoplot.png") 

# Exercise 3.2
# The TP53 gene codes for tumor protein p53.
# This protein, and therefore the gene, regulates cell division by keeping it from growing/dividing too fast and uncontrollably.


# Exercise 3.3
clinic = maf_object@clinical.data
young_patient_ids = clinic$Tumor_Sample_Barcode[clinic$age_category == "young"]
young_maf = subsetMaf(maf = maf_object, young_patient_ids)

old_maf = subsetMaf(maf = maf_object, clinic$Tumor_Sample_Barcode[clinic$age_category == "old"])

# Exercise 3.4
coOncoplot(m1 = young_maf,
           m2 = old_maf,
           m1Name = "Young Patients         ",
           m2Name = "Old Patients")
ggsave("/Users/peterchong/Documents/qbio_data_analysis_peterc/week7_maf/coOncoplot.png")

# Exercise 3.5
# Not all of the genes are more mutated in one more group than the other
# 4 of 6 are more mutated in old patients, while 2 of them are more mutated in young patients.
# I expected old patients to have more mutations across the board as mutations accrue the older an individual is. 

# Exercise 3.6
lollipopPlot(maf_object, gene = "TP53")
ggsave("/Users/peterchong/Documents/qbio_data_analysis_peterc/week7_maf/lollipopTP53.png")

# Exercise 3.7
lollipopPlot2(m1 = young_maf,
           m2 = old_maf,
           m1_name = "Young Patients         ",
           m2_name = "Old Patients",
           gene = "TP53")
# The gene is more mutated in old patients and most of the mutations occur in the P53 site, typically missense mutations.

ggsave("/Users/peterchong/Documents/qbio_data_analysis_peterc/week7_maf/lollipopPlot2TP53.png")

# Exercise 3.8
# There are 10 samples that have neither a mutation in Gene A or Gene B. This is a 10% of samples.

# Exercise 3.9
# b = 7, c = 2, d = 35, e = 37, f = 42
# 6   7   13
# 2   35  37
# 8   42  50
# The mutations in gene A and gene B are not independent.
# Knowing that there is a mutation in gene B significantly increases the chance that there is a mutation in gene A.
# This is proven by the fact that 6 of 8 samples with a mutation in gene B have a mutation in gene A as well.

# Exercise 3.10
# geneA is TP53
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

# geneB is KRAS
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")

# Exercise 3.11
geneA_maf
geneB_maf
# subsetMaf() creates smaller Maf dataframes that contain a specific variable that we want.
maf_object@data
maf_object@clinical.data
# I do not think that that there is only one mutation per sample in geneA.
# This is because there are many many more mutations than patients, as can be seen in the vastly different
# number of rows in data versus clinical.data
geneA_maf@data #222
geneA_maf@clinical.data # 213
geneB_maf@data #165
geneB_maf@clinical.data # 163
# No the number of samples in each data section is not the same as the number of samples in the clinical.data section.
# This is because some patients have more than one mutation in their genes; thus they pop more more than once in data.

# Exercise 3.12
# 1. Access the barcodes of the patients with mutations in genes A and B
# bc stands for barcode
mut_bc_geneA = geneA_maf@clinical.data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@clinical.data$Tumor_Sample_Barcode

# 2. Get the lengths of these two vectors
num_mut_geneA = length(mut_bc_geneA)
num_mut_geneA
num_mut_geneB = length(mut_bc_geneB)
num_mut_geneB

# 3. Fill in the intersect here! Then get the nubmer of patients
mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)
num_mut_geneAB

# Exercise 3.13
num_mut_geneA_only = 135
num_mut_geneB_only = 85

# Exercise 3.14
num_neither_mutation = 78

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)

contig_table

# Exercise 3.15
fe_results <- fisher.test(contig_table)
fe_results
# The p-value is extremely low, therefore we reject the null hypothesis, saying that it is very likely
# that the mutations are not independent and therefore dependent.