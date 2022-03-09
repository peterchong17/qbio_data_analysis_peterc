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


# BOX PLOT
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
jpeg("/Users/peterchong/Documents/qbio_data_analysis_peterc/Midterm Project/gender_boxplot.jpg") # saves the volcano plot
box=boxplot(gene_counts ~ gender_no_NAs, xlab = "Gender", ylab = "Counts of CFTR")
box$stats # summarizes the stats for the two boxplots to find the medians. For female, it is 5376 and for male it is 5742
dev.off()

#clean counts of the samples missing gender information
counts=assays(sum_exp)$"HTSeq - Counts"
patient_data=colData(sum_exp)
patient_data_na_mask = is.na(colData(sum_exp)$gender)
counts = counts[ , !(patient_data_na_mask) ]
patient_data = patient_data[ !(patient_data_na_mask) , ]

if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("row names changed")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
counts_row_sums = rowSums(counts)
low_counts_mask = counts_row_sums < 10
counts = counts[ !low_counts_mask, ]
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = patient_data, 
                             design = ~gender)

patient_data$gender=factor(patient_data$gender, levels=c("female", "male"))



# VOLCANO PLOT

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run
  
# get the female vs. male comparison
results = results(dds_obj, format = "DataFrame", contrast = c("gender", "female", "male"))

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05
log2FoldChange_threshold = log2(fc_threshold) # sets a threshold by calculating log2 of the fc_threshold

padj_mask = results$padj < p_threshold # creates a boolean mask of the results that are significant (padj < 0.05)
log2FoldChange_mask_up = results$log2FoldChange > log2FoldChange_threshold # creates a boolean mask of the results that are overexpressed in females; underexpressed in males
log2FoldChange_mask_down = results$log2FoldChange < -log2FoldChange_threshold # creates a boolean mask of the results that are underexpressed in females; overexpressed in males
namask = !is.na(results$padj)

sub_results_up = results [ padj_mask & log2FoldChange_mask_up & namask, ] # creates a table of the results that are overexpressed (in females)
sub_results_down = results [ padj_mask & log2FoldChange_mask_down & namask, ] # creates a table of the results that are underexpressed
write.csv(sub_results_up, "/Users/peterchong/Documents/qbio_data_analysis_peterc/Midterm Project/sub_results_up.csv") # saves the overexpressed table
write.csv(sub_results_down, "/Users/peterchong/Documents/qbio_data_analysis_peterc/Midterm Project/sub_results_down.csv") # saves the underexpressed table

jpeg("/Users/peterchong/Documents/qbio_data_analysis_peterc/Midterm Project/gender_volplot.jpg") # saves the volcano plot
plot(x = results$log2FoldChange, y= -log10(results$padj), xlab = "log2(FoldChange)", ylab = "-log10(padj)")
  abline(v=c(log2FoldChange_threshold, -log2FoldChange_threshold), h = c(-log10(p_threshold)), col = "green")
  
dev.off()



# KAPLAN MEIER OF EXPRESSION

ensembl_CFTR=rowData(sum_exp)[rowData(sum_exp)$external_gene_name=="CFTR",]$ensembl_gene_id # finds the gene ensemble name of CFTR
counts_CFTR=counts[ensembl_CFTR,] # creates a subset of counts that only contains CFTR info
summary(counts_CFTR) # found that the median count value is 5707
patient_data=colData(sum_exp)[!bool_gender_na, ]

#patient_data$epression=ifelse(counts_CFTR>5707, "overexpressed", 'underexpressed') # create new column of over or underexpressed
patient_data$expression[counts_CFTR>10772] = "overexpressed"
patient_data$expression[counts_CFTR<2396] = "underexpresed"
patient_data$expression[counts_CFTR<10772 & counts_CFTR>2396] = "normal"

# if there is no days_to_death information, it uses days_to_last_follow_up as an estimate
patient_data$days_to_death = ifelse(is.na(patient_data$days_to_death), patient_data$days_to_last_follow_up, patient_data$days_to_death)
sum(is.na(patient_data$days_to_death)) #checking that there is no NAs in days_to_death
patient_data = patient_data[!is.na(patient_data$days_to_death),] # removes rows in patient data that do not have info

# create new column death event by turning the numerical vital_status into boolean values
patient_data$death_event = ifelse(patient_data$vital_status=="Dead", T, F)

# initialize a survival object for the time frame that we want
surv_object <- Surv(time = patient_data$days_to_death, 
                    event = patient_data$death_event)

# initialize a fit object for over vs underexpression 
cftr_fit <- surv_fit( surv_object ~ patient_data$expression, data = patient_data )

# create the age survival plot, while inputting proper margins and putting the legend to the right
survplot = ggsurvplot(cftr_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# tweaking the formatting of the age plot
expr_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
# calling and showing the age survival plot
jpeg("/Users/peterchong/Documents/qbio_data_analysis_peterc/Midterm Project/expression_kapmeier.jpg") # saves the expression kaplan meier plot
expr_surv_plot
dev.off()



# KAPLAN MEIER OF GENDER

# initialize a fit object for age 
gender_fit <- surv_fit( surv_object ~ patient_data$gender, data = patient_data )

# create the age survival plot, while inputting proper margins and putting the legend to the right
survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# tweaking the formatting of the age plot
gender_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
# calling and showing the age survival plot
jpeg("/Users/peterchong/Documents/qbio_data_analysis_peterc/Midterm Project/gender_kapmeier.jpg") # saves the gender kaplan meier plot
gender_surv_plot
dev.off()
