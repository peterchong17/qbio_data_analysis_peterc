if(!require(BiocManager)) install.packages("BiocManager")
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
getwd()
setwd("Documents/qbio_data_analysis_peterc/analysis_data")
getwd()
clin_query <- GDCquery(project = "TCGA-COAD", 
                       data.category = "Clinical",
                       file.type = "xml")
#GDCdownload(clin_query) 
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
names(clinic)[names(clinic)=="days_to_last_followup"]<-"days_to_last_follow_up"

# exercise 1.1
str(clinic)
head(clinic)
#str() looks at the overall structure, while head() gives you the first 6 lines. I think the days_to_birth variable is interesting in that they chose to use this format to say how old the patients are.

# exercise 1.2
colnames(clinic)
clinic$vital_status

# exercise 2.1
plot(x=clinic$age_at_initial_pathologic_diagnosis,y=clinic$weight,
     xlab="Initial Diagnosis Age", ylab="Weight")

# exercise 2.2
unique(clinic$race_list)
par(mar=c(10,1,1,1))
boxplot(clinic$age_at_initial_pathologic_diagnosis~clinic$race_list, las=2, cex.axis=0.5, xlab="Race", ylab="Age")

# exercise 2.3
sum(clinic$race_list=="")
clinic$race_list[clinic$race_list==""]="No data"

clinic$race_list = as.character(clinic$race_list)
clinic$race_list

# exercise 2.4
min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
summary(clinic$age_at_initial_pathologic_diagnosis)

# exercise 2.5
numYoung=sum(clinic$age_at_initial_pathologic_diagnosis<50)
numOld=sum(clinic$age_at_initial_pathologic_diagnosis>=50)
numYoung+numOld==nrow(clinic)

# exercise 2.6
boolAge=clinic$age_at_initial_pathologic_diagnosis<50
young_patient_ids=clinic$bcr_patient_uuid[clinic$age_at_initial_pathologic_diagnosis[boolAge]]
old_patient_ids=clinic$bcr_patient_uuid[clinic$age_at_initial_pathologic_diagnosis[!boolAge]]
length(young_patient_ids)==numYoung
length(old_patient_ids)==numOld

# exercise 2.7
clinic$age_category=ifelse(clinic$age_at_initial_pathologic_diagnosis<50, "young", "old")
head(clinic)

# exercise 2.8
clinic[1,1] #top left entry of dataframe
clinic[1,] #first row of dataframe
clinic[2:5,] #rows 2 through 5 of the dataframe
clinic[,3] #3rd column

# exercise 2.9
young_clinic=clinic[clinic$age_category=="young",]
old_clinic=clinic[clinic$age_category=="old",]

# exercise 2.10
young_clinic_one_line=clinic[clinic$age_at_initial_pathologic_diagnosis<50,]
head(young_clinic_one_line)
identical(dim(young_clinic), dim(young_clinic_one_line))

# exercises 3
install.packages("survival")
install.packages("survminer")
library(survival)
library(survminer)

# exercise 3.1
clinic$days_to_death[is.na(clinic$days_to_death)] <- clinic$days_to_last_follow_up[is.na(clinic$days_to_death)] 
head(clinic)

# exercise 3.2
clinic$death_event=ifelse(clinic$vital_status=="Alive", 0, 1)
head(clinic)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("../week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)

# exercise 3.3
# African American patients have a decreased probability of surviving over time than whites, with Asians marked as a 100% survival probability
# However, the time allotted for the Asian data points are not very long and indicate that it is not finished.)
# Additionally, there seems to be no or very limited data on American Indian or Alaska Natives.
# One question is why do African Americans have lower survival rates than the other races? What biological reasons form this?

write.csv(clinic, "/Users/peterchong/Documents/qbio_data_analysis_peterc/week4_homework/coad_clinical_data.csv", row.names = F)
