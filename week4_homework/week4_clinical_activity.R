#1. A categorical variable is a variable whose observations are categorized into groups. An example is favorite colors.
# Discrete variables have data whos numeric values are discrete (i.e. whole numbers). An example is population.
# A continuous variable has data that is continuous (can be any number) within a range. An example is height.

#2.
clinic <- read.csv("/Users/peterchong/Documents/qbio_data_analysis_peterc/week4_homework/coad_clinical_data.csv")
names(clinic)
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))
# we chose age as our variable

#3.
#Age is measured by collecting personal information from the patient. The variable is continuous (but the data provided in clinic is discrete).

#4.
#"Age and Cancer Risk" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544764/
# This article talks about how to modify age-related "symptoms" in order to address cancer progression.
#"The Biology of Aging and Cancer: A Brief Overview of Shared and Divergent Molecular Hallmarks" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5614326/
# This article talks about how age and cancer are linked. It looks at specific markers of aging and how they affect cancer, such as telomere attrition.

#5.
sum(is.na(clinic$height))
# we chose whether the patients' height as our  second variable
# Height is the distane length of a person from head to toe. It is measured via rulers.

#6.
# 1) The older a patient is, the shorter the average patient becomes.
# 2) Older age is correlated with a lower survival rate in cancer patients.
# 3) Shorter people are more likely to have decreases survival rates in cancer patients.

#7.
# From the first plot relating age and height, we see that young people are on average taller than older people,
# but also that older people have a larger spread of heights. From the survival plot of age, young people tend 
# have overall higher survival. From the height survival plot, short people tend to live longer after the beginning
# period of time of colorectal cancer.

#CODING SEGMENT
#1.
# It is easier to see the trend if we use age categories rather than the age directly.
# Thus, the continuous variable age turns into the categorical variable age_category 

# create a boolean mask where NAs in height are marked as FALSE and everything else is TRUE
mask=!is.na(clinic$height)
# create cleaned dataframe where only rows with heights are included
cleaned_clinic <- clinic[mask,]

#create a boxplot where the height ranges of the age categories of young and old (<50 and >=50) are graphed
boxplot(cleaned_clinic$height~cleaned_clinic$age_category,
        xlab="age",
        ylab='height', main="Age vs Height", ylim=c(120,200))

library(survival)
library(survminer)

#2.
# initialize a survival object for the time frame that we want
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# initialize a fit object for age 
age_fit <- surv_fit( surv_object ~ clinic$age_category, data = clinic )

# create the age survival plot, while inputting proper margins and putting the legend to the right
survplot = ggsurvplot(age_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# tweaking the formatting of the age plot
age_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
# calling and showing the age survival plot
age_surv_plot


#3.

# finding the center of the cleaned clinic data height in order to define tall as above average and short as below average
median(cleaned_clinic$height)

# initializing a column height_category to turn height into a categorical variable so that it can be put into a survival plot
cleaned_clinic$height_category = ifelse(cleaned_clinic$height<=170, "short", "tall")

# initialize a survival object for the time frame that we want
surv_object <- Surv(time = cleaned_clinic$days_to_death, 
                    event = cleaned_clinic$death_event)

# initialize a fit object for height 
height_fit <- surv_fit( surv_object ~ cleaned_clinic$height_category, data = cleaned_clinic )

# create the height survival plot, while inputting proper margins and putting the legend to the right
survplot = ggsurvplot(height_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# tweaking the formatting of the height plot
height_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
# calling and showing the height survival plot
height_surv_plot
