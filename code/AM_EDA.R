#install.packages("caret")
#install.packages("psych")
#install.packages("Hmisc")
#install.packages("readr")
#install.packages("stringr")

#library(Hmisc)
#library(caret)
#library(psych)
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

patients <- read_csv("patients.csv")
microarray <- read_csv("microarray.csv")


#Replacing spaces and dashes by underscores
colnames(patients) <- str_replace_all(colnames(patients), " ", "_")
colnames(patients) <- str_replace_all(colnames(patients), "-", "_")
colnames(patients)
#We can link these two data sets by comparing the DLBCL sample (LYM number) in the “patients” 
#with the LYM number in column names of the “microarray” data.

#Creating a modified dataset to enhance interpretability and linkability
microarray_LYM <- microarray

#Renaming the columns for the corresponding LYM number
for (i in 3:276){
  current_name = colnames(microarray_LYM)[i]
  first_pos_of_LYM = str_locate(current_name, "LYM")[1]
  colnames(microarray_LYM)[i] = substring(current_name, first_pos_of_LYM+3, first_pos_of_LYM+5)
}

#merging of the two tables (each row of the table micro_array table corresponds to a sample
#in order to merge it, it is necessary to transpose the table. The name of each variable is the
#concatenation of the first two rows of the table)
microarray_LYM$concatenated <- paste(microarray_LYM$UNIQID, "-", microarray_LYM$NAME) #creating the column corresponding to the name of each variable of the sample
microarray_LYM_transposed <- data.frame(t(microarray_LYM[ , 3:(ncol(microarray_LYM)-1)])) #transposing of the table in order to have 293 observations with 7293 variables
colnames(microarray_LYM_transposed) <- microarray_LYM$concatenated #definition of the name of each variable
microarray_LYM_transposed$`DLBCL_sample_(LYM_number)` <- rownames(microarray_LYM_transposed)#definition of each patients Lym_number as a variable to be later used on the merging

final_merged_table = merge(patients, microarray_LYM_transposed, by = "DLBCL_sample_(LYM_number)") #inner join of the patients and the microarray tables

#missing values
sum(is.na(final_merged_table))
#Distribution of Predicted Outcome by Group of Lymphoma
#ggplot(patients, aes(Outcome_predictor_score, colour = Subgroup)) +
       geom_freqpoly(binwidth = 1) + labs(title="-")


###########################################EDA#########################
#bar charts for categorical

ggplot(patients, aes(x=Analysis_Set))+geom_bar(fill='lightskyblue4')+labs(x='Analysis set')
ggplot(patients, aes(x=`Status_at_follow-up`))+geom_bar(fill='lightskyblue4')+labs(x='Status at follow-up')
ggplot(patients, aes(x=patients$Subgroup))+geom_bar(fill='lightskyblue4')+labs(x='Subgroup')
ggplot(patients, aes(x=patients$IPI_Group))+geom_bar(fill='lightskyblue4')+labs(x='IPI group')

require(gridExtra)
grid.arrange(ggplot(patients, aes(x=Analysis_Set))+geom_bar(fill='lightskyblue4')+labs(x='Analysis set'),
             ggplot(patients, aes(x=`Status_at_follow-up`))+geom_bar(fill='lightskyblue4')+labs(x='Status at follow-up'),
             ggplot(patients, aes(x=patients$Subgroup))+geom_bar(fill='lightskyblue4')+labs(x='Subgroup'),
             ggplot(patients, aes(x=patients$IPI_Group))+geom_bar(fill='lightskyblue4')+labs(x='IPI group'), ncol=2)



#ggplot outcome predictor com density 
ggplot(patients, aes(x=Outcome_predictor_score)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="blue")+
  geom_density(alpha=.5, colour="red")

#histograms for continuous independent variables

ggplot(patients, aes(x=patients$`Follow-up_(years)`))  +labs(x='Follow-up(years)')+
  geom_histogram(aes(y=..density..), colour="black", fill="grey")

ggplot(patients, aes(x=patients$Germinal_center_B_cell_signature)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey")

ggplot(patients, aes(x=patients$Lymph_node_signature)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey")

ggplot(patients, aes(x=patients$Proliferation_signature)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey")

ggplot(patients, aes(x=patients$BMP6)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey")

ggplot(patients, aes(x=patients$MHC_class_II_signature)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey")

require(gridExtra)
grid.arrange(ggplot(patients, aes(x=patients$`Follow-up_(years)`)) + labs(x='Follow-up(years)')+
               geom_histogram(aes(y=..density..), colour="black", fill="grey"),
             ggplot(patients, aes(x=patients$Germinal_center_B_cell_signature)) +labs(x='Germinal center B cell signature')+
               geom_histogram(aes(y=..density..), colour="black", fill="grey"),
             ggplot(patients, aes(x=patients$Lymph_node_signature)) + labs(x='Lymph node signature')+
               geom_histogram(aes(y=..density..), colour="black", fill="grey"),
             ggplot(patients, aes(x=patients$Proliferation_signature)) + labs(x='Proliferation signature')+
               geom_histogram(aes(y=..density..), colour="black", fill="grey"),
             ggplot(patients, aes(x=patients$BMP6)) + labs(x='BMP6')+labs(x='BMP6')+
               geom_histogram(aes(y=..density..), colour="black", fill="grey"),
             ggplot(patients, aes(x=patients$MHC_class_II_signature)) + labs(x='MHC class II signature')+
               geom_histogram(aes(y=..density..), colour="black", fill="grey"))

#ggplots em relacao com a target - categoricas

a <- ggplot(patients, aes(x=Outcome_predictor_score, fill=IPI_Group, color=IPI_Group)) +
  geom_histogram(binwidth = 0.1) + labs(title="-")
a + theme_bw()

b <- ggplot(patients, aes(x=Outcome_predictor_score, fill=Subgroup, color=Subgroup)) +
  geom_histogram(binwidth = 0.1) + labs(title="-")
b + theme_bw()
c <- ggplot(patients, aes(x=Outcome_predictor_score, fill=Status_at_follow_up, color=Status_at_follow_up)) +
  geom_histogram(binwidth = 0.1) + labs(title="-")
c + theme_bw()

#ggplots em relacao com a target - numericas

ggplot(patients,x=aes(BMP6, color=patients$Outcome_predictor_score))

#colnames(patients)
patients[, c('Analysis_Set')] <- list(NULL)
patients

#correlation chart
library(PerformanceAnalytics)
chart.Correlation(patients[,c(3,7,8,9,10,11,12)], histogram=TRUE, pch=19)

data.cpca <- prcomp(patients[,1:9], scale. = FALSE, retx=TRUE)
print(data.cpca)

#report
library(DataExplorer)  #tabela de correlacao etc
create_report(patients)

par(mfrow=c(2,2))
#boxplots
boxplot(Outcome_predictor_score~Status_at_follow_up,
        data=final_merged_table,
        col="floralwhite",
        border="sienna"
)

boxplot(Outcome_predictor_score ~Subgroup,
        data=final_merged_table,
        col="floralwhite",
        border="sienna"
)
boxplot(Outcome_predictor_score~IPI_Group,
        data=final_merged_table,
        col="floralwhite",
        border="sienna"
)

# Kolmogorov-Smirnov normality test (since n>=50)
#nrow(patients)
My_list <- split(patients, f = list(colnames(patients)))
loop_Shapiro2 <- lapply(My_list, function(x) ks.test(x$Outcome_predictor_score, "pnorm"))
print(loop_Shapiro2)

final_merged_table

#pca
results <- prcomp(final_merged_table, scale = TRUE)

#reverse the signs
results$rotation <- -1*results$rotation

#display principal components
results$rotation


final_merged_table[, c("Subgroup",'Analysis_Set','Status_at_follow_up',"IPI_Group")] <- list(NULL)
final_merged_table
str(final_merged_table)
str(final_merged_table)
cor(final_merged_table)

library(tibble)
library(tidyr)



