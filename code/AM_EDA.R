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

patients <- read_csv("patients.csv")
microarray <- read_csv("microarray.csv")

#Replacing spaces and dashes by underscores
colnames(patients) <- str_replace_all(colnames(patients), " ", "_")
colnames(patients) <- str_replace_all(colnames(patients), "-", "_")

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

#Distribution of Predicted Outcome by Group of Lymphoma
ggplot(patients, aes(Outcome_predictor_score, colour = Subgroup)) +
       geom_freqpoly(binwidth = 1) + labs(title="-")

ggplot(patients, aes(Outcome_predictor_score, colour = Status_at_follow_up)) +
  geom_freqpoly(binwidth = 1) + labs(title="-")

c <- ggplot(patients, aes(x=Outcome_predictor_score, fill=Status_at_follow_up, color=Status_at_follow_up)) +
  geom_histogram(binwidth = 0.1) + labs(title="-")
c + theme_bw()




