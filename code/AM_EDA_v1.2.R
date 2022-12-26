#library(Hmisc)
#library(psych)

{
  library(readr)
  library(stringr)
  library(ggplot2)
  library(dplyr)
  library(caret)
  library(tibble)
  library(tidyr)
  library(tseries)
  library(rrcov)
  
}

patients <- read_csv("patients.csv")
microarray <- read_csv("microarray.csv")

#Original Datasets available at https://llmpp.nih.gov/DLBCL/
#patients2 <- read.csv("/Users/bernardo/Desktop/MMAC/3sem/AM/DLBCL_patient_data_NEW.csv", sep="\t")
#microarray2 <- read.csv("/Users/bernardo/Desktop/MMAC/3sem/AM/NEJM_Web_Fig1data.csv", sep="\t")

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
  colnames(microarray_LYM)[i] = str_remove(substring(current_name, first_pos_of_LYM+3, first_pos_of_LYM+5),"^0+")
}

#Storing those weird gene names column in a separate dataframe and eliminate that column from the following manipulations
gene_database <- microarray_LYM[1:2]
microarray_LYM <- subset(microarray_LYM, select = -c(2)) 

#merging of the two tables (each row of the table micro_array table corresponds to a sample
#in order to merge it, it is necessary to transpose the table. The name of each variable is the
#concatenation of the first two rows of the table)
microarray_LYM_transposed <- data.frame(t(microarray_LYM[ , 2:(ncol(microarray_LYM))])) #transposing of the table in order to have 292 observations with 7292 variables
colnames(microarray_LYM_transposed) <- microarray_LYM$UNIQID #definition of the name of each variable as the unique ID of a gene
microarray_LYM_transposed$`DLBCL_sample_(LYM_number)` <- rownames(microarray_LYM_transposed)#definition of each patients Lym_number as a variable to be later used on the merging

merged_data_v1 = merge(patients, microarray_LYM_transposed, by = "DLBCL_sample_(LYM_number)") #inner join of the patients and the microarray tables

#missing values
sum(is.na(merged_data_v1))
#We have 177352 NAs

#LYM NUMBERS EM FALTA ANTES:
d <- c(1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,11 ,12 ,13 , 14, 15, 16, 17, 18, 19, 22, 24, 25, 26, 27, 30, 32, 34, 35, 38,
       + 40, 41, 43, 44, 45, 48, 50, 51, 54, 57, 58, 59, 61, 62, 63, 64, 65, 66, 68, 69, 72, 74, 76, 78, 79, 80, 81, 83,
       + 86, 88, 89, 91, 92, 93, 94, 95, 97, 98, 99)

#Distribution of the % of Missing values by number of features
hist(colMeans(is.na(merged_data_v1))*100,  col = 1:100, breaks=c(0:100), 
     xlab = '% of NAs', ylab= 'Number of Features', main='') 
#Distribution of the % of Missing values by number of entries
hist(rowMeans(is.na(merged_data_v1))*100,  col = 1:100, breaks=c(0:100), 
     xlab = '% of NAs', ylab= 'Number of Entries', main='') 

#From the first 12 features of the patients dataset, only the "IPI_group" has missing values and has few (6.936416 %)
#Thus, we want to focus on the features that correspond to specific genes (the features 13 to 7303) which have a lot more NAs

#We have 7303-12=7291 genes in this dataset and in the microarray dataset we have 7291 genes --- nice, checks out!
#Mean of % of NAs in genes features
mean(colMeans(is.na(merged_data_v1[13:7303])))

#Average of maximum values in genes features
mean(mapply(max,merged_data_v1[13:7303], na.rm = TRUE))
#Average of minimum values in genes features
mean(mapply(min,merged_data_v1[13:7303], na.rm = TRUE))

#Maximum Expression Level over all gene features
max(merged_data_v1[13:7303],na.rm = TRUE)
#Minimum Expression Level over all gene features
min(merged_data_v1[13:7303],na.rm = TRUE)

#Eliminate features with more than 20% of NAs
per_NA <- colMeans(is.na(merged_data_v1))
merged_data_v2 <- merged_data_v1[-which(per_NA >= 0.1)]

#With this cut we end up with 4949 features

#Hence, we eliminated 7303-4949=2354 gene features that had 10% of more of NAs
sum(is.na(merged_data_v2))
#Now we have 32375 NAs (we reduced (177352-32375)*100/177352 = 81.75% of NAs)


#Removing highly correlated gene features
#https://stats.stackexchange.com/questions/50537/should-one-remove-highly-correlated-variables-before-doing-pca
{
  patients_features <- merged_data_v2[1:12]
  genes_features <- merged_data_v2[13:4949]
  #Creating matrix of correlations
  #aux = cor(genes_data, use = "pairwise.complete.obs")
  hc = findCorrelation(cor(genes_features, use = "pairwise.complete.obs"), cutoff=0.95) #Eliminate features with more than 0.95 correlation with others
  print(length(hc))
  hc = sort(hc)
  reduced_genes_features = genes_features[,-c(hc)]
  merged_data_v3 <- cbind(patients_features, reduced_genes_features)
  }
#With this reduction of gene features, we end up with 4578 gene features, and 4590 features in total
#Still a lot of features but we started with 7303!

#Checking how many NAs at this point
sum(is.na(merged_data_v3))
#With this version we have 31669 NAs (17.85% of the initial amount)



#Replacing the Nas with their median value in collumn

merged_data_v4 <- merged_data_v3 %>% mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE))))

sum(is.na(merged_data_v4))
sum(is.na(patients_features))

#We replaced every Na except for the IPI subgroup hence the is.numeric command

genes_final = merged_data_v4[13:4590]

sum(is.na(genes_final))

########################################### EDA #########################

#Distribution of Predicted Outcome by Group of Lymphoma
ggplot(patients, aes(Outcome_predictor_score, colour = Subgroup)) +
  geom_freqpoly(binwidth = 1) + labs(title="-")




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

ggplot(patients,x=aes(BMP6, color=patients$Outcome_predictor_score))4


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
        data=merged_data_v3,
        col="floralwhite",
        border="sienna"
)

boxplot(Outcome_predictor_score ~Subgroup,
        data=merged_data_v3,
        col="floralwhite",
        border="sienna"
)
boxplot(Outcome_predictor_score~IPI_Group,
        data=merged_data_v3,
        col="floralwhite",
        border="sienna"
)


# Normality tests
 
# Shapiro-Wilk normality test
#if p < 0.05, we don't believe that our variables follow a normal distribution

#First on a small scale first to make sure everything runs smoothly
 

s<-shapiro.test(patients_features$Outcome_predictor_score)
print(s$p.value)

for (x in 7:12){
  s<-shapiro.test(patients_features[,x])
  print(s$p.value)
}

sum_s=0.0
for (x in 7:12){
  s<-shapiro.test(patients_features[,x])
  if (s$p.value<0.05) {sum_s =sum_s+1}
  }
print(sum_s)


sum_s = 0.0
for (x in 1:4578){
  s<-shapiro.test(genes_final[,x])
  if (s$p.value<0.05){
    sum_s =sum_s+1
  }
}
print (sum_s)
cat(round(sum_s/4578*100,2),"%")

#We have then only 35.85% of variables which pass the Shapiro-Wilk normality  test


sum_s = 0.0
for (x in 1:4578){
  s<-shapiro.test(genes_features[,x])
  if (s$p.value<0.05){
    sum_s =sum_s+1
  }
}
print(sum_s)
cat(round(sum_s/4578*100,2),"%")


#If we were to check without setting NAs to the median values the percentage increases to 42.22% 



# Jarque-Bera test

j<-jarque.bera.test(patients_features$Outcome_predictor_score)
print(j$p.value)

for (x in 7:12){
  j<-jarque.bera.test(patients_features[,x])
  print(j$p.value)
}

sum_j = 0.0
for (x in 7:12){
  j<-jarque.bera.test(patients_features[,x])
  if (j$p.value<0.05) {sum_j =sum_j+1}
}
print(sum_j)


sum_j = 0.0
for (x in 1:4578){
  j<-jarque.bera.test(genes_final[,x])
  if (j$p.value<0.05){
    sum_j =sum_j +1
  }
}

print(sum_j)
cat(round(sum_j/4578*100,2),"%")

#We now get 35.76% of variables which pass the Jarque-Bera normality test, only a difference of 4 features in total



# Kolmogorov-Smirnov normality test 


My_list <- split(patients, f = list(colnames(patients)))
loop_Shapiro2 <- lapply(My_list, function(x) ks.test(x$Outcome_predictor_score, "pnorm"))
print(loop_Shapiro2)

#we tried using this code for testing the normality of our features however we got very diferent results depending on how we chose the normal distributions to which we compare our features
#We also got lots of warnings (over 50)

sum_k=0.0
for (x in 1:4578){
  k<-ks.test(genes_final[,x],"pnorm",mean=mean(genes_final[,x]),sd=sd(genes_final[,x]) )
  if (k$p.value<0.05){
    sum_k =sum_k+1
  }
}
print(sum_k)
cat(round(sum_k/4578*100,2),"%")

sum_k=0.0
for (x in 1:4578){
  k<-ks.test(genes_final[,x],"pnorm" )
  if (k$p.value<0.05){
    sum_k =sum_k+1
  }
}
print(sum_k)
cat(round(sum_k/4578*100,2),"%")

#It's also easy to see how the results difer from the results obtained above with other tests




#PCA
results <- prcomp(genes_final, scale = TRUE)

#reverse the signs
results$rotation <- -1*results$rotation

#display principal components
results$rotation

summary (results)

#We need to retain 75 PCA's to explain 80% of the variability

#PCA Robust

pc.Grid<- PcaGrid(genes_final, scale=FALSE)
summary(pc.Grid)

#With PP-estimators for PCA using the grid search algorithm we need to retain 113 PCA's to explain 80% of the variability

pc.Grid$eigenvalues

pc.ROBPCA <- PcaHubert(genes_final,kmax=63, scale=FALSE)
summary(pc.ROBPCA)

#With ROBPCA we only need 63 PCA's to explain 80% of the variability

pc.ROBPCA$eigenvalues
pc.ROBPCA$loadings[,1]
pc.ROBPCA$loadings[,2]
View(pc.ROBPCA$scores)























