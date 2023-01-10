#-----------------------------------------------------------------------------------------------------
# Multivariate Analysis - Project
#
# Samuel Ramos de Pina (78882) & Miguel Bernardo G. Pereira (102211) &
# Joana Ogura () & Thomas Gaehtgens () & Cristiano Mendonça ()
# 
# DATASET USED: Patients & Microarray
# ORIGINAL PAPER: The use of molecular profiling to predict survival after chemotherapy for diffuse large-b-cell lymphoma
#                 New England Journal of Medicine, vol. 346, no. 25, pp. 1937–1947, 2002
# SOURCE: https://llmpp.nih.gov/DLBCL/
#-----------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
#Installation and Loading of Packages
{
#packages instalation
packages <- c("readr", "stringr", "ggplot2", "dplyr", "caret", "tibble", "tidyr", 
              "tseries", "rrcov", "viridis", "ggpmisc", "PerformanceAnalytics", "DataExplorer",
              "creditmodel", "corrplot", "gplots", "factoextra")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
}
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#Part 1 - Data Manipulation

#Data loading
patients <- read_csv("patients.csv")
microarray <- read_csv("microarray.csv")

#Replacing spaces and dashes by underscores
colnames(patients) <- str_replace_all(colnames(patients),c(" " = "_", "-" = "_" ))

#There was an entry in IPI_group with the value "missing", which is the same as NA
patients[patients == "missing"] <- NA

{
#We can link these two data sets by comparing the DLBCL sample (LYM number) in the “patients” 
#with the LYM number in column names of the “microarray” data.
#Creating a modified dataset to enhance interpretability and linkability
microarray_LYM <- microarray
#Renaming the columns for the corresponding LYM number
for (i in 3:276){
  current_name = colnames(microarray_LYM)[i] 
  first_pos_of_LYM = str_locate(current_name, "LYM")[1]
  colnames(microarray_LYM)[i] = str_remove(substring(current_name, first_pos_of_LYM+3, first_pos_of_LYM+5),"^0+")}

#Storing those weird gene names column in a separate dataframe and eliminate that column from the following manipulations
gene_database <- microarray_LYM[1:2]
microarray_LYM <- subset(microarray_LYM, select = -c(2)) 

#merging of the two tables (each row of the table micro_array table corresponds to a sample
#in order to merge it, it is necessary to transpose the table. The name of each variable is the
#concatenation of the first two rows of the table)
microarray_LYM_transposed <- data.frame(t(microarray_LYM[ , 2:(ncol(microarray_LYM))])) #transposing of the table in order to have 292 observations with 7292 variables
colnames(microarray_LYM_transposed) <- microarray_LYM$UNIQID #definition of the name of each variable as the unique ID of a gene
microarray_LYM_transposed$`DLBCL_sample_(LYM_number)` <- rownames(microarray_LYM_transposed)#definition of each patients Lym_number as a variable to be later used on the merging
}

#inner join of the patients and the microarray tables
merged_data_v1 = merge(patients, microarray_LYM_transposed, by = "DLBCL_sample_(LYM_number)") 
patients_features_v1 = merged_data_v1[1:12]
genes_features_v1 = merged_data_v1[13:length(merged_data_v1)]
  
#missing values
sum(is.na(merged_data_v1))
#We have 177352 NAs

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
{
per_NA <- colMeans(is.na(merged_data_v1))
merged_data_v2 <- merged_data_v1[-which(per_NA >= 0.1)]
patients_features_v2 <- merged_data_v2[1:12]
genes_features_v2 <- merged_data_v2[13:length(merged_data_v2)]
}

#With this cut we end up with 4949 features

#Hence, we eliminated 7303-4949=2354 gene features that had 10% of more of NAs
sum(is.na(merged_data_v2))
#Now we have 32375 NAs (we reduced (177352-32375)*100/177352 = 81.75% of NAs)

#Removing highly correlated gene features
#https://stats.stackexchange.com/questions/50537/should-one-remove-highly-correlated-variables-before-doing-pca
{
  #Creating matrix of correlations
  hc = findCorrelation(cor(genes_features_v2, use = "pairwise.complete.obs"), cutoff=0.95) #Eliminate features with more than 0.95 correlation with others
  #Number of Features to Eliminate
  print(length(hc))
  hc = sort(hc)
  reduced_genes_features = genes_features_v2[,-c(hc)]
}

merged_data_v3 <- cbind(patients_features_v2, reduced_genes_features)
patients_features_v3 <- patients_features_v2
genes_features_v3 <- reduced_genes_features

#With this reduction of gene features, we end up with 4578 gene features, and 4590 features in total
#Still a lot of features but we started with 7303!

#Checking how many NAs at this point
sum(is.na(merged_data_v3))
#With this version we have 31669 NAs (17.85% of the initial amount)

#Replacing the NAs with their median value in column
merged_data_v4 <- merged_data_v3 %>% mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE))))
patients_features_v4 <- merged_data_v4[1:12]
genes_features_v4 <- merged_data_v4[13:length(merged_data_v4)]

sum(is.na(merged_data_v4))
sum(is.na(patients_features_v4))

#We replaced every NA except for the IPI subgroup, hence the is.numeric command
sum(is.na(genes_features_v4))

#-----------------------------------------------------------------------------------------------------
#Part 2 - Exploratory Data Analysis
#-----------------------------------------------------------------------------------------------------
#2.1 - PATIENTS FEATURES

#Distribution of NAs across Rows and columns
{
  #Distribution of the % of Missing values by number of features
  hist(colMeans(is.na(merged_data_v1))*100,  col = 1:100, breaks=c(0:100), 
       xlab = '% of NAs', ylab= 'Number of Features', main='') 
  #Distribution of the % of Missing values by number of entries
  hist(rowMeans(is.na(merged_data_v1))*100,  col = 1:100, breaks=c(0:100), 
       xlab = '% of NAs', ylab= 'Number of Entries', main='') 
}

#Minimum, Mean and Maximum of Genes Expressions
{
  max_vals <- mapply(max,merged_data_v1[13:7303], na.rm = TRUE)
  mean_vals <- mapply(mean,merged_data_v1[13:7303], na.rm = TRUE)
  min_vals <- mapply(min,merged_data_v1[13:7303], na.rm = TRUE)
  
  gene_dists <- data.frame(max_vals, mean_vals, min_vals)
  
  boxplot(gene_dists,
          xlab = "Gene Expression Level", ylab = "",
          horizontal = TRUE,
          col=c("red", "orange", "blue"),
          main = "Minimum/Mean/Maximum values of Gene Expression",
          ylim = c(-9, 9))
  axis(1,at=-9:9)
  stripchart(gene_dists, method = "jitter", pch = 19, add = TRUE, col="#8B89890D")
}

#Plot IPI Group distribution by status_at_follow_up(Dead/Alive)
{
  ipi_vs_status <- data.frame(patients_features_v1[,c(4,6)])
  ipi_vs_status <- as.matrix(as.data.frame.matrix(table(ipi_vs_status)))
  barplot(ipi_vs_status,        
          col = c("#1b98e0", "#353436"),
          beside = TRUE)
  legend("topright",                      
         legend = c("Alive", "Dead"),
         fill = c("#1b98e0", "#353436"))
}

#Plot Subgroup of Lymphoma distribution by Status_at_follow_up(Dead/Alive)
{
  subgroup_vs_status <- data.frame(patients_features_v1[,c(4,5)])
  subgroup_vs_status <- as.matrix(as.data.frame.matrix(table(subgroup_vs_status)))
  barplot(subgroup_vs_status,        
          col = c("#1b98e0", "#353436"),
          beside = TRUE)
  legend("topright",                      
         legend = c("Alive", "Dead"),
         fill = c("#1b98e0", "#353436"))
}

#stacked bar chart for IPI_GROUP vs Status
subset(merged_data_v4, IPI_Group %in% c("Medium", "Low", "High")) %>%  
  group_by(Status_at_follow_up, IPI_Group) %>%  
  summarize(Count = n()) %>% 
  ggplot(aes(x=IPI_Group, y=Count, fill=Status_at_follow_up)) + 
  geom_bar(stat='identity', position= "fill")+
  scale_x_discrete(limits = c("Low", "Medium", "High"))+
  scale_fill_manual(values = c("#1b98e0", "#353436"))
#scale_fill_manual(values = c("#b9e38d", "#eb8060"))

{
quantile_predictor = quantile(merged_data_v4$Outcome_predictor_score, probs = c(1/3, 2/3))

predictor_discretized = vector(,dim(merged_data_v4)[1]); #the same as "vec = vector(length = 10);"
predictor_discretized[]="Medium"
predictor_discretized[merged_data_v4$Outcome_predictor_score<=quantile_predictor[1]]="Low"
predictor_discretized[merged_data_v4$Outcome_predictor_score>quantile_predictor[2]]="High"

#stacked bar chart for Predictor(OPS) Discretized vs Status
subset(cbind(merged_data_v4, predictor_discretized), predictor_discretized %in% c("Medium", "Low", "High")) %>%  
  group_by(Status_at_follow_up, predictor_discretized) %>%  
  summarize(Count = n()) %>% 
  ggplot(aes(x=predictor_discretized, y=Count, fill=Status_at_follow_up)) + 
  geom_bar(stat='identity', position= "fill")+
  scale_x_discrete(limits = c("Low", "Medium", "High"))+
  scale_fill_manual(values = c("#1b98e0", "#353436"))
#scale_fill_manual(values = c("#b9e38d", "#eb8060"))
}

{require(gridExtra)
  grid.arrange(subset(merged_data_v4, IPI_Group %in% c("Medium", "Low", "High")) %>%  
                 group_by(Status_at_follow_up, IPI_Group) %>%  
                 summarize(Count = n()) %>% 
                 ggplot(aes(x=IPI_Group, y=Count, fill=Status_at_follow_up)) + 
                 geom_bar(stat='identity', position= "fill")+
                 scale_x_discrete(limits = c("Low", "Medium", "High"))+
                 scale_fill_manual(values = c("#1b98e0", "#353436"))+
                 theme(legend.position = "none"),
               
               subset(cbind(merged_data_v4, predictor_discretized), predictor_discretized %in% c("Medium", "Low", "High")) %>%  
                 group_by(Status_at_follow_up, predictor_discretized) %>%  
                 summarize(Count = n()) %>% 
                 ggplot(aes(x=predictor_discretized, y=Count, fill=Status_at_follow_up)) + 
                 geom_bar(stat='identity', position= "fill")+
                 scale_x_discrete(limits = c("Low", "Medium", "High"))+
                 scale_fill_manual(values = c("#1b98e0", "#353436"))+
                 theme(legend.position = "none"),
               ncol=2
  )}

#Distribution of Predicted Outcome by Group of Lymphoma
ggplot(patients, aes(Outcome_predictor_score, colour = Subgroup, fill=Subgroup)) +
  geom_histogram() + labs(title="-")

#Boxplots for 
{require(gridExtra)
  grid.arrange(boxplot(Outcome_predictor_score~Analysis_Set,
                       data=merged_data_v3,
                       col="floralwhite",
                       border="sienna"),
               boxplot(Outcome_predictor_score~Status_at_follow_up,
                       data=merged_data_v3,
                       col="floralwhite",
                       border="sienna"),
               boxplot(Outcome_predictor_score ~Subgroup,
                       data=merged_data_v3,
                       col="floralwhite",
                       border="sienna"),
               boxplot(Outcome_predictor_score~IPI_Group,
                       data=merged_data_v3,
                       col="floralwhite",
                       border="sienna"), ncol=2)}

#Barplots for Categorical Features in Patients
{require(gridExtra)
grid.arrange(ggplot(patients, aes(x=patients$Analysis_Set, y=..count.., fill=stat(count)))+geom_bar()+labs(x='Analysis set')+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$Status_at_follow_up, y=..count.., fill=stat(count)))+geom_bar()+labs(x='Status at follow-up')+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$Subgroup, y=..count.., fill=stat(count)))+geom_bar()+labs(x='Subgroup')+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$IPI_Group, y=..count.., fill=stat(count)))+geom_bar()+labs(x='IPI group')+theme(legend.position = "none"),#+scale_fill_viridis(), 
             ncol=2)}

#Histograms for Continuous Features in Patients
{require(gridExtra)
grid.arrange(ggplot(patients, aes(x=patients$`Follow_up_(years)`, y=..density.., fill=stat(count))) + labs(x='Follow-up(years)')+
               geom_histogram(colour="black")+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$Germinal_center_B_cell_signature, y=..density.., fill=stat(count))) +labs(x='Germinal center B cell signature')+
               geom_histogram(colour="black")+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$Lymph_node_signature, y=..density.., fill=stat(count))) + labs(x='Lymph node signature')+
               geom_histogram(colour="black")+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$Proliferation_signature, y=..density.., fill=stat(count))) + labs(x='Proliferation signature')+
               geom_histogram(colour="black")+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$BMP6, y=..density.., fill=stat(count))) +labs(x='BMP6')+
               geom_histogram(colour="black")+theme(legend.position = "none"),#+scale_fill_viridis(),
             ggplot(patients, aes(x=patients$MHC_class_II_signature, y=..density.., fill=stat(count))) + labs(x='MHC class II signature') +
               geom_histogram(colour="black")+theme(legend.position = "none"), ncol=2)}

#Correlation(Pearson, I think) Matrix for Continuous Variables in Patients
patients_cont <- patients[,c(7:12)]
colnames(patients_cont) <- c("GCBCs","LNs","Ps","BMP6","MHCIIs", "OPS")
plot_correlation(patients_cont, type = c("continuous")) + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y = element_blank())

#Correlation Chart between Continuous Features
chart.Correlation(patients[,c(7:12)], histogram=TRUE, pch=19)

#Correlation (Cramer's V) Matrix for Nominal Categorical Features
M <- as.matrix(char_cor(patients[,c(2,4,5)]))
corrplot(M, method="color",
         order="hclust", 
         addCoef.col = "black", 
         tl.col="black", 
         tl.srt=45)

#Testing Correlation between Nominal and Ordinal Categorical Features

chisq.test(patients$IPI_Group, patients$Status_at_follow_up)
#Since the p-value is 1.512e-07<<0.05, we reject the null hypothesis that the 
#IPI_group of the patients is not associated with the Status_at_follow_up.
chisq.test(patients$IPI_Group, patients$Analysis_Set)
#Since the p-value is 0.976>>0.05, we DONT reject the null hypothesis that the 
#IPI_group of the patients is not associated with the Analysis_Set.
chisq.test(patients$IPI_Group, patients$Subgroup)
#Since the p-value is 0.4377>>0.05, we DONT reject the null hypothesis that the 
#IPI_group of the patients is not associated with the Subgroup.

#Density of Outcome Predictor Score
ggplot(patients, aes(x=Outcome_predictor_score, fill=stat(density)))+
  geom_histogram(aes(y=..density..), colour="black")+
  #scale_fill_viridis()+
  theme(legend.position = "none")

#ggplots em relacao com a target - categoricas
{require(gridExtra)
grid.arrange(ggplot(patients, aes(x=Outcome_predictor_score, fill=IPI_Group, color=IPI_Group)) +
              geom_histogram(binwidth = 0.1) + labs(title="-")+ theme_bw(),
            ggplot(patients, aes(x=Outcome_predictor_score, fill=Subgroup, color=Subgroup)) +
              geom_histogram(binwidth = 0.1) + labs(title="-")+theme_bw(),
            ggplot(patients, aes(x=Outcome_predictor_score, fill=Status_at_follow_up, color=Status_at_follow_up)) +
              geom_histogram(binwidth = 0.1) + labs(title="-")+theme_bw(), ncol=1)}

#ggplots em relacao com a target - numericas
{require(gridExtra)
grid.arrange(ggplot(patients,aes(x = patients$`Follow_up_(years)`, y = patients$Outcome_predictor_score))+
               geom_point()+
               geom_smooth(method=lm)+
               stat_poly_eq() +
               labs(x="Follow_up_(years)",
                    y="Outcome_predictor_score"),
             ggplot(patients,aes(x = patients$Germinal_center_B_cell_signature, y = patients$Outcome_predictor_score))+
               geom_point()+
               geom_smooth(method=lm)+
               stat_poly_eq() +
               labs(x="Germinal_center_B_cell_signature",
                    y="Outcome_predictor_score"),
             ggplot(patients,aes(x = patients$Lymph_node_signature, y = patients$Outcome_predictor_score))+
               geom_point()+
               geom_smooth(method=lm)+
               stat_poly_eq() +
               labs(x="Lymph_node_signature",
                    y="Outcome_predictor_score"),
             ggplot(patients,aes(x = patients$Proliferation_signature, y = patients$Outcome_predictor_score))+
               geom_point()+
               geom_smooth(method=lm)+
               stat_poly_eq() +
               labs(x="Proliferation_signature",
                    y="Outcome_predictor_score"),
             ggplot(patients,aes(x=BMP6,y=Outcome_predictor_score))+
               geom_point()+
               geom_smooth(method=lm)+
               stat_poly_eq() +
               labs(x="BMP6",
                    y="Outcome_predictor_score"),
             ggplot(patients,aes(x = patients$MHC_class_II_signature, y = patients$Outcome_predictor_score))+
               geom_point()+
               geom_smooth(method=lm)+
               stat_poly_eq() +
               labs(x="MHC_class_II_signature",
                    y="Outcome_predictor_score"),ncol=2)}

#-----------------------------------------------------------------------------------------------------
#2.2 - GENES FEATURES

#VISUALIZING GENES V1
genes_matrix <- data.matrix(genes_features_v1, rownames.force = NA)
genes_matrix[genes_matrix < -2 ] <- -2
genes_matrix[genes_matrix > 2 ] <- 2

Samples = 1:nrow(genes_matrix)
Genes = 1:ncol(genes_matrix)

image(Samples,
      Genes,
      genes_matrix,
      col = colorpanel(25, "green", "black", "red"))

genes_values <- data.frame(log_expressivity = c(genes_matrix))
genes_values$expressivity <- 2**genes_values$log_expressivity

ggplot(genes_values, aes(x=genes_values$expressivity, y=..count.., fill=stat(count)))+ 
  labs(x='expressivity')+
  geom_histogram(colour="black", bins=100)+
  xlim(0, 7)+theme(legend.position = "none")

#VISUALIZING GENES V3
genes_matrix <- data.matrix(genes_features_v3, rownames.force = NA)
genes_matrix[genes_matrix < -2 ] <- -2
genes_matrix[genes_matrix > 2 ] <- 2

Samples = 1:nrow(genes_matrix)
Genes = 1:ncol(genes_matrix)

image(Samples,
      Genes,
      genes_matrix,
      col = colorpanel(25, "green", "black", "red"))

genes_values <- data.frame(log_expressivity = c(genes_matrix))
genes_values$expressivity <- 2**genes_values$log_expressivity

ggplot(genes_values, aes(x=genes_values$expressivity, y=..count.., fill=stat(count)))+ 
  labs(x='expressivity')+
  geom_histogram(colour="black", bins=100)+
  xlim(0, 7)+theme(legend.position = "none")

#VISUALIZING GENES V4
genes_matrix <- data.matrix(genes_features_v4, rownames.force = NA)
genes_matrix[genes_matrix < -2 ] <- -2
genes_matrix[genes_matrix > 2 ] <- 2

Samples = 1:nrow(genes_matrix)
Genes = 1:ncol(genes_matrix)

image(Samples,
      Genes,
      genes_matrix,
      col = colorpanel(25, "green", "black", "red"))

genes_values <- data.frame(log_expressivity = c(genes_matrix))
genes_values$expressivity <- 2**genes_values$log_expressivity

ggplot(genes_values, aes(x=genes_values$expressivity, y=..count.., fill=stat(count)))+ 
  labs(x='expressivity')+
  geom_histogram(colour="black", bins=100)+
  xlim(0, 7)+theme(legend.position = "none")


#-----------------------------------------------------------------------------------------------------
#2.3 -  NORMALITY TESTS

# Shapiro-Wilk normality test
#if p < 0.05, we don't believe that our variables follow a normal distribution

for (x in 7:12){
  s<-shapiro.test(patients_features_v4[,x])
    if (s$p.value>=0.05){
    print(colnames(patients_features_v4)[x])}
}
# Lymph node signature, Proliferation signature and Outcome predictor score pass the test
{
  sum_s = 0.0
  for (x in 1:4578){
    s<-shapiro.test(genes_features_v4[,x])
    if (s$p.value>=0.05){
      sum_s =sum_s+1}}
  print (sum_s)
  cat(round(sum_s/4578*100,2),"%")
}
# We have then only 35.85% of variables which pass the Shapiro-Wilk normality  test

{
  sum_s = 0.0
  for (x in 1:4578){
    s<-shapiro.test(genes_features_v3[,x])
    if (s$p.value>=0.05){
      sum_s =sum_s+1}}
  print(sum_s)
  cat(round(sum_s/4578*100,2),"%")
}
# If we were to check without setting NAs to the median values the percentage increases to 41.94% 


# Jarque-Bera test
#if p < 0.05, we don't believe that our variables follow a normal distribution

for (x in 7:12){
  j<-jarque.bera.test(patients_features_v4[,x])
  if (j$p.value>=0.05){
    print(colnames(patients_features_v4)[x])}
}
# Again Lymph node signature, Proliferation signature and Outcome predictor score pass the test
{
  sum_j = 0.0
  for (x in 1:4578){
    j<-jarque.bera.test(genes_features_v4[,x])
    if (j$p.value>=0.05){
      sum_j =sum_j +1}}
  print(sum_j)
  cat(round(sum_j/4578*100,2),"%")
}
#We now get 35.76% of variables which pass the Jarque-Bera normality test, only a difference of 4 features in total




#-----------------------------------------------------------------------------------------------------
#Part 3 - ROBUST PCA
#-----------------------------------------------------------------------------------------------------

pc.Grid<- PcaGrid(genes_features_v4, scale=FALSE)
summary(pc.Grid)
#With PP-estimators for PCA using the grid search algorithm we need to retain 113 PCA's to explain 80% of the variability


pc.ROBPCA <- PcaHubert(genes_features_v4,kmax=63, scale=FALSE)
summary(pc.ROBPCA)
#With ROBPCA we only need 63 PCA's to explain 80% of the variability


#CLUSTERING

df <- genes_features_v4
res.dist <- dist(df, method = "euclidean")
as.matrix(res.dist)[1:6, 1:6]
res.hc <- hclust(d = res.dist, method = "ward.D2")
#possible methods: “ward.D”, “ward.D2”, “single”,
#“complete”, “average”, “mcquitty”, “median” or “centroid”
fviz_dend(res.hc, cex = 0.5)
# Compute cophentic distance
res.coph <- cophenetic(res.hc)
# Correlation between cophenetic distance and # the original distance
cor(res.dist, res.coph)
res.hc2 <- hclust(res.dist, method = "average") 
cor(res.dist, cophenetic(res.hc2))
# Cut tree into 3 groups
grp <- cutree(res.hc, k = 3) 
head(grp, n = 3)
# Number of members in each cluster
table(grp)
# Get the names for the members of cluster 1
rownames(df)[grp == 1]

# Cut in 3 groups and color by groups 
fviz_dend(res.hc, k = 3, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"), color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

fviz_cluster(list(data = df, cluster = grp),
             palette = c("#2E9FDF", "#E7B800", "#FC4E07"), ellipse.type = "convex", # Concentration ellipse
             repel = TRUE, # Avoid label overplotting (slow) 
             show.clust.cent = FALSE, ggtheme = theme_minimal())

#VISUALIZING GENES V4 AFTER CLUSTERING
genes_w_cluster <- genes_features_v4
genes_w_cluster$cluster <- as.vector(grp)
genes_w_cluster <- genes_w_cluster %>% arrange(cluster)

genes_w_cluster_matrix <- data.matrix(genes_w_cluster, rownames.force = NA)
genes_w_cluster_matrix[genes_w_cluster_matrix < -2 ] <- -2
genes_w_cluster_matrix[genes_w_cluster_matrix > 2 ] <- 2

Samples = 1:nrow(genes_w_cluster_matrix)
Genes = 1:ncol(genes_w_cluster_matrix)

image(Samples,
      Genes,
      genes_w_cluster_matrix,
      col = colorpanel(25, "green", "black", "red"))

######################################
#             K-MEDOIDS              #
######################################
{
  library("fpc")
  library("dbscan")
  library("factoextra")
  library("cluster")
  library("ggplot2")
  library(plyr)
}

data<-read.csv('patients_clean_standardized.csv')

status<-data$Status.at.follow.up
# Remove  Columns in List
data <- data[,!names(data) %in% c("Status.at.follow.up")]
data$IPI.Group = as.numeric(data$IPI.Group) 

fviz_nbclust(data, pam, method = "wss")

#perform k-medoids clustering with k = 2 clusters
kmed <- pam(data, k = 2)

#view results
kmed

#plot results of final k-medoids model
fviz_cluster(kmed, data = data)

{
kmed_mad<-mapvalues(kmed$clustering, from = c(1, 2), to = c("Dead", "Alive"))
table(kmed_mad, status)
}

{
  kmed_mad<-mapvalues(kmed$clustering, from = c(1, 2), to = c("Alive", "Dead"))
  table(kmed_mad, status)
}
