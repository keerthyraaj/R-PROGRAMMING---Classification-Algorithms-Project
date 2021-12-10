##################################################
###                                             ##
### RStudio - Assignment 5                      ## 
##################################################
#                                               ##
##################################################
# Written by Keerthy Raaj Shanmugam
# ID: 8779954
#
##################################################
### Basic Set Up                                ##
##################################################

# Clear all plots
if(!is.null(dev.list())) dev.off()

# Clear entire console
cat("\014") 

# Clean and clear theworkspace
rm(list=ls())

#Set work directory to an appropriate location
setwd("C:/Users/keert/Documents/Data")

options(scipen=9)

##################################################
### Remove Packages Installed                   ##
##################################################

##################################################
### Install Libraries                           ##
##################################################

#If the library is not already downloaded, download it

#For ROC Curves

if(!require(pROC)){install.packages("pROC")}
library(pROC)

#For N-B Analysis

if(!require(klaR)){install.packages("klaR")}
library("klaR")

# For LDA

if(!require(MASS)){install.packages("MASS")}
library("MASS")

#Read the Tumor_21F data which is in CSV format
Tum_dat <- read.csv("Tumor_21F.csv",header=TRUE,sep=",")
head(Tum_dat)
str(Tum_dat)

#Rename Variables 
names(Tum_dat) <- c("Age_KS", "Sex_KS", "Bone_KS", "Marrow_KS", "Lung_KS",
                    "Pleura_KS", "Liver_KS", "Brain_KS", "Skin_KS",
                    "Neck_KS","Supra_KS", "Axil_KS", "Media_KS", "Out_KS")
str(Tum_dat)

#####################  Part-A  ######################

#1.1 Preliminary Data Preparation

#Apply the Missing Value Filter to remove appropriate columns of data.
summary(Tum_dat)

#Apply the Low Variance Filter to remove appropriate columns of data.
stat.desc(Tum_dat) #Consider coef of var

#Apply the High Correlation Filter to remove appropriate columns of data.
hcf <- Tum_dat[,c("Age_KS", "Sex_KS", "Bone_KS", "Marrow_KS", "Lung_KS",
                  "Pleura_KS", "Liver_KS", "Brain_KS", "Skin_KS",
                  "Neck_KS","Supra_KS", "Axil_KS", "Media_KS", "Out_KS")]

cor(hcf,method="spearman")

#Comment on any outliers you see and deal with them appropriately.
summary(Tum_dat)

#2 Exploratory Analysis
#2.1 Correlations: Create numeric correlations (as demonstrated) and comment on what you see. Are there co-linear variables? 
Tum_dat_Corr <- Filter(is.numeric, Tum_dat)  # Only take numeric
str(Tum_dat_Corr)

Tumdat_Corr <- cor(Tum_dat_Corr, method="spearman")
round(Tumdat_Corr, 2) 

#2.2 Identify the two most significant predictors of tumors and provide statistical evidence (in addition to the correlation coefficients) that suggest they are associated with tumors (Think of the contingency tables we did in class). 

#Brain
Tbl_Brain <- table(Tum_dat$Out_KS, Tum_dat$Brain_KS, dnn=list("Out","Brain"))
Tbl_Brain
prop.table(Tbl_Brain, 2) # col percentages

#Check the Chi Squared Test - NOTE Removal of Yate's Continuity Correction

chisq_Brain <- chisq.test(Tum_dat$Out_KS, Tum_dat$Brain_KS, correct=FALSE)      
chisq_Brain

chisq_Brain$observed   # What we observed
chisq_Brain$expected   # If there were no relationship


#Marrow
Tbl_Marrow <- table(Tum_dat$Out_KS, Tum_dat$Marrow_KS, dnn=list("Out","Marrow"))
Tbl_Marrow
prop.table(Tbl_Marrow, 2) # col percentages

#Check the Chi Squared Test - NOTE Removal of Yate's Continuity Correction

chisq_Marrow <- chisq.test(Tum_dat$Out_KS, Tum_dat$Marrow_KS, correct=FALSE)      
chisq_Marrow

chisq_Marrow$observed   # What we observed
chisq_Marrow$expected   # If there were no relationship

#3.1.Forward selection model.
#stepwise model
Tum_glm = glm(Out_KS ~ Age_KS + Sex_KS + Bone_KS + Marrow_KS + Lung_KS
              + Pleura_KS + Liver_KS + Brain_KS + Skin_KS +Neck_KS + 
              Supra_KS + Axil_KS + Media_KS,
              family="binomial",data=Tum_dat, na.action=na.omit)

#3.1.Forward selection model.
Tum_fsm <- step(Tum_glm,direction="forward", details=TRUE)
summary(Tum_fsm)


#3.2. Two additional models using variables that you select based on the above output (recall lecture slides on variable selection). We will refer to these models as "User Model 1" and "User Model 2".
#User Model 1
Tum_User_Model_1 = glm(Out_KS ~ Age_KS + Sex_KS + Marrow_KS + Lung_KS + Liver_KS + 
                   Brain_KS + Skin_KS +Neck_KS + Supra_KS + Axil_KS + Media_KS,
                   family="binomial",data=Tum_dat, na.action=na.omit)

summary(Tum_User_Model_1)

#User Model 2
Tum_User_Model_2 = glm(Out_KS ~ Age_KS + Sex_KS + Marrow_KS + Lung_KS  
                       + Brain_KS + Neck_KS + Supra_KS + Media_KS,
                    family="binomial",data=Tum_dat, na.action=na.omit)

summary(Tum_User_Model_2)


#4.1 For User Model 1 and User Model 2, create and evaluate the confusion matrices. Set the default predictive level to 50% for "success"

# Check the User Model 1
resp_UM1 <- predict(Tum_User_Model_1, type="response")   # creates probabilities
head(resp_UM1,20)
Class_UM1 <- ifelse(resp_UM1 > 0.5,1,0)           # Classifies probablities (i.e. >50% then likely to have Tumor)
head(Class_UM1)
True_log_1 <- Tum_dat$Out_KS                        #Creates a vector of the true outcomes
T1 <- table(True_log_1, Class_UM1, dnn=list("Act Tumor","Predicted") )  # Creates a Contingency Table
T1

# Check the User Model 2
resp_UM2 <- predict(Tum_User_Model_2, type="response")   # creates probabilities
head(resp_UM2,20)
Class_UM2 <- ifelse(resp_UM2 > 0.5,1,0)           # Classifies probablities (i.e. >50% then likely to have Tumor)
head(Class_UM2)
True_log_2 <- Tum_dat$Out_KS                        #Creates a vector of the true outcomes
T2 <- table(True_log_2, Class_UM2, dnn=list("Act Tumor","Predicted") )  # Creates a Contingency Table
T2

#ROC Curve (and Area Under the Curve)

#User Model 1
plot(roc(Tum_dat$Out_KS ,resp_UM1, direction="<"),
     col="red", lwd=2, main='ROC Curve for Logistic, Tumor Prediction')

auc(Tum_dat$Out_KS, resp_UM1)

#User Model 2
plot(roc(Tum_dat$Out_KS ,resp_UM2, direction="<"),
     col="red", lwd=2, main='ROC Curve for Logistic, Tumor Prediction')

auc(Tum_dat$Out_KS, resp_UM2)



#####################  Part-B  ######################

#stepwise

#1.1 As above, use the forward option in the glm function to fit the model NOTE - These results should match the output from 3.1.
start_time <- Sys.time()

Tum_glms = glm(Out_KS ~ Age_KS + Sex_KS + Bone_KS + Marrow_KS + Lung_KS
              + Pleura_KS + Liver_KS + Brain_KS + Skin_KS +Neck_KS + 
                Supra_KS + Axil_KS + Media_KS,
              family="binomial",data=Tum_dat, na.action=na.omit)
stp_Tum_glm <- step(Tum_glms, direction="forward", details=TRUE)

end_time <- Sys.time()

SW_Time <- end_time - start_time

summary(stp_Tum_glm)

#1.2 Summarize the results in a Confusion Matrix .
resp_SM <- predict(stp_Tum_glm, type="response")   # creates probabilities
head(resp_SM,20)
Class_SM <- ifelse(resp_SM > 0.5,1,0)           # Classifies probablities (i.e. >50% then likely to have Tumor)
head(Class_SM)
True_log_SM <- Tum_dat$Out_KS                        #Creates a vector of the true outcomes
T_SM <- table(True_log_SM, Class_SM, dnn=list("Act Tumor","Predicted") )  # Creates a Contingency Table
T_SM

#1.3 As demonstrated in class, calculate the time (in seconds) it took to fit the model and include this in your summary.
SW_Time


#NAIVE BAYES

#2.1 As demonstrated in class, transform the variables as necessary for N-B classification.
str(Tum_dat)

Tum_dat$Out_KS <- as.factor(Tum_dat$Out_KS)

str(Tum_dat)


#2.2.Use all the variables in the dataset to fit a Naïve-Bayesian classification model. 
start_time <- Sys.time()

Tumor_Naive <- NaiveBayes(Out_KS ~ Age_KS + Sex_KS + Bone_KS + Marrow_KS + Lung_KS
                          + Pleura_KS + Liver_KS + Brain_KS + Skin_KS +Neck_KS + 
                            Supra_KS + Axil_KS + Media_KS,
                          data=Tum_dat, na.action=na.omit)
end_time <- Sys.time()

NB_Time <- end_time - start_time

#Classifies
pred_Tum_bay <- predict(Tumor_Naive,Tum_dat)
head(pred_Tum_bay$class,5)
head(pred_Tum_bay$posterior,5)


#2.3.Summarize the results in a Confusion Matrix.
CF_NB <- table(Actual=Tum_dat$Out_KS, Predicted=pred_Tum_bay$class)
CF_NB

#2.4.As demonstrated in class, calculate the time (in seconds) it took to fit the model and include this in your summary. 
NB_Time


#3.1. As demonstrated in class, transform the variables as necessary for LDA classification.
Tum_dat$Out_KS <- as.factor(Tum_dat$Out_KS)

str(Tum_dat)

#3.2. Use all the variables in the dataset to fit an LDA classification model.
start_time <- Sys.time()

Tumor_Discrim <- lda(Out_KS ~ Age_KS + Sex_KS + Bone_KS + Marrow_KS + Lung_KS
                     + Pleura_KS + Liver_KS + Brain_KS + Skin_KS +Neck_KS + 
                       Supra_KS + Axil_KS + Media_KS,
                     data=Tum_dat, na.action=na.omit)

end_time <- Sys.time()

LDA_Time <- end_time - start_time


#Classifies
pred_Tum_dis <- predict(Tumor_Discrim, Tum_dat)
head(pred_Tum_dis$class,5)
head(pred_Tum_dis$posterior,5)


#3.3. Summarize the results in a Confusion Matrix.
#Confusion Matrix
CF_LDA <- table(Actual=Tum_dat$Out_KS, Predicted=pred_Tum_dis$class)
CF_LDA

#3.4. As demonstrated in class, calculate the time (in seconds) it took to fit the model and include this in your summary.
LDA_Time

#4.1. Which classifier is most accurate? (provide evidence)
#4.2. Which classifier is most suitable when processing speed is most important?
#4.3. Which classifier minimizes Type 1 errors?
#4.4. Which classifier minimizes Type 2 errors?
#4.5. Which classifier is best overall?

#Confusion Matrix
T_SM
CF_NB
CF_LDA

#time (in seconds)
SW_Time
NB_Time
LDA_Time

#4.6. How do these classifiers compare to the best model you built in Part 1?
T2
