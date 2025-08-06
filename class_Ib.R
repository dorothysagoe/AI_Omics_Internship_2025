#Set a working directory
setwd("C:\\Users\\Dorothy\\Documents\\AI_Omics_Internship_2025")

#Organize files
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

#loading patient_info data set
patient_info_raw_data <- read.csv("raw_data/patient_info.csv")

#Viewing raw data csv file
View(patient_info_raw_data)

#Checking data structure of patient_info data set
str(patient_info_raw_data)

#Converting smoker column to a factor
patient_info_raw_data$smoker_fac <- as.factor(patient_info_raw_data$smoker)

#checking the structure of gender_fac
str(patient_info_raw_data$smoker_fac)

#Converting gender column to numeric
patient_info_raw_data$smoker_num <- ifelse(patient_info_raw_data$smoker_fac =="Yes",1, 0)

#Viewing data
View(patient_info_raw_data)

#Dropping smoker and smoker columns
patient_info_raw_data$smoker <- NULL

#Renaming dataset after data cleaning
patient_info_clean <- patient_info_raw_data

#Save cleaned data set into csv file
write.csv(patient_info_clean, "clean_data/patient_info_clean.csv", row.names = FALSE)

#View cleaned dataset
View(patient_info_clean)



