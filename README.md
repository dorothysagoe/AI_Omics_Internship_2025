# AI_Omics_Internship_2025
Task 1: Scripting in R
# Intro to R for Bioinformatics 🚀

This repository documents my weekly progress in learning **R for bioinformatics programming** as part of the AI Omics Research Internship (2025).  

---

## 📂 Lecture 1 (Class Ib): Getting Started with R  

### 🔑 Topics Covered
1. Setting the working directory properly  
2. Creating and organizing project folders  
3. How R code works: functions, syntax, and execution  
4. Variables and data types in R (numeric, integer, character, factor, logical)  
5. Importing CSV files and working with categorical data  
6. Saving scripts, outputs, and the R workspace  

🎥 [Lecture Recording](https://youtu.be/YybbWfD_VjE?feature=shared)  
📌 [Course GitHub Repo](https://github.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025)  

---

### 📖 Key Learnings
- **Working Directory**: how to set and use working folders so R knows where to look for files and save outputs.  
- **Project Organization**: creating structured subfolders (`data/`, `scripts/`, `results/`) for reproducible research.  
- **R Basics**:  
  - Functions (`mean()`, `plot()`, `hist()`, etc.)  
  - Variables and assignment (`<-`)  
  - Simple data visualizations (scatterplot, histogram, barplot)  
- **Data Types in R**:  
  - Numeric vs Integer  
  - Character / String  
  - Factors for categorical variables  
  - Logical data (`TRUE`/`FALSE`)  
- **Data Handling**:  
  - Importing `.csv` files with `read.csv()`  
  - Checking structure with `str()`  
  - Converting variables into factors or numeric codes (`as.factor()`, `ifelse()`)  
- **Saving Outputs**:  
  - Export cleaned datasets with `write.csv()`  
  - Save workspace and objects (`save()`, `save.image()`)  

---

### 🧬 Assignment / Tasks
1. **Set Working Directory**  
   - Create a new folder `AI_Omics_Internship_2025`.  

2. **Create Project Folder**  
   - In RStudio, make a new project called `Module_I`.  
   - Inside, create subfolders: `raw_data/`, `clean_data/`, `scripts/`, `results/`, `plots/`.  

3. **Data Cleaning Task**  
   - Download `patient_info.csv` from GitHub.  
   - Import the dataset into R.  
   - Inspect structure (`str()`).  
   - Identify variables with incorrect data types.  
   - Convert them to appropriate formats (e.g., factors, numeric).  

4. **Feature Engineering**  
   - Create a new binary variable for smoking status:  
     - `1 = Yes`  
     - `0 = No`  

5. **Save Outputs**  
   - Save cleaned dataset as `clean_data/patient_info_clean.csv`.  
   - Save script as `scripts/class_Ib.R`.  
   - Upload both into this GitHub repository.  

---

### 📌 Example Code Snippets

```r
# Set working directory
setwd("C:/Users/YourName/Documents/AI_Omics_Internship_2025")

# Import CSV
data <- read.csv("raw_data/patient_info.csv")

# Inspect structure
str(data)

# Convert gender to factor
data$gender_fac <- as.factor(data$gender)

# Create binary smoking variable
data$smoking_binary <- ifelse(data$smoking == "Yes", 1, 0)

# Save cleaned dataset
write.csv(data, file = "clean_data/patient_info_clean.csv", row.names = FALSE)



# Intro to R for Bioinformatics 🚀

This repository documents my weekly progress in learning **R for bioinformatics programming** as part of the AI Omics Internship (2025).  

## 📂 Contents
- **Lecture Notes & Scripts**: R scripts from weekly lessons.  
- **Assignments**: My solutions to assignments with explanations.  
- **Projects**: Applications of R in bioinformatics data analysis.  

## 📖 This Week's Focus
### Topic: Differential Expression Analysis & Gene Classification
- Learned how to:
  - Define and use **functions in R**.
  - Apply logical conditions to classify genes as *Upregulated*, *Downregulated*, or *Not Significant*.
  - Handle **missing data** (`NA`) using replacement strategies.
  - Add new columns to data frames (`$status`) for classification results.
  - Save and organize results into a dedicated folder (`Results/`).
  - Summarize results using `table()` to count gene categories.

### Assignment
Classify genes based on `logFC` and `padj` values:
- **Upregulated**: `logFC > 1 & padj < 0.05`
- **Downregulated**: `logFC < -1 & padj < 0.05`
- **Not Significant**: otherwise  

📌 Example function implemented:  

```r
classify_gene <- function(logFC, padj){
  if (logFC > 1 & padj < 0.05){
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05){
    return("Down regulated")
  } else {
    return("Not significant")
  }
}
