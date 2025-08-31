# AI_Omics_Internship_2025
Task 1: Scripting in R


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
