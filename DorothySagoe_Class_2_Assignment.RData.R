# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Assignment II: Getting Started with R
# -------------------------------------------------------------------
# ===================================================================


# Printing current working directory
getwd()

# Changing working directory to desired

setwd("C://Users//Dorothy//Documents//AI_Omics_Internship_2025//Assignment_II")


#Creating folder directories for Raw_Data, Results and Plots

dir.create("Raw_Data")
dir.create("Results")

#Moving files from downloads to Raw_Data folder

from1 <- "C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/DEGs_Data_1.csv"
to1   <- "C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/Assignment_II/Raw_Data/DEGs_Data_1.csv"

file.rename(from1, to1)

from2 <- "C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/DEGs_Data_2.csv"
to2   <- "C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/Assignment_II/Raw_Data/DEGs_Data_2.csv"

file.rename(from2, to2)

setwd("C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/Assignment_II/Raw_Data")

#viewing my dataset
DEGs_Data_1 <- read.csv("DEGs_Data_1.csv", header=TRUE)

head(DEGs_Data_1)

DEGs_Data_2 <- read.csv("DEGs_Data_2.csv", header=TRUE)

head(DEGs_Data_2)


# Defining classify_gene function
classify_gene <- function(logFC, padj){
  
  if (logFC > 1 & padj < 0.05){
    return("Upregulated")
  }
  
  else if (logFC < -1 & padj <0.05){
    return("Down regulated")
  }  
  else{
    return("Not significant")
  
  }
}


# Calling classify_gene function for LogFC = 2 and padj = 0.04
classify_gene(2, 0.04)


# Viewing NA 
sum(is.na(DEGs_Data_1$padj))

# Replacing missing padj values with 1
DEGs_Data_1$padj[is.na(DEGs_Data_1$padj)] <- 1
DEGs_Data_2$padj[is.na(DEGs_Data_2$padj)] <- 1


# Adding a new column called status
DEGs_Data_1$status <- mapply(classify_gene, DEGs_Data_1$logFC, DEGs_Data_1$padj)
DEGs_Data_2$status <- mapply(classify_gene, DEGs_Data_2$logFC, DEGs_Data_2$padj)

# Viewing changes in first 6 rows
head(DEGs_Data_1)

#Creating results directory
Results <- "C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/Assignment_II/Results"
if (!dir.exists(Results)) dir.create(Results)



# Saving processed files into Results folder
# Save DEGs_Data_1
write.csv(DEGs_Data_1,
          file = file.path(Results, "DEGs_Data_1_results.csv"),
          row.names = FALSE)

# Save DEGs_Data_2
write.csv(DEGs_Data_2,
          file = file.path(Results, "DEGs_Data_2_results.csv"),
          row.names = FALSE)


# Summary counts for DEGs_Data_1
cat("DEGs_Data_1 summary:\n")
print(table(DEGs_Data_1$status))

# Summary counts for DEGs_Data_2
cat("\nDEGs_Data_2 summary:\n")
print(table(DEGs_Data_2$status))



