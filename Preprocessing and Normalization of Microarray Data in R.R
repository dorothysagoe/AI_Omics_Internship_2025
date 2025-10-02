# Assignment III

################################################################

### Pre-processing and Normalization of Microarray Data in R ###

################################################################



# Objectives of this assignment

# Perform quality control before and after normalization and check whether any arrays are flagged as outliers. (note down how many you found before and after normalization)
# Normalize the data and then apply filtering to remove low-intensity probes and note how many transcripts remain.
# Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)

# Setting working directory
setwd('C:\\Users\\Dorothy\\Documents\\AI_Omics_Internship_2025\\Assignment III')

#Creating directory for Raw_Data file
Raw_Data = dir.create("C:\\Users\\Dorothy\\Documents\\AI_Omics_Internship_2025\\Assignment III\\Raw_Data")

#importing libraries
library(arrayQualityMetrics)
library(ArrayExpress)
library(affy)
library(GEOquery)

# Downloading files from array express
gse_data <- getGEO("GSE79973", GSEMatrix = TRUE)

# Extract expression data
expression_data <- exprs(gse_data$GSE79973_series_matrix.txt.gz)

# Extract feature data 
feature_data <- fData(gse_data$GSE79973_series_matrix.txt.gz)

#Extract phenotype data 

phenotype_data <- pData(gse_data$GSE79973_series_matrix.txt.gz)


# Reading raw data CEL files
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

# Data quality metrics before normalization
quality_metrics <- arrayQualityMetrics(expressionset=raw_data,outdir = "Results/GC_Metrics_Raw")

# Normalization of raw data
normalized_data <- rma(raw_data)

# Data quality metrics after normalization
arrayQualityMetrics(expressionset=normalized_data,
                    outdir = "Results/GC_Metrics_Normalized",
                    force=TRUE)

# Extraction of expression of normalized values
processed_data <- as.data.frame(exprs(normalized_data ))

# Inspecting the dimension of extracted file
dim(processed_data)

############################################################################
###Filtration of Low transcript reads###
############################################################################

# Computing median to determine median process
row_median <- rowMedians(as.matrix(processed_data))

read.table(row_median)

# Plotting row_median histogram plot to visualize distribution

Median_Intensity_Distribution<-hist(row_median,
     breaks=100,
     freq=FALSE,
     main="Median Intensity Distribution")

dir.create('Plots')


# Saving plotted files as pdf
pdf(file = "Plots/Median Intensity Distribution.pdf", width = 8, height = 6)

plot(hist(row_median,
          breaks=100,
          freq=FALSE,
          main="Median Intensity Distribution"))

# Closing the device to finalize saving
dev.off()
  

# Setting threshold for filtration for probes to be 3.5
threshold <- 3.5

# Indicating Threshold on plot
abline(v=threshold, col="red", lw=2)

indx <- row_median > threshold

filtered_data <- processed_data[indx, ]

dim(filtered_data)

#Renaming columns

# Verifying if length of rows and column names match
length(colnames(filtered_data))
length(rownames(phenotype_data))

#Overwriting column names with filtered_data
colnames(filtered_data) <- rownames(phenotype_data)

# Viewing newly formatted filtered_data df
View(filtered_data)

# Passing filtered_data df into "Processed_data"
Processed_data <- filtered_data

dim(Processed_data)

############################################################################
###Phenotype processing###
############################################################################

#Inspecting class of phenotype source name
class(phenotype_data$sources_name_ch1)

#Converting it to a factor 
groups <- factor(phenotype_data$sources_name_ch1,
                 levels = c("gastric adenocarcinoma", "gastric mucosa"),
                 labels = c("normal", "cancer"))

class(groups)



##################### End of Preprocessing and Normalization ################
