# Assignment III

################################################################

### Differential Gene Expression Analyses of Microarray Data in R ###

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
library(AnnotationDbi)
library(hgu133plus2.db)
library(limma)
library(tibble)
library(ggplot2)
library(pheatmap)
library(dplyr)

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

#read.table(row_median)

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

#Saving processed data into a file

save(Processed_data, file="GSE79973.RData")

##################### End of Pre-processing and Normalization ################


#Check annotation slot of the data set
annotation(raw_data)

#View objects available in the package
ls("package:hgu133plus2.db")

columns(hgu133plus2.db)
keytypes(hgu133plus2.db)


#Extract probe IDs from processed microarray data
probe_ids <- rownames(processed_data)

#Mapping probe IDs gene symbols using the platform annotation database
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multivals = "first"
)

#Convert mapping to a data frame and rename columns
gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBE_ID") %>%
  dplyr::rename(SYMBOL = 2)

#Address many-to-one by average or summarize probe signals to maintain one row per gene

#Summarize number of probes per gene symbol
duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))
  
# Identify gene associated with many probes
  duplicate_genes <- duplicate_summary %>%
    filter(probes_per_gene >1)

sum(duplicate_genes$probes_per_gene)


#Verify
all(gene_map_df$PROBE_ID == row.names(processed_data))

#Merge annotation (SYMBOL)
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)


#Remove probes without valid gene symbol anotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

#select only numerical data
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

# -------------------------------------------------------------
# Collapse multiple probes per gene using average expression
# -------------------------------------------------------------
# limma::avereps() computes the average for probes representing the same gene
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

dim(averaged_data)

# Convert averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix

#Differential Gene Expression Analysis

# Define sample groups based on phenotype data
# Adjust group labels according to dataset annotation
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)

# -------------------------------------------------------------
# Create design matrix for linear modeling
# -------------------------------------------------------------
# Using no intercept (~0 + groups) allows each group to have its own coefficient
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit linear model to expression data
fit_1 <- lmFit(data, design)

# Define contrast to compare cancer vs normal samples
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# -------------------------------------------------------------
# Extract list of differentially expressed genes (DEGs)
# -------------------------------------------------------------
deg_results <- topTable(fit_2,
                        coef = "cancer_vs_normal",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction

# -------------------------------------------------------------
# Classify DEGs into Upregulated, Downregulated, or Not Significant
# -------------------------------------------------------------
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))

# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

# Combine both sets of DEGs
deg_updown <- rbind(upregulated, downregulated)

write.csv(deg_results, file = "Results/DEGs_Results.csv")
write.csv(upregulated, file = "Results/Upregulated_DEGs.csv")
write.csv(downregulated, file = "Results/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Results/Updown_DEGs.csv")



############################################################################
###Data Visualization###
############################################################################

# -------------------------------------------------------------
# Volcano Plot: visualizes DEGs by logFC and adjusted p-values
# -------------------------------------------------------------
# Note: x-axis = log2 fold change, y-axis = -log10 adjusted p-value

# Save volcano plot as PNG
dir.create("Results_Plots")
png("Results_Plots/volcano_plot.png", width = 2000, height = 1500, res = 300)

volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")

print(volcano_plot)

dev.off()

# -------------------------------------------------------------
# Heatmap of Top Differentially Expressed Genes
# -------------------------------------------------------------

# Select top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 10)

# Subset averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names

# Save heatmap as PNG
png("Results_Plots/heatmap_top50_DEGs.png", width = 2000, height = 1500, res = 300)

# Generate heatmap without additional scaling
pheatmap(
  heatmap_data,
  scale = "none", # for already normalized data
  cluster_rows = FALSE,              # Disable row clustering
  cluster_cols = TRUE,               # Cluster samples
  show_rownames = TRUE,              # Display gene names
  show_colnames = TRUE,              # Display sample labels
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 10 Differentially Expressed Genes"
)

dev.off()



