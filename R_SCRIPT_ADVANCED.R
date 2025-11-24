
#Name:Kayyari Deeksha
#R_ADVANCED_ASSIGNMENT





# ============================================================
# LOAD REQUIRED PACKAGES
# ============================================================
library(readxl)   # for Excel files
library(readr)    # for CSV/TSV imports
library(dplyr)    # for data wrangling
library(tidyr)    # for splitting columns
library(stringr)  # for string cleanup


# Your folder path
path <- "C:/Users/deeksha/OneDrive - Indiana University/Documents"

# ========================
# ROBUSTNESS: FILE CHECKS
# ========================
required_files <- c(
  "Gene_Expression_Data.xlsx",
  "Gene_Information.csv",
  "Sample_Information.tsv"
)

for (f in required_files) {
  f_path <- file.path(path, f)
  if (!file.exists(f_path)) {
    stop(paste("ERROR: File not found →", f_path))
  }
}


# ============================================================
# 1. LOAD GENE EXPRESSION DATA (EXCEL)
# ============================================================
gene_expression <- read_excel(file.path(path, "Gene_Expression_Data.xlsx"))

# ------------------ QC CHECKS -------------------------------
cat("QC: Gene Expression Data\n")
cat("  Rows:", nrow(gene_expression), "\n")
cat("  Columns:", ncol(gene_expression), "\n")
str(gene_expression)
cat("\n------------------------------------------------\n\n")

# ============================================================
# 2. LOAD GENE INFORMATION DATA (CSV)
# ============================================================
gene_info <- read_csv(file.path(path, "Gene_Information.csv"))

# ------------------ QC CHECKS -------------------------------
cat("QC: Gene Information Data\n")
cat("  Rows:", nrow(gene_info), "\n")
cat("  Columns:", ncol(gene_info), "\n")
str(gene_info)
cat("\n------------------------------------------------\n\n")

# ============================================================
# 3. LOAD SAMPLE INFORMATION DATA (TSV)
# ============================================================
sample_info <- read_tsv(file.path(path, "Sample_Information.tsv"))

# ------------------ QC CHECKS -------------------------------
cat("QC: Sample Information Data (Before Cleaning)\n")
cat("  Rows:", nrow(sample_info), "\n")
cat("  Columns:", ncol(sample_info), "\n")
str(sample_info)
cat("\n------------------------------------------------\n\n")

# ============================
# ROBUSTNESS: SAFE FILE READER
# ============================
safe_read <- function(expr) {
  tryCatch(expr, error = function(e) {
    cat("ERROR:", e$message, "\n")
    stop("Stopping due to read failure.")
  })
}


# ============================================================
# 4. CLEAN SAMPLE INFORMATION TSV
# ============================================================
sample_info_clean <- sample_info %>%
  # Split on the tab inside the patient column
  separate(patient, into = c("tissue", "patient_raw"), sep = "\t", fill = "right") %>%
  
  # Clean up text
  mutate(
    tissue = str_trim(tissue),
    patient = str_replace(patient_raw, "patient:\\s*", ""),
    patient = as.numeric(patient)
  ) %>%
  
  select(group, tissue, patient)

# ------------------ QC CHECKS -------------------------------
cat("QC: Sample Information Data (AFTER Cleaning)\n")
cat("  Rows:", nrow(sample_info_clean), "\n")
cat("  Columns:", ncol(sample_info_clean), "\n")
str(sample_info_clean)
cat("\n------------------------------------------------\n\n")

# ------------------------------------------------------------
# 1.b Rename samples in gene_expression based on phenotype
# ------------------------------------------------------------

# Create a lookup table: GSM ID → new name
sample_info_clean <- sample_info_clean %>%
  mutate(new_name = paste(tissue, patient, sep = "_"))

# Show mapping QC
print(sample_info_clean)

# Extract original sample names (exclude Probe_ID)
original_samples <- colnames(gene_expression)[-1]

# Create a named vector for renaming
rename_vector <- sample_info_clean$new_name
names(rename_vector) <- sample_info_clean$group

# Apply renaming to gene_expression
gene_expression_renamed <- gene_expression

colnames(gene_expression_renamed)[-1] <- rename_vector[colnames(gene_expression)[-1]]

# QC: show new column names
cat("QC: New Column Names in Gene Expression Data:\n")
print(colnames(gene_expression_renamed))


# ------------------------------------------------------------
# 1.c Split gene expression into tumor vs normal
# ------------------------------------------------------------

# Tumor columns + Probe_ID
tumor_data <- gene_expression_renamed %>%
  select(Probe_ID, starts_with("tumor"))

# Normal columns + Probe_ID
normal_data <- gene_expression_renamed %>%
  select(Probe_ID, starts_with("normal"))

# ------------------ QC CHECKS -------------------------------
cat("QC: Tumor Data Dimensions:\n")
print(dim(tumor_data))

cat("\nQC: Normal Data Dimensions:\n")
print(dim(normal_data))

cat("\nTumor Columns:\n")
print(colnames(tumor_data))

cat("\nNormal Columns:\n")
print(colnames(normal_data))

# ------------------------------------------------------------
# 1.d Compute average expression values for each gene
# ------------------------------------------------------------

# Average tumor expression
tumor_avg <- tumor_data %>%
  mutate(tumor_mean = rowMeans(select(., starts_with("tumor")), na.rm = TRUE)) %>%
  select(Probe_ID, tumor_mean)

# Average normal expression
normal_avg <- normal_data %>%
  mutate(normal_mean = rowMeans(select(., starts_with("normal")), na.rm = TRUE)) %>%
  select(Probe_ID, normal_mean)

# ------------------ QC CHECKS -------------------------------
cat("QC: Tumor Avg Dimensions:\n")
print(dim(tumor_avg))

cat("\nQC: Normal Avg Dimensions:\n")
print(dim(normal_avg))

cat("\nTumor Avg Preview:\n")
print(head(tumor_avg))

cat("\nNormal Avg Preview:\n")
print(head(normal_avg))

# ------------------------------------------------------------
# 1.e Compute log2 fold change using the assignment formula
# log2((Tumor - Control) / Control)
# ------------------------------------------------------------

# Merge tumor and normal averages by Probe_ID
merged_avg <- tumor_avg %>%
  inner_join(normal_avg, by = "Probe_ID")

# Compute log2 fold change (NaNs expected when tumor < normal)
fold_change <- merged_avg %>%
  mutate(log2FC = log2((tumor_mean - normal_mean) / normal_mean))

# OPTIONAL: Replace NaN with NA (not required, but cleaner)
# fold_change$log2FC[is.nan(fold_change$log2FC)] <- NA

# ------------------ QC CHECKS -------------------------------
cat("QC: Fold Change Dimensions:\n")
print(dim(fold_change))

cat("\nFold Change Preview:\n")
print(head(fold_change))

cat("\nExplanation: NaNs occur when tumor_mean < normal_mean, making the log undefined.\n")

# ------------------------------------------------------------
# 1.f Identify genes with |log2FC| > 5 and merge gene info
# ------------------------------------------------------------

# Remove NaN values (cannot compare NaN > 5)
fold_change_clean <- fold_change %>%
  filter(!is.nan(log2FC))

# Filter genes with absolute fold change > 5
fc_filtered <- fold_change_clean %>%
  filter(abs(log2FC) > 5)

# Merge with gene information
fc_annotated <- fc_filtered %>%
  inner_join(gene_info, by = "Probe_ID")

# ------------------ QC CHECKS -------------------------------
cat("QC: Number of genes with |log2FC| > 5:\n")
print(nrow(fc_annotated))

cat("\nPreview of annotated significant genes:\n")
print(head(fc_annotated))

# ------------------------------------------------------------
# 1.g Add a column indicating whether gene is higher in Tumor or Normal
# ------------------------------------------------------------

fc_annotated_direction <- fc_annotated %>%
  mutate(Expression_High_in = case_when(
    tumor_mean > normal_mean ~ "Tumor",
    normal_mean > tumor_mean ~ "Normal",
    TRUE ~ "Equal"
  ))

# ------------------ QC CHECKS -------------------------------
cat("QC: Dimensions of annotated data with direction column:\n")
print(dim(fc_annotated_direction))

cat("\nPreview:\n")
print(head(fc_annotated_direction))


# Keep only valid human chromosomes
valid_chromosomes <- c(as.character(1:22), "X", "Y", "MT")

fc_annotated_direction_filtered <- fc_annotated_direction_clean %>%
  filter(Chromosome %in% valid_chromosomes)

# ============================================================
# 2.c Histogram of DEGs by chromosome segregated by sample type
# ============================================================

library(ggplot2)
library(dplyr)

# Count DEGs by chromosome AND direction (Tumor vs Normal)
deg_by_chr_type <- fc_annotated_direction_filtered %>%
  count(Chromosome, Expression_High_in)

# QC: Print counts
cat("QC: DEGs per chromosome by sample type:\n")
print(deg_by_chr_type)

# Grouped bar plot
ggplot(deg_by_chr_type, aes(x = Chromosome, y = n, fill = Expression_High_in)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(
    title = "DEGs by Chromosome Separated by Sample Type",
    x = "Chromosome",
    y = "Number of DEGs",
    fill = "Higher in"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Count DEGs per valid chromosome
deg_by_chr_filtered <- fc_annotated_direction_filtered %>%
  count(Chromosome)

# QC: Print filtered counts
cat("QC: Number of DEGs per REAL chromosome:\n")
print(deg_by_chr_filtered)

# Histogram (bar plot)
ggplot(deg_by_chr_filtered, aes(x = Chromosome, y = n)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "DEGs by Chromosome (Valid Human Chromosomes Only)",
    x = "Chromosome",
    y = "Number of DEGs"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# 2.d Bar Chart of Upregulated vs Downregulated DEGs
# ============================================================

library(dplyr)
library(ggplot2)

# Count DEGs by category
deg_percent <- fc_annotated_direction %>%
  count(Expression_High_in) %>%
  mutate(
    Percentage = round(100 * n / sum(n), 2)
  )

# QC: Show percentages
cat("QC: DEG Percentages:\n")
print(deg_percent)

# Bar plot of percentages
ggplot(deg_percent, aes(x = Expression_High_in, y = Percentage, fill = Expression_High_in)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Percentage of DEGs Upregulated vs Downregulated in Tumor",
    x = "DEG Category",
    y = "Percentage (%)",
    fill = "Higher In"
  ) +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Tumor" = "red", "Normal" = "blue"))

library(dplyr)
library(pheatmap)
library(tidyr)

# ============================================================
# 2.e Heatmap using TOP 500 most variable genes (raw data)
# ============================================================

# Use the renamed raw expression matrix from part 1b
expr_raw <- gene_expression_renamed

# Convert Probe_ID to rownames
expr_matrix <- expr_raw %>%
  column_to_rownames("Probe_ID") %>%
  as.matrix()

# ------------------------------------------------------------
# 1. Compute variance for each gene
# ------------------------------------------------------------
gene_variances <- apply(expr_matrix, 1, var)

# Select top 500 most variable genes
top500_genes <- names(sort(gene_variances, decreasing = TRUE))[1:500]

expr_top500 <- expr_matrix[top500_genes, ]

# ------------------------------------------------------------
# 2. Create sample annotation (Tumor vs Normal)
# ------------------------------------------------------------
sample_type <- ifelse(grepl("tumor", colnames(expr_top500)), "Tumor", "Normal")

annotation_col <- data.frame(
  SampleType = sample_type
)
rownames(annotation_col) <- colnames(expr_top500)

# ------------------------------------------------------------
# 3. Heatmap
# ------------------------------------------------------------
pheatmap(
  expr_top500,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Heatmap of Top 500 Most Variable Genes (Raw Expression)"
)



# ============================================================
# 2.f IMPROVED CLUSTERMAP with Blue–Green Color Palette
# ============================================================

library(pheatmap)
library(dplyr)

# ------------------ Create custom blue–green palette ------------------
blue_green_palette <- colorRampPalette(
  c("#08306B", "#2171B5", "#4DAF4A", "#A1D99B")
)(100)

# ------------------ Improved Clustermap -------------------------------
pheatmap(
  expr_top500,
  scale = "row",                        # Z-score normalization
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  annotation_col = annotation_col,
  color = blue_green_palette,           # <--- CUSTOM PALETTE
  show_rownames = FALSE,
  main = "Improved Clustermap (Blue–Green Palette)"
)






# ============================================================
# 2.g SHORT SUMMARY OF FINDINGS
# ============================================================

# - Tumor and normal samples showed clear differences in average gene expression.
#
# - Several genes had extremely high |log2FC| values using the assignment-specified
#   formula log2((Tumor - Normal) / Normal), indicating strong differential expression.
#
# - Chromosome-level plots showed that DEGs are distributed across many chromosomes,
#   with some (e.g., Chr 1, 11, 12, 17) containing more DEGs.
#
# - Both tumor-upregulated and normal-upregulated genes were identified.
#
# - The heatmap (top 500 most variable genes) showed clear separation between
#   tumor and normal samples based on global expression patterns.
#
# - The improved clustermap (using correlation distance, Ward.D2 clustering,
#   and a blue–green color palette) further highlighted strong grouping of
#   tumor vs. normal samples.
#
# Overall, the analysis confirms significant transcriptional differences
# between tumor and normal tissues.
