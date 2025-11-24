# B573_Assignment 5  
**Name**: Deeksha Kayyari  
**Programming Language**: R  
**Date**:  11/23/2025


## Description
This project performs a multi-step analysis of gene expression data using R. It integrates expression data, gene annotation, and sample phenotype information to compute fold changes, identify differentially expressed genes (DEGs), and generate exploratory visualizations, including heatmaps and clustermaps. The project applies data wrangling, 
statistical calculations, and bioinformatics-focused visualization to evaluate gene expression differences between tumor and normal samples.
## Required Files
Gene_Expression_Data.xlsx — Raw gene expression matrix

Gene_Information.csv — Gene metadata (symbols, loci, chromosome, etc.)

Sample_Information.tsv — Sample phenotype labels

R_SCRIPT_ADVANCED.R — R script containing all analysis steps

## Required packages
readxl – Load Excel data

readr – Load CSV/TSV data

dplyr – Data wrangling

tidyr – Data cleaning

stringr – String operations

ggplot2 – Visualization

pheatmap – Heatmap + clustermap visualization

## steps for execution
1.Run in RStudio

2.Open RStudio

3.Place all required files in the same directory

4.Open the script R_SCRIPT_ADVANCED.R

5.Set working directory:
Ctrl + Shift + Enter
or run each section individually

 ## output summary

 -Running the script produces the following outputs:

-Data Processing

-Cleaned sample phenotype table

-Renamed columns in gene expression matrix based on tumor/normal status

-Split expression data into tumor and normal groups

-Average expression per probe calculated for both groups

-Statistical Analysis

-Fold change computed using assignment formula:
log2((Tumor – Normal) / Normal)

-Identification of DEGs using the cutoff:
|log₂FC| > 5

-Annotation of DEGs with gene metadata

-Determination of whether expression is higher in Tumor or Normal samples

## Visualizations

-Histogram of DEGs by chromosome

-Histogram of DEGs by chromosome split by Tumor/Normal category

-Bar chart showing percentage of DEGs upregulated in Tumor vs Normal

-Heatmap of the Top 500 variable genes

-Clustermap (with correlation distance + Ward.D2) showing sample clustering

-Clear separation observed between tumor and normal profiles

 ### Notes
 -Phenotype-based renaming ensures clear Tumor vs Normal identification.

-Fold change formula follows the assignment requirement:

-log2((Tumour – Control) / Control)

-DEGs filtered with abs(log₂FC) > 5, as instructed.

-Heatmap uses the top 500 most variable genes to improve readability.

-Clustermap uses correlation distance, Ward.D2 clustering, and a blue–green palette, providing a more detailed and visually distinct structure compared to the heatmap.

-The analysis demonstrates significant transcriptomic differences between tumor and normal samples.
