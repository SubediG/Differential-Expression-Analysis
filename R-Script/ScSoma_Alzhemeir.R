# Loading required libraries
library(here)
library(DESeq2)
library(Matrix)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)


# Loading dataset
file <- readRDS("ScSoma_Alzhemeir.rds")


# Extracting raw counts and filter low-expression genes
raw_counts <- GetAssayData(file, layer = "counts")


#Histogram to identify target threshold for raw counts
hist(rowSums(raw_counts), breaks = 100, main = "Distribution of Gene Counts", xlab = "Counts per Gene")


# Defining filtering thresholds according to histogram as most count is 0 and setting count 5 as min across 10% of the sample.
min_counts <- 5
min_samples <- floor(0.10 * ncol(raw_counts)) 
filtered_counts <- raw_counts[rowSums(raw_counts > min_counts) >= min_samples, ]
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0, ] 
cat("Filtered genes from", nrow(raw_counts), "to", nrow(filtered_counts), "based on low-expression filtering.\n")


# Extracting and preparing metadata
metadata <- file@meta.data
metadata$disease <- gsub("[^a-zA-Z0-9_.]", "_", metadata$disease)
metadata$disease <- as.factor(metadata$disease)

if (!all(colnames(filtered_counts) %in% rownames(metadata))) {
  stop("Mismatch between counts column names and metadata row names.")
}



# Validating library size ratio for normalization strategy
library_sizes <- colSums(filtered_counts)
library_size_ratio <- max(library_sizes) / min(library_sizes)
if (library_size_ratio > 3) {
  cat("Library size ratio exceeds 3-fold. Consider using voom-limma for better normalization.\n")
}



# Creating DGE list
dge <- DGEList(counts = filtered_counts)
dge$targets <- metadata
dge <- calcNormFactors(dge, method = "TMM")


# Design matrix and setting reference
metadata$disease <- relevel(metadata$disease, ref = "normal")
design <- model.matrix(~disease, data = metadata)


# Applying voom transformation
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)


# Extracting Differentially Expressed Genes(DEGs) with thresholds
results <- topTable(fit, coef = 2, number = Inf)
DEG_thresholds <- list(adj_p_val = 0.05, logFC = 1)
DEGs <- results[results$adj.P.Val < DEG_thresholds$adj_p_val & abs(results$logFC) > DEG_thresholds$logFC, ]


# Printing DEGs summary
cat("Identified", nrow(DEGs), "differentially expressed genes (DEGs).\n")


# Preparing data for volcano plot
results$Siinggnificance <- ifelse(results$adj.P.Val < DEG_thresholds$adj_p_val & 
                                 abs(results$logFC) > DEG_thresholds$logFC, 
                               "Significant", "Not Significant")
top_genes <- results[results$Significance == "Significant", ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ][1:min(10, nrow(top_genes)), ]


# Creating volcano plot with volcano plot
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.6, size = ifelse(results$Significance == "Significant", 2.5, 1)) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3.5, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log Fold Change", y = "-log10(P-value)") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

print(volcano_plot)


# Saving results to CSV
output_file <- here("DEGs_results.csv")
write.csv(DEGs, output_file, row.names = TRUE)
cat("DEG results saved to", output_file, "\n")
