#1. Install packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "clusterProfiler", "org.Mm.eg.db"))
options(repos = BiocManager::repositories())

#2. Load the libraries
library(DESeq2)
library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)  

########################################Create metadata table

library(readxl)
# 3. Read gene count data from Excel
file_path <- file.choose(gene_count) ##The file choose will change with each analysis
counts_data <- read_excel(file_path, sheet = 1)
head(counts_data)

##Get it to pull from gene name instead of gene ID
gene_id <- counts_data$gene_id     
gene_names <- counts_data$gene_name 

# 4. Extract raw count data columns (assume column 1 = gene IDs, columns 2:24 = count data)
raw_counts <- counts_data[, 2:24]
#Make sure no duplicate genes
if(any(duplicated(gene_names))) {
  gene_names <- make.unique(as.character(gene_names))
}

#5. Convert the count data to numeric (replace "-" with 0 if appropriate)
raw_counts <- as.data.frame(lapply(raw_counts, function(x) as.numeric(gsub("-", "0", x))))
rownames(raw_counts) <- gene_names  # Set gene names as rownames

##Ensure dimensions are aligned
length(gene_names)
nrow(raw_counts)

#6. Check the conversion - make sure values are not 0
summary(raw_counts)
min(as.matrix(raw_counts))

#7. Create a metadata table; adjust the condition vector based on your experimental design
sample_names <- colnames(counts_data)[2:24]
metadata <- data.frame(
  sample = sample_names,
  condition = c("WN", "WN", "CN", "CN", "WI", "WI", "CI", "CI",
                "WN", "WN", "CN", "CN", "WI", "WI", "CI", "CI",
                "WN", "WN", "CN", "CN", "WI", "WI", "CI")
)
rownames(metadata) <- metadata$sample

#Ensure dimenions of metadata table match
length(sample_names)  # This should return however many columns you have for your data
nrow(metadata)         # should equal number_of_samples

# Preview the metadata table
print(metadata)

####################################################################################################################################################################################

# 9. Filter genes: throws out genes that don't reach count of 10 in at least half of the samples
keep <- rowSums(raw_counts >= 10) >= (0.5 * ncol(raw_counts))
filtered_counts <- raw_counts[keep, ]
dim(filtered_counts)  # Check dimensions after filtering

##DESEq########################################
#Ensure all are factors and not characters 
metadata$condition <- factor(metadata$condition)

library(DESeq2)
# 10. Create the DESeq2 dataset using the filtered counts and metadata
dds <- DESeqDataSetFromMatrix(countData = as.matrix(filtered_counts),
                              colData = metadata,
                              design = ~ condition)
dds <- DESeq(dds)

###Setting Comparisons##
## DESeq2 uses negative binomial model with wALD Test to compare conditions
####Log2 Fold Change (M): How much the gene’s expression changes between WN and WI. P-value and Adjusted P-value (padj): The statistical significance of that change. By default, DESeq2 flags genes with an adjusted p-value below 0.1 as significant, although you can apply stricter cutoffs (e.g., padj < 0.05) when filtering the results.

res_WN_WI <- results(dds, contrast = c("condition", "WN", "WI"))
summary(res_WN_WI)
res_CN_CI <- results(dds, contrast = c("condition", "CN", "CI"))
summary(res_CN_CI)
res_WI_CI <- results(dds, contrast = c("condition", "WI", "CI"))
summary(res_WI_CI)
res_WN_CN <- results(dds, contrast = c("condition", "WN", "CN"))
summary(res_WN_CN)

# Export DESeq2 results for comparisons ith the clean data and the group labels, the code uses a tool called DESeq2. This tool compares the groups (like WN vs. WI, CN vs. CI, etc.) to see which genes are turned on or off differently between them.
write.csv(as.data.frame(res_WN_WI), file = "DEG_WN_WI.csv", row.names = TRUE)
write.csv(as.data.frame(res_CN_CI), file = "DEG_CN_CI.csv", row.names = TRUE)
write.csv(as.data.frame(res_WN_CN), file = "DEG_WN_CN.csv", row.names = TRUE)
write.csv(as.data.frame(res_WI_CI), file = "DEG_WI_CI.csv", row.names = TRUE)

####MA Plots#### automatically highlights genes considered significant (typically those with an adjusted p-value < 0.1). Genes outside this significance cutoff are shown in a default, muted color (usually grey), whereas significant genes are often highlighted in red.
plotMA(res_WN_WI, main = "MA Plot: WN vs WI", ylim = c(-2, 2))
plotMA(res_WI_CI, main = "MA Plot: WI vs CI", ylim = c(-2, 2))
plotMA(res_CN_CI, main = "MA Plot: CN vs CI", ylim = c(-2, 2))
plotMA(res_WN_CN, main = "MA Plot: WN vs CN", ylim = c(-2, 2))



################################################################################################################################################################
###KEGG Pathway Analysis###

# Convert gene symbols (gene names) to Entrez IDs
gene.df <- bitr(rownames(raw_counts), 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Mm.eg.db)

# Remove duplicate SYMBOL entries to retain one mapping per gene
gene.df <- gene.df[!duplicated(gene.df$SYMBOL), ]


# You can then proceed with your KEGG enrichment analysis using the converted IDs
##KEGG over ALL conditions
kegg_results <- enrichKEGG(gene         = gene.df$ENTREZID,
                           organism     = 'mmu',   # 'hsa' for human
                           pvalueCutoff = 0.05) #Pathways with p- value less than 0.05 considered significant

dotplot(kegg_results) #For all conditions 
barplot(kegg_results)

# 1. Filter DEGs from the WN vs WI comparison
# (e.g., adjusted p-value < 0.05) and, optionally, with a meaningful fold change (e.g., absolute log₂ fold change > 1). 
deg_WN_WI <- res_WN_WI[!is.na(res_WN_WI$padj) & res_WN_WI$padj < 0.05 & abs(res_WN_WI$log2FoldChange) > 1, ]
gene_symbols_WN_WI <- rownames(deg_WN_WI)

# 2. Convert filtered gene symbols to Entrez IDs (using org.Mm.eg.db for mouse)
gene.df_WN_WI <- bitr(gene_symbols_WN_WI,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

# 3. Perform KEGG enrichment analysis using the Entrez IDs
kegg_WN_WI <- enrichKEGG(gene         = gene.df_WN_WI$ENTREZID,
                         organism     = 'mmu',         # 'mmu' for mouse; use 'hsa' for human
                         pvalueCutoff = 0.05)

# 4. View the KEGG results
print(kegg_WN_WI)
dotplot(kegg_WN_WI)
barplot(kegg_WN_WI)

# Filter DEGs for CN vs CI comparison
deg_CN_CI <- res_CN_CI[!is.na(res_CN_CI$padj) & res_CN_CI$padj < 0.05 & abs(res_CN_CI$log2FoldChange) > 1, ]
gene_symbols_CN_CI <- rownames(deg_CN_CI)

# Convert gene symbols to Entrez IDs
gene.df_CN_CI <- bitr(gene_symbols_CN_CI,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

# Perform KEGG enrichment analysis
kegg_CN_CI <- enrichKEGG(gene         = gene.df_CN_CI$ENTREZID,
                         organism     = 'mmu',
                         pvalueCutoff = 0.05)

# View the results
print(kegg_CN_CI)
dotplot(kegg_CN_CI)
barplot(kegg_CN_CI)

# Filter DEGs for WI vs CI comparison
deg_WI_CI <- res_WI_CI[!is.na(res_WI_CI$padj) & res_WI_CI$padj < 0.05 & abs(res_WI_CI$log2FoldChange) > 1, ]
gene_symbols_WI_CI <- rownames(deg_WI_CI)

# Convert gene symbols to Entrez IDs
gene.df_WI_CI <- bitr(gene_symbols_WI_CI,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

# Perform KEGG enrichment analysis
kegg_WI_CI <- enrichKEGG(gene         = gene.df_WI_CI$ENTREZID,
                         organism     = 'mmu',
                         pvalueCutoff = 0.05)

# View the results
print(kegg_WI_CI)
dotplot(kegg_WI_CI)
barplot(kegg_WI_CI)

# Filter DEGs for WN vs CN comparison
deg_WN_CN <- res_WN_CN[!is.na(res_WN_CN$padj) & res_WN_CN$padj < 0.05 & abs(res_WN_CN$log2FoldChange) > 1, ]
gene_symbols_WN_CN <- rownames(deg_WN_CN)

# Convert gene symbols to Entrez IDs
gene.df_WN_CN <- bitr(gene_symbols_WN_CN,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

# Perform KEGG enrichment analysis
kegg_WN_CN <- enrichKEGG(gene         = gene.df_WN_CN$ENTREZID,
                         organism     = 'mmu',
                         pvalueCutoff = 0.05)

# View the results
print(kegg_WN_CN)
dotplot(kegg_WN_CN)
barplot(kegg_WN_CN)

##########################################################################################################################
# Function to create a volcano plot using base R
volcano_plot <- function(res, title, topN = 5) {
  df <- as.data.frame(res)
  # Create the base plot
  plot(df$log2FoldChange, -log10(df$padj),
       pch = 20, main = title,
       xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted P-value",
       col = "gray")
  
  
  # Identify significant points
  significant <- which(df$padj < 0.05 & abs(df$log2FoldChange) > 1)
  
  # Add colored points for significant genes
  points(df$log2FoldChange[significant],
         -log10(df$padj)[significant],
         pch = 20, col = "pink")
  
  # Only label the topN significant genes sorted by lowest padj
  if (length(significant) > 0) {
    ordered_sig <- significant[order(df$padj[significant])]
    top_genes <- head(ordered_sig, topN)
    
    text(df$log2FoldChange[top_genes],
         -log10(df$padj)[top_genes],
         labels = rownames(df)[top_genes],
         pos = 3, cex = 0.8, col = "purple")
  }
}



# Create volcano plots for each comparison
volcano_plot(res_WN_WI, "Volcano Plot: WN vs WI", topN = 10)
volcano_plot(res_CN_CI, "Volcano Plot: CN vs CI", topN = 10)
volcano_plot(res_WI_CI, "Volcano Plot: WI vs CI", topN = 10)
volcano_plot(res_WN_CN, "Volcano Plot: WN vs CN" , topN = 10)

###########################################################################################################################################################
##GO analysis 

go_WN_WI <- enrichGO(gene          = gene.df_WN_WI$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",         # Choose "BP", "CC", or "MF"
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)

# Visualize the GO results
dotplot(go_WN_WI, showCategory = 10)
barplot(go_WN_WI, showCategory = 10)

# ---- NEW CODE: GO Analysis for CN vs CI ----

go_CN_CI <- enrichGO(gene          = gene.df_CN_CI$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)
dotplot(go_CN_CI, showCategory = 10)
barplot(go_CN_CI, showCategory = 10)

# ---- NEW CODE: GO Analysis for WI vs CI ----

go_WI_CI <- enrichGO(gene          = gene.df_WI_CI$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)
dotplot(go_WI_CI, showCategory = 10)
barplot(go_WI_CI, showCategory = 10)

# ---- NEW CODE: GO Analysis for WN vs CN ----

go_WN_CN <- enrichGO(gene          = gene.df_WN_CN$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)
dotplot(go_WN_CN, showCategory = 10)
barplot(go_WN_CN, showCategory = 10)







##################################################################################################################################################
##Reactome pathway analysis

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

# Load libraries
library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)

# Load your DEG data (adjust the file path if needed)
deg <- read.csv("DEG_WI_CI.csv", header = TRUE, stringsAsFactors = FALSE)

# Optional: rename the gene column for clarity (assumes gene symbols are in the first column)
names(deg)[1] <- "GeneSymbol"

# Inspect the data
head(deg)

# Filter significant DEGs based on adjusted p-value and log2 fold change thresholds
# Adjust these thresholds as needed (here, padj < 0.05 and |log2FoldChange| > 1)
sig_deg <- deg[deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, ]

# Extract the gene symbols from the significant DEGs
gene_list <- sig_deg$GeneSymbol

# Convert gene symbols to Entrez IDs using the org.Mm.eg.db annotation database
gene_df <- bitr(gene_list,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# Run Reactome pathway enrichment analysis (for all significant DEGs)
reactome_results <- enrichPathway(gene         = gene_df$ENTREZID,
                                  organism     = "mouse",
                                  pvalueCutoff = 0.05,
                                  readable     = TRUE)

# View the top enriched pathways
head(reactome_results)

# Visualize the results using a dotplot
dotplot(reactome_results)

# Optionally, perform separate analyses for upregulated and downregulated genes:
# Upregulated genes only:
up_deg <- deg[deg$padj < 0.05 & deg$log2FoldChange > 1, ]
up_gene_list <- up_deg$GeneSymbol
up_gene_df <- bitr(up_gene_list,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
reactome_up <- enrichPathway(gene         = up_gene_df$ENTREZID,
                             organism     = "mouse",
                             pvalueCutoff = 0.05,
                             readable     = TRUE)
dotplot(reactome_up)

# Downregulated genes only:
down_deg <- deg[deg$padj < 0.05 & deg$log2FoldChange < -1, ]
down_gene_list <- down_deg$GeneSymbol
down_gene_df <- bitr(down_gene_list,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
reactome_down <- enrichPathway(gene         = down_gene_df$ENTREZID,
                               organism     = "mouse",
                               pvalueCutoff = 0.05,
                               readable     = TRUE)
dotplot(reactome_down)

# Exporting as a PDF file
pdf("reactome_dotplot.pdf", width = 8, height = 6)
dotplot(reactome_results)  # Your plotting function call
dev.off()  # Close the device and write the file

# Alternatively, exporting as a PNG file
png("reactome_dotplot.png", width = 800, height = 600, res = 150)
dotplot(reactome_results)
dev.off()




















































