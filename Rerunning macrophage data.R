#1. Install packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "clusterProfiler", "org.Mm.eg.db"))
options(repos = BiocManager::repositories())

#2. Load the libraries
library(clusterProfiler)
library(org.Mm.eg.db)  

if (!requireNamespace("conflicted", quietly = TRUE)) install.packages("conflicted")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::mutate)

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

#Ensure that the count data and metadata table information match
if (!identical(sort(sample_names), sort(rownames(metadata)))) {
  stop("Error: Mismatch between sample names in count data and metadata!")
}

####################################################################################################################################################################################

# # Define filtering thresholds
min_count <- 10
min_samples <- 0.5 * ncol(raw_counts)

# Apply filtering
keep <- rowSums(raw_counts >= min_count) >= min_samples
filtered_counts <- raw_counts[keep, ]
message("After filtering, ", nrow(filtered_counts), " genes remain.")

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

# Plot dispersion estimates
plotDispEsts(dds)


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



###########################################################################################
##Looping KEGG and GO analysis for all comparisons
####################################################################################################
### — Batch GO & KEGG for all comparisons — ###

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# Define comparisons and thresholds
comparisons <- c("WN_WI", "CN_CI", "WI_CI", "WN_CN")
padj.cut   <- 0.05
lfc.cut    <- 1
out.dir    <- "GO_KEGG_plots"
dir.create(out.dir, showWarnings = FALSE)

for (cmp in comparisons) {
  # 1) grab results object
  res.obj <- get(paste0("res_", cmp))
  
  # 2) filter DEGs
  deg     <- res.obj[!is.na(res.obj$padj) & res.obj$padj < padj.cut & abs(res.obj$log2FoldChange) > lfc.cut, ]
  genes   <- rownames(deg)
  
  # 3) map to Entrez
  gene.df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  # 4a) KEGG enrichment
  kegg.enr <- enrichKEGG(
    gene         = gene.df$ENTREZID,
    organism     = "mmu",
    pvalueCutoff = padj.cut
  )
  # 4b) GO enrichment (BP only)
  go.enr <- enrichGO(
    gene          = gene.df$ENTREZID,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = padj.cut,
    qvalueCutoff  = 0.2
  )
  
  # 5) plotting + titles
  p_kegg_dot <- dotplot(kegg.enr) + ggtitle(paste(cmp, "KEGG dotplot"))
  p_kegg_bar <- barplot(kegg.enr) + ggtitle(paste(cmp, "KEGG barplot"))
  p_go_dot   <- dotplot(go.enr, showCategory=10) + ggtitle(paste(cmp, "GO dotplot"))
  p_go_bar   <- barplot(go.enr, showCategory=10) + ggtitle(paste(cmp, "GO barplot"))
  
  # 6) save
  ggsave(file.path(out.dir, paste0(cmp, "_KEGG_dotplot.png")), plot=p_kegg_dot, width=8, height=6, dpi=300)
  ggsave(file.path(out.dir, paste0(cmp, "_KEGG_barplot.png")), plot=p_kegg_bar, width=8, height=6, dpi=300)
  ggsave(file.path(out.dir, paste0(cmp, "_GO_dotplot.png")),   plot=p_go_dot,   width=8, height=6, dpi=300)
  ggsave(file.path(out.dir, paste0(cmp, "_GO_barplot.png")),   plot=p_go_bar,   width=8, height=6, dpi=300)
}


############################################################################################
##Volcano PLots
# Install necessary packages if not already installed
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggrepel)) install.packages("ggrepel")  # Optional, for better labels

# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)


# Read your DEG data
deg <- read.csv("DEG_WI_CI.csv", header = TRUE, stringsAsFactors = FALSE)
names(deg)[1] <- "Gene"  # Rename the gene column for clarity

# Create a new column "Direction" to classify genes:
# Upregulated: padj < 0.05, log2FoldChange > 1
# Downregulated: padj < 0.05, log2FoldChange < -1
# Not Significant: otherwise
deg <- deg %>% mutate(Direction = case_when(
  padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
  padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
  TRUE ~ "Not Significant"
))

# Select the top 15 significant genes (based on lowest padj) for labeling
top15 <- deg %>% 
  filter(Direction != "Not Significant") %>% 
  arrange(padj) %>% 
  head(15)

# Create the volcano plot
volcano_plot <- ggplot(deg, aes(x = log2FoldChange, y = -log10(pvalue), color = Direction)) +
  geom_point(alpha = 0.8, size = 1.5) +
  # Assign colors: red for upregulated, blue for downregulated, grey for non-significant
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_bw() +
  xlab("Log₂ Fold Change") +
  ylab("-Log₁₀ p-value") +
  ggtitle("Differential Gene Expression of WT v CR Infected Macrophages") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # Add labels for the top 15 genes only
  geom_text_repel(data = top15,
                  aes(label = Gene),
                  size = 3,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  )

# Display the plot
print(volcano_plot)

# Save the plot to a file (e.g., PNG format)
ggsave("volcano_plot_colored.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)


######Saving all files
#Looping for Reactome and Volcano Plots for all comparisons
############################################################################################################
##How to do a loop of volcano plots and Reactome for all comparisons 

# Load libraries
library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)


comparisons <- c("WI_CI","WN_CN","WN_WI","CN_CI")


# choose an output folder
output_dir <- rstudioapi::selectDirectory("Select folder to save all plots")


# Confirm where we’ll save everything
cat("👉 All plots will be written to:\n", output_dir, "\n\n")


for(comp in comparisons) {
  message("=== Running for comparison: ", comp, " ===")
  
  # 1) Read DEGs
  deg <- read.csv(paste0("DEG_", comp, ".csv"), stringsAsFactors=FALSE)
  names(deg)[1] <- "Gene"
  
  # 2) Reactome enrichment on all sig. DEGs
  sig <- deg[deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, ]
  entrez <- bitr(sig$Gene, "SYMBOL","ENTREZID", OrgDb=org.Mm.eg.db)
  react_res <- enrichPathway(entrez$ENTREZID, organism="mouse", pvalueCutoff=0.05, readable=TRUE)
  
  # 3) Build & title the dotplot
  p_react <- dotplot(react_res) +
    ggtitle(paste("Reactome pathways —", comp)) +
    theme(plot.title=element_text(hjust=0.5))
  print(p_react)
  
  # 4) Save it as PDF & PNG
  ggsave(
    filename = file.path(output_dir, paste0("reactome_dotplot_", comp, ".pdf")),
    plot     = p_react,
    width    = 8, height = 6
  )
  ggsave(
    filename = file.path(output_dir, paste0("reactome_dotplot_", comp, ".png")),
    plot     = p_react,
    width    = 8, height = 6, dpi = 150
  )
  
  # 5) Volcano plot
  deg <- deg %>% 
    mutate(Direction = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  top15 <- deg %>% filter(Direction!="Not Significant") %>% arrange(padj) %>% head(15)
  
  p_volc <- ggplot(deg, aes(log2FoldChange, -log10(pvalue), color=Direction)) +
    geom_point(alpha=0.8, size=1.5) +
    scale_color_manual(values=c("Upregulated"="red","Downregulated"="blue","Not Significant"="grey")) +
    theme_bw() +
    xlab("Log₂ Fold Change") + ylab("-Log₁₀ p-value") +
    ggtitle(paste("Volcano plot —", comp)) +
    theme(plot.title=element_text(hjust=0.5)) +
    geom_text_repel(data=top15, aes(label=Gene), size=3,
                    box.padding=0.3, point.padding=0.3)
  print(p_volc)
  
  ggsave(
    filename = file.path(output_dir, paste0("volcano_plot_", comp, ".png")),
    plot     = p_volc,
    width    = 8, height = 6, dpi = 300
  )
}

# … after you build p_react …
print(p_react)      # sends it to the active graphics device (e.g. RStudio Plots pane)

# … after you build p_volc …
print(p_volc)       # same for the volcano plot

cat("\n✅ Done! Files in output directory:\n")
print(list.files(output_dir, pattern="\\.(png|pdf)$"))







###Individual running for comparisons ###############

######################################################################################################
##KEGG and GO individal runs for comparisons 
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
                           organism     = 'mmu',   
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

### — Save all generated KEGG & GO plots — ###

# Re‑compute and collect plots into a named list
plots <- list(
  kegg_all_dot    = dotplot(kegg_results),
  kegg_all_bar    = barplot(kegg_results),
  
  WN_WI_kegg_dot  = dotplot(kegg_WN_WI),
  WN_WI_kegg_bar  = barplot(kegg_WN_WI),
  
  CN_CI_kegg_dot  = dotplot(kegg_CN_CI),
  CN_CI_kegg_bar  = barplot(kegg_CN_CI),
  
  WI_CI_kegg_dot  = dotplot(kegg_WI_CI),
  WI_CI_kegg_bar  = barplot(kegg_WI_CI),
  
  WN_CN_kegg_dot  = dotplot(kegg_WN_CN),
  WN_CN_kegg_bar  = barplot(kegg_WN_CN),
  
  WN_WI_go_dot    = dotplot(go_WN_WI, showCategory = 10),
  WN_WI_go_bar    = barplot(go_WN_WI, showCategory = 10),
  
  CN_CI_go_dot    = dotplot(go_CN_CI, showCategory = 10),
  CN_CI_go_bar    = barplot(go_CN_CI, showCategory = 10),
  
  WI_CI_go_dot    = dotplot(go_WI_CI, showCategory = 10),
  WI_CI_go_bar    = barplot(go_WI_CI, showCategory = 10),
  
  WN_CN_go_dot    = dotplot(go_WN_CN, showCategory = 10),
  WN_CN_go_bar    = barplot(go_WN_CN, showCategory = 10)
)

# Loop over the list and save each plot
for (nm in names(plots)) {
  ggsave(
    filename = paste0(nm, ".png"),
    plot     = plots[[nm]],
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}

##############################################################################################################
#Running the Reactome one at a time without looping
#########################################################################

##Reactome pathway analysis Do the rest of this code individually for each comparison

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
# Alternatively, exporting as a PNG file
png("reactome_dotplot.png", width = 800, height = 600, res = 150)
dotplot(reactome_up)
dev.off()


##############################################################################################################
##How to look at genes responsible for this

# Convert the enrichment results to a data frame
reactome_df <- as.data.frame(reactome_down)

# Iterate over each enriched pathway and print the pathway description and its genes
for(i in 1:nrow(reactome_df)) {
  cat("Pathway:", reactome_df$Description[i], "\n")
  # Split the geneID string to obtain individual gene symbols
  pathway_genes <- strsplit(reactome_df$geneID[i], "/")[[1]]
  cat("Genes:", paste(pathway_genes, collapse = ", "), "\n\n")
}

# Alternatively, store the genes in a list for further analysis
pathway_genes_list <- lapply(reactome_df$geneID, function(x) strsplit(x, "/")[[1]])
names(pathway_genes_list) <- reactome_df$Description

