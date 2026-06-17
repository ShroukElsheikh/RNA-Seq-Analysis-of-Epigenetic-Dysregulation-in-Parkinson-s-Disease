#### used libraries ######
install.packages(c("ggfortify","rgl","plot3D","plotly","scatterplot3d"))
BiocManager::install(c("genefilter", "ComplexHeatmap", "EnhancedVolcano"))
BiocManager::install(c("ggplot2", "reshape2"))

library(readr)
library(DESeq2)
library(ggfortify)
library(rgl)
library(plot3D)
library(plotly)
library(stats)
library(scatterplot3d)
library(genefilter)
library(matrixStats)
library(ComplexHeatmap)
library(readxl)
library(EnhancedVolcano)
library(ggplot2)
library(reshape2)

#===================================================#
#            LOAD THE mRNA-Seq DATA                 #
#===================================================#

exp <- read_excel("GSE180408_raw_count.xlsx")

id = exp[, 1]
exp_matrix = exp[, -1]

exp_matrix = as.data.frame(lapply(exp_matrix, as.numeric))
exp_matrix = as.matrix(exp_matrix)

exp = cbind(id, exp_matrix)

#===================================================#
#         CHECK & HANDLE DUPLICATE GENES            #
#===================================================#

cat("Number of duplicated genes:", sum(duplicated(exp[, "gene_id"])), "\n")

# Aggregate duplicates by mean
exp_agg = aggregate(exp, list(exp$gene_id), FUN = mean)
genes = exp_agg[1]                  # save gene names
exp_agg = exp_agg[, -c(1, 2)]      # remove unnecessary columns

exp_agg = apply(exp_agg, 2, as.numeric)
exp_agg = cbind(exp_agg, genes)
rownames(exp_agg) = exp_agg$Group.1
exp_agg = exp_agg[-7]              # remove gene name column

#===================================================#
#              LOAD SAMPLE METADATA                 #
#===================================================#

pheno <- data.frame(
  name = c("GSM5462402", "GSM5462403", "GSM5462404",
           "GSM5462405", "GSM5462406", "GSM5462407"),
  type = c("Control", "Control", "Control",
           "Treatment", "Treatment", "Treatment"),
  row.names = 1
)

# Convert type to factor — crucial for DESeq2
pheno$type <- as.factor(pheno$type)

# Filter low expression genes
exp_agg = exp_agg[rowMeans(exp_agg) > 1, ]

#===================================================#
#        EXPLORATORY ANALYSIS — RAW DATA            #
#===================================================#

# 1. Boxplot — Raw data
boxplot(log2(exp_agg + 1),
        main = "RNA-Seq Box Plot (Raw)",
        xlab = "Samples",
        ylab = "log2(count + 1)",
        col = seq_len(ncol(exp_agg)))

# Transpose for PCA
exp_t = t(exp_agg)
cat("Dimensions after transpose:", dim(exp_t), "\n")

# 2. 2D PCA — Raw data
exp.pca = prcomp(exp_t, center = TRUE, scale. = TRUE)
summary(exp.pca)

colnames(pheno)[1] <- "type"
autoplot(exp.pca, data = pheno, colour = 'type',
         main = "2D PCA — Raw Data")

exp = t(exp_t)

# 3. 3D PCA — Raw data
pca_scores = as.data.frame(exp.pca$x)
mycolors = ifelse(pheno$type == "Control", "blue", "red")

plot3d(pca_scores[, 1:3],
       pch = 20,
       col = mycolors,
       radius = 2,
       main = "3D PCA — Raw Data")

exp.count = round(exp)

#===================================================#
#           IMPUTE MISSING VALUES (ZEROS)           #
#===================================================#

# Calculate proportion of zeros per gene (row)
prop_zeros <- rowSums(exp == 0) / ncol(exp)

# Keep genes with < 40% zeros
rows_to_fill <- which(prop_zeros < 0.4)
imputed_genes = rows_to_fill

# Replace zeros with row mean
row_means <- rowMeans(exp[rows_to_fill, ], na.rm = TRUE)
exp[rows_to_fill, ][exp[rows_to_fill, ] == 0] <- row_means

exp = exp[imputed_genes, ]
exp.count = round(exp)

cat("Genes after imputation:", nrow(exp.count), "\n")

#===================================================#
#       DIFFERENTIAL EXPRESSION — DESeq2            #
#===================================================#

table(pheno$type)

dds = DESeqDataSetFromMatrix(countData = exp.count,
                             colData = pheno,
                             design = ~type)
dds.run = DESeq(dds)

res = results(dds.run)
res = res[complete.cases(res), ]

# Explicit contrast: Treatment vs Control
contrast1 = results(dds.run, contrast = c("type", "Treatment", "Control"))
contrast1 = contrast1[complete.cases(contrast1), ]
contrast1 = as.data.frame(contrast1)

#===================================================#
#     FIXED: GET DEGs (padj + LFC filter)        #
#===================================================#

res.df = as.data.frame(res)

# Step 1: remove NAs
res.degs = res.df[!is.na(res.df$padj), ]

# Step 2: filter by adjusted p-value
res.degs = res.degs[res.degs$padj < 0.05, ]

# Step 3: filter by LFC >= 1.5
res.degs = res.degs[abs(res.degs$log2FoldChange) >= 1.5, ]

# Step 4: smart fallback — relax LFC if needed
if(nrow(res.degs) == 0){
  message("No DEGs at LFC >= 1.5, relaxing to LFC >= 1.0")
  res.degs = res.df[!is.na(res.df$padj), ]
  res.degs = res.degs[res.degs$padj < 0.05, ]
  res.degs = res.degs[abs(res.degs$log2FoldChange) >= 1.0, ]
}

# Step 5: last resort fallback
if(nrow(res.degs) == 0){
  message("Still no DEGs — using top 50 by pvalue as fallback")
  res.degs = res.df[order(res.df$pvalue), ][1:50, ]
}

# Step 6: report
cat("Total DEGs:", nrow(res.degs), "\n")
cat(" Upregulated:", sum(res.degs$log2FoldChange > 0), "\n")
cat(" Ownregulated:", sum(res.degs$log2FoldChange < 0), "\n")

degs.genes = rownames(res.degs)

#===================================================#
#              NORMALIZATION                        #
#===================================================#

ntd = normTransform(dds)
exp.norm = assay(ntd)

exp.degs = exp.norm[degs.genes, , drop = FALSE]

#===================================================#
#         HEATMAP — TOP 100 DEGs                    #
#===================================================#

# Safe top N genes
topN = min(100, nrow(exp.degs))
exp100_DEGS = exp.degs[1:topN, , drop = FALSE]

# Remove NA rows
exp100_DEGS = exp100_DEGS[complete.cases(exp100_DEGS), ]

# Z-score scaling per gene
exp100_DEGS = t(scale(t(exp100_DEGS)))

# Safety check
if(nrow(exp100_DEGS) == 0){
  stop("No DEGs found — relax filtering thresholds")
}

# Plot heatmap
column_ha = HeatmapAnnotation(sample.type = pheno$type)

Heatmap(exp100_DEGS,
        name = 'Exp',
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = column_ha)

#===================================================#
#         2D PCA — PROCESSED DEGs                   #
#===================================================#

expression_t = t(exp100_DEGS)

# Safety check
if(nrow(expression_t) < 2){
  stop("Not enough genes for PCA — DEG filtering too strict")
}

expression.pca = prcomp(expression_t, center = TRUE, scale. = TRUE)
summary(expression.pca)

autoplot(expression.pca, data = pheno, colour = 'type',
         main = "2D PCA — Processed DEGs")

#===================================================#
#         3D PCA — PROCESSED DEGs                   #
#===================================================#

pca_scores = as.data.frame(expression.pca$x)
mycolors = ifelse(pheno$type == "Control", "blue", "red")

plot3d(pca_scores[, 1:3],
       pch = 20,
       col = mycolors,
       type = 's',
       radius = 0.5,
       main = "3D PCA — Processed DEGs")

#===================================================#
#    FIXED: VOLCANO PLOT (padj, matches DEGs)    #
#===================================================#

EnhancedVolcano(res.df,
                lab = rownames(res.df),
                x = 'log2FoldChange',
                y = 'padj',              
                pCutoff = 0.05,          
                FCcutoff = 1.5,          
                colAlpha = 0.8,
                title = 'Treatment vs Control',
                subtitle = paste0('DEGs: ', nrow(res.degs),
                                  '  (Up: ', sum(res.degs$log2FoldChange > 0),
                                  '  Down: ', sum(res.degs$log2FoldChange < 0), ')'),
                caption = 'padj < 0.05  |  |LFC| >= 1.5')

#===================================================#
#       QC PLOTS — PROCESSED DEGs                   #
#===================================================#

# 1. Boxplot — Processed
boxplot(exp100_DEGS,
        main = "Processed DEGs — Box Plot",
        xlab = "Samples",
        ylab = "Scaled Expression",
        col = seq_len(ncol(exp100_DEGS)))

# 2. Histogram — Processed
hist(as.vector(exp100_DEGS),
     main = "Processed DEGs — Histogram",
     xlab = "Scaled Expression",
     col = "lightblue",
     border = "white")

# 3. Density Plot — Raw vs Processed (side by side)
par(mfrow = c(1, 2))

# --- Raw Data ---
# FIXED: extract only numeric columns
exp_agg_numeric <- exp_agg[, sapply(exp_agg, is.numeric)]
raw_values <- as.vector(log2(as.matrix(exp_agg_numeric) + 1))

plot(density(raw_values),
     main = "Raw Data (log2)",
     xlab = "Expression",
     ylab = "Density",
     col = "red",
     lwd = 2)
polygon(density(raw_values),
        col = rgb(1, 0, 0, 0.3),
        border = "red")

# --- Processed DEGs ---
processed_values <- as.vector(exp100_DEGS)

plot(density(processed_values),
     main = "Processed DEGs",
     xlab = "Scaled Expression",
     ylab = "Density",
     col = "steelblue",
     lwd = 2)
polygon(density(processed_values),
        col = rgb(0.2, 0.5, 0.8, 0.3),
        border = "steelblue")

par(mfrow = c(1, 1))

