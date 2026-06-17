#### used libraries ######
install.packages(c( "ggfortify","rgl","plot3D","plotly","scatterplot3d"))
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

##### load the mRNA-Seq data #####

exp <- read_excel("GSE180408_raw_count.xlsx")

id = exp[, 1]
exp_matrix = exp[, -1]

exp_matrix = as.data.frame(lapply(exp_matrix, as.numeric))
exp_matrix = as.matrix(exp_matrix)

exp = cbind(id , exp_matrix)


######## check duplicates of gene symbol #######
sum(duplicated(exp[, "gene_id"]))
######## aggregate of duplicated proteins ######3
exp_agg= aggregate(exp, list(exp$gene_id),FUN=mean)
genes = exp_agg[1] # save gene names
exp_agg = exp_agg[,-c(1,2)] # remove unnecessary columns

exp_agg = apply(exp_agg , 2, as.numeric) # apply the function to each column of the matrix
exp_agg = cbind(exp_agg , genes) # add the gene column back to the dataset
rownames(exp_agg) = exp_agg$Group.1 # set gene names as row names 
exp_agg = exp_agg[-7]

########### load mRNA sample sheet #######
# Create the metadata table
pheno <- data.frame(
  name = c("GSM5462402", "GSM5462403", "GSM5462404", "GSM5462405", "GSM5462406", "GSM5462407"),
  type = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
  row.names = 1 # This makes the 'name' column the row names
)

# Convert 'type' to a factor (Crucial for DESeq2)
pheno$type <- as.factor(pheno$type)

exp_agg = exp_agg[rowMeans(exp_agg)>1,]

############## exploratory analysis on raw expression matrix ###########
# 1. boxplot
boxplot(log2(exp_agg + 1),
        main = "RNA-Seq Box Plot",
        col = seq_len(ncol(exp_agg)))

# transpose data to prepare it to PCA
exp_t = t(exp_agg)
dim(exp_t)

# 2. 2D PCA
exp.pca = prcomp(exp_t, center = TRUE, scale. = TRUE)
summary(exp.pca)
colnames(pheno)[colnames(pheno) == "Sample_Id"] <- "type"
colnames(pheno)[1] <- "type"
autoplot(exp.pca, data = pheno, colour = 'type')

exp = t(exp_t)

# 3. 3D PCA
pca_scores = as.data.frame(exp.pca$x)

mycolors = ifelse(pheno$type == "Control", "blue", "red")

plot3d(pca_scores[,1:3],
       pch = 20,
       col = mycolors,
       radius = 2)

exp.count = round(exp)


############### impute the missing values mean ##############
# Calculate the proportion of zeros in each row
prop_zeros <- rowSums(exp == 0) / ncol(exp)
# Identify rows with less than 40% zeros
rows_to_fill <- which(prop_zeros < 0.4)
imputed_genes = rows_to_fill
# Calculate the row means for these rows
row_means <- rowMeans(exp[rows_to_fill, ], na.rm = TRUE)
# Replace the zeros in these rows with the row means and remove others rows
exp[rows_to_fill, ][exp[rows_to_fill, ] == 0] <- row_means
exp = exp[imputed_genes,]
exp.count = round(exp)

############### Differential expression analysis using DESeq2 #########
table(pheno$type)

dds = DESeqDataSetFromMatrix(countData = exp.count , colData = pheno , design = ~type)
dds.run = DESeq(dds)

res=results(dds.run)
res=res[complete.cases(res),]

########### make contrasts (comparisons between all conditions) ########
contrast1 = results(dds.run, contrast = c("type", "Treatment", "Control"))
contrast1 = contrast1[complete.cases(contrast1), ]
contrast1 = as.data.frame(contrast1)

############# get the DEGs based on adj pval , LFC ############
res.df = as.data.frame(res)

res.degs = res.df[!is.na(res.df$padj), ]

# filter by adjusted p-value
res.degs = res.degs[res.degs$padj < 0.05, ]

# filter by absolute LFC >= 1.5
res.degs = res.degs[abs(res.degs$log2FoldChange) >= 1.5, ]

# fallback if still empty
if(nrow(res.degs) == 0){
  res.degs = res.df[order(res.df$pvalue), ]
}

# fallback if still empty
if(nrow(res.degs) == 0){
  res.degs = res.df[order(res.df$pvalue), ]
}

degs.genes = rownames(res.degs)

############# do normalization for all exp data to further analysis ########
ntd = normTransform(dds)
exp.norm = assay(ntd)

exp.degs = exp.norm[degs.genes, , drop = FALSE]

############### creating a heatmap for the top 100 DEG genes #####
exp.degs = exp.norm[degs.genes, , drop = FALSE]

# SAFE top genes
topN = min(100, nrow(exp.degs))
exp100_DEGS = exp.degs[1:topN, , drop = FALSE]

# remove NA
exp100_DEGS = exp100_DEGS[complete.cases(exp100_DEGS), ]

# scaling
exp100_DEGS = t(scale(t(exp100_DEGS)))


# FIX: sample alignment
exp.degs = exp.norm[degs.genes, , drop = FALSE]

# prevent empty matrix
if(nrow(exp.degs) == 0){
  stop("No DEGs found — relax filtering thresholds")
}


column_ha = HeatmapAnnotation(sample.type = pheno$type)

Heatmap(exp100_DEGS,
        name = 'Exp',
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = column_ha)


############### 2D PCA #############

expression_t = t(exp100_DEGS)

# safety check (VERY IMPORTANT)
if(nrow(expression_t) < 2){
  stop("Not enough genes for PCA — DEG filtering too strict or empty set")
}

expression.pca = prcomp(expression_t, center = TRUE, scale. = TRUE)

summary(expression.pca)
autoplot(expression.pca, data = pheno, colour = 'type')


############### 3D PCA ##########
pca_scores = as.data.frame(expression.pca$x)

mycolors = ifelse(pheno$type == "Control", "blue", "red")

plot3d(pca_scores[,1:3],
       pch = 20,
       col = mycolors,
       type = 's',
       radius = 0.5)

############### Volcano plot ##########

EnhancedVolcano(res.df,
                lab = rownames(res.df),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.8)

# 1. boxplot
boxplot(exp100_DEGS,
        main = "processed genes Box Plot",
        col = seq_len(ncol(exp100_DEGS)))

# 2. histogram
hist(as.vector(exp100_DEGS),
     main = "processed genes Histogram")

# 3. # 3. Density plot side by side
par(mfrow = c(1, 2))

# --- Raw Data ---
exp_agg_numeric <- exp_agg[, sapply(exp_agg, is.numeric)]
raw_values <- as.vector(log2(as.matrix(exp_agg_numeric) + 1))

plot(density(raw_values),
     main = "Raw Data (log2)",
     xlab = "Expression",
     col = "red",
     lwd = 2)
polygon(density(raw_values),
        col = rgb(1, 0, 0, 0.3),
        border = "red")

# --- Processed DEGs ---
processed_values <- as.vector(exp100_DEGS)

plot(density(processed_values),
     main = "Processed DEGs",
     xlab = "Expression",
     col = "steelblue",
     lwd = 2)
polygon(density(processed_values),
        col = rgb(0.2, 0.5, 0.8, 0.3),
        border = "steelblue")

par(mfrow = c(1, 1))
