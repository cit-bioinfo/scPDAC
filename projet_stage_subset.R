#Loading packages----
library(Seurat)
library(ggraph)
library(ggplot2)
library(clustree)

#Subset of epithelial cells ----
#It concerns cells of clusters highly rich in tumor cells : n°2, 9, 14, 18, 20, 24, 25. 

#Add metadata to Seurat object : CNA column
metadata_pool_subset <- metadata_pool_subset[, c(1, 10, 11)]
metadata_pool <- merge(x=metadata_pool_object, y=metadata_pool_subset, by="cell", all.x=TRUE)

infer_CNV_N_colMeans$cell <- rownames(infer_CNV_N_colMeans)
metadata_pool <- merge(x=metadata_pool, y=infer_CNV_N_colMeans, by="cell", all.x=TRUE)
metadata_pool$infer_CNV_N <- metadata_pool$infer_CNV_N_colMeans
metadata_pool$infer_CNV_N_colMeans <- NULL
metadata_pool <- metadata_pool[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 11)]

metadata_pool$infer_CNV <- metadata_pool$infer_CNV_T
metadata_pool$infer_CNV[!is.na(metadata_pool$infer_CNV_N)] = metadata_pool$infer_CNV_N[!is.na(metadata_pool$infer_CNV_N)]
metadata_pool[, c(10, 11)] <- NULL
metadata_pool <- metadata_pool[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10)]

summary(metadata_pool)
str(metadata_pool)

load("/home/rose/Documents//projet_stage/data/data/pool_object.RData")
load("/home/rose/Documents//projet_stage/data/analysis_pool_object/metadata_pool.RData")

rownames(metadata_pool) <- metadata_pool$cell
metadata_pool$cell <- NULL
metadata_pool_CNA <- metadata_pool[, c(9, 10)]
pool <- AddMetaData(pool, metadata=metadata_pool_CNA)
pool

#Creation of a new seurat object : subset of cells from specific clusters
Idents(pool) <- "seurat_clusters"
pool_subset <- subset(x=pool, idents=c("2", "9", "14", "18", "20", "24", "25"))
pool_subset
#11467 cells 

#Seurat Pipeline ----

## **Standard pre-processing workflow - Data filtering**
### Calculate the percentage of mitochondrial genes 
pool_subset[["percent.mt"]] <- PercentageFeatureSet(object = pool_subset, pattern = "^MT-")
head(pool_subset@meta.data, 5)

## **Normalizing data - logNormalize**
### Normalizing the data
pool_subset <- NormalizeData(pool_subset)
#pool_subset[["RNA"]]@data

### Identification of highly variable features
pool_subset <- FindVariableFeatures(object=pool_subset, selection.method = "vst", nfeatures = 2000)
top10_variablefeatures_subset <- head(VariableFeatures(pool_subset), 10)
top10_variablefeatures_subset

## **Scaling data**
pool_subset <- ScaleData(pool_subset)

## **PCA**
pool_subset <- RunPCA(pool_subset, features = VariableFeatures(object = pool_subset))
DimHeatmap(pool_subset, dims=10:20, cells=500, balanced=TRUE)

pool_subset <- JackStraw(pool_subset, num.replicate = 100)
pool_subset <- ScoreJackStraw(pool_subset, dims=1:20)
JackStrawPlot(pool_subset, dims=1:20)
#0 : p-value is smaller than can be represented in R
ElbowPlot(pool_subset, ndims = 50, reduction = "pca")
#gap between pca 17 and 18, dimension = 1:17

## **Cluster the cell**
pool_subset <- FindNeighbors(pool_subset, reduction='pca', dims = 1:17)
pool_subset <- FindClusters(pool_subset, resolution = 0.2) 
head(pool_subset@meta.data, 5)
#0.1 : 13 communities
#0.2 : 16 communities

clustree(pool_subset, prefixe="RNA_snn_res.0.1")

#pool_subset[["RNA_snn_res.1.5"]] <- NULL 

## **Run non-linear dimensional reduction**
### Uniform Manifold Approximation and Projection (UMAP) 
pool_subset <- RunUMAP(pool_subset, dims = 1:17)
plot2 <- DimPlot(pool_subset, reduction="umap", pt.size=0.6, label=T)
plot2
plot3 <- DimPlot(pool_subset, reduction="umap", pt.size=0.6, group.by="CNA", cols=c("springgreen3", "red3", "black"))
plot3

#Barplots----
##Data treatment
write.table(pool_subset@meta.data, file="~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset.tsv", sep="\t", row.names=T, col.names=T, quote=F)
metadata_pool_subset_analysis <- read.table("~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset.tsv", sep="\t", header=F, fill=T)
metadata_pool_subset_analysis <- metadata_pool_subset[-1,]
colnames(metadata_pool_subset_analysis) <- c("cell", "orig_ident", "nCount_RNA", "nFeature_RNA", "peng_all_annotations", "infer_CNV", "CNA", "percent_mt", "seurat_clusters", "RNA_snn_res.0.1", "RNA_snn_res.0.2")
metadata_pool_subset_analysis <- metadata_pool_subset_analysis[, c(1, 2, 3, 4, 5, 12, 6, 7, 8, 9, 10, 11)]
metadata_pool_subset_analysis <- merge(metadata_pool_subset_analysis, metadata_pool_object, by="cell")
metadata_pool_subset_analysis[, c(15, 14, 13, 12)] <- NULL
colnames(metadata_pool_subset_analysis) <- c("cell", "nCount_RNA", "nFeature_RNA",  "percent_mt", "orig_ident", "orig_ident_N_T", "seurat_clusters",  "RNA_snn_res.0.1", "RNA_snn_res.0.2", "peng_all_annotations", "infer_CNV", "CNA")
head(metadata_pool_subset_analysis)

write.table(metadata_pool_subset_analysis, file="~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset.tsv", quote=FALSE, sep="\t", row.names=F)
load("~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset_analysis.RData")

##Barplots

