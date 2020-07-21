#Loading dataset and packages----
library(Seurat)
library(ggraph)
library(ggplot2)
library(clustree)

load("/home/rose/Documents//projet_stage/data/analysis_pool_subset/pool_subset.RData")
load("/home/rose/Documents//projet_stage/data/analysis_pool_object/metadata_pool_subset.RData") #object

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
pool_subset <- FindClusters(pool_subset, resolution = 0.5)
head(pool_subset@meta.data, 5)

clustree(pool_subset, prefixe="RNA_snn_res.0.1")

#pool_subset[["RNA_snn_res.1.5"]] <- NULL 

## **Run non-linear dimensional reduction**
### Uniform Manifold Approximation and Projection (UMAP) 
pool_subset <- RunUMAP(pool_subset, dims = 1:17)
plot2 <- DimPlot(pool_subset, reduction="umap", pt.size=0.6, label=T)
plot2+ labs(title="UMAP res 1 / 31 communities")
plot3 <- DimPlot(pool_subset, reduction="umap", pt.size=0.6, group.by="CNA", cols=c("springgreen3", "red3", "black"))
plot3+ labs(title="UMAP CNA")


Idents(pool_subset) <- "orig.ident"
experiment.merged <- RenameIdents(object = pool_subset, 'N11'='N1', 'N10'='N1', 'N9'='N1', 'N8'='N1', 'N7'='N1', 'N6'='N1', 'N5'='N1', 'N4'='N1', 'N3'='N1', 'N2'='N1', 'T2'='T1', 'T3'='T1','T4'='T1','T5'='T1','T6'='T1','T7'='T1','T8'='T1','T9'='T1','T10'='T1','T11'='T1','T12'='T1','T13'='T1','T14'='T1','T15'='T1','T16'='T1','T17'='T1','T18'='T1','T19'='T1','T20'='T1','T21'='T1','T22'='T1','T23'='T1','T24'='T1')
plot6 <-DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "umap")
new.orig_ident.ids <- c("Control samples", "Cancer samples")
names(new.orig_ident.ids) <- levels(experiment.merged)
pool_subset <- RenameIdents(experiment.merged, new.orig_ident.ids)

summary(Idents(pool_subset))
#43 control samples
#11424 cancer samples
plot4 <- DimPlot(pool_subset, reduction="umap", pt.size=0.6, cols=c("springgreen3", "red3"))
plot4+ labs(title="UMAP N/T")


#Barplots----
##Data treatment
metadata_pool_subset_analysis_4 <- read.table("~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset_4.tsv", sep="\t", header=F, fill=T)
write.table(pool_subset@meta.data, file="~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset_4.tsv", sep="\t", row.names=T, col.names=T, quote=F)
metadata_pool_subset_analysis_4 <- metadata_pool_subset_analysis_4[-1,]
colnames(metadata_pool_subset_analysis_4) <- c("cell","RNA_snn_res.1", "seurat_clusters", "RNA_snn_res.0.5", "RNA_snn_res.0.1", "RNA_snn_res.0.2", "RNA_snn_res.0.3", "RNA_snn_res.0.4", "RNA_snn_res.0.6", "RNA_snn_res.0.7", "RNA_snn_res.0.8", "RNA_snn_res.0.9")
metadata_pool_subset_analysis <- metadata_pool_subset_analysis[, c(1,13, 14, 15, 16, 17, 18, 5, 6, 7, 8, 4, 9, 10, 11, 12, 2, 3, 19, 20)]
metadata_pool_subset_analysis <- merge(metadata_pool_subset_analysis_4, metadata_pool_subset_analysis, by="cell")
metadata_pool_subset_analysis <- metadata_pool_subset_analysis[, -(c(18, 19, 20))]
colnames(metadata_pool_subset_analysis) <- c("cell", "nCount_RNA", "nFeature_RNA",  "percent_mt", "orig_ident", "orig_ident_N_T", "seurat_clusters",  "RNA_snn_res.0.1", "RNA_snn_res.0.2", "peng_all_annotations", "infer_CNV", "CNA")
head(metadata_pool_subset_analysis)

metadata_pool_subset_analysis$CNA[is.na(metadata_pool_subset_analysis$CNA)]<- "NA"
metadata_pool_subset_analysis$infer_CNV[is.na(metadata_pool_subset_analysis$infer_CNV)]<- "NA"

write.table(metadata_pool_subset_analysis, file="~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset.tsv", quote=FALSE, sep="\t", row.names=F)
load("~/Documents/projet_stage/data/analysis_pool_subset/metadata_pool_subset.RData") #subset

##Barplots
head(pool_subset@meta.data, 5)
summary(metadata_pool_subset_analysis)

#clusters / CNA
clusters_table_CNA <- table(metadata_pool_subset_analysis$CNA, factor(metadata_pool_subset_analysis$RNA_snn_res.0.4, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")))
clusters_table_CNA 
clusters_barplot_CNA <- barplot(clusters_table_CNA, col=c("springgreen3", "red3", "black"), width=.3, beside=TRUE, ylim=c(0, 1000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_CNA", col.main="gray45")
text(clusters_barplot_CNA, clusters_table_CNA, paste(clusters_table_CNA), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.02,-0.05), legend= c("CNA-", "CNA+", "NA"), col =c("springgreen3", "red3", "black"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#prop_table 
addmargins(clusters_table_CNA)
clusters_proptable_CNA <- prop.table(clusters_table_CNA, 2)*100
clusters_proptable_CNA <- round(clusters_proptable_CNA, 3)
clusters_proptable_CNA
clusters_propbarplot_CNA <- barplot(clusters_proptable_CNA, col=c("springgreen3", "red3", "black"), width=.3, beside=TRUE, ylim=c(0, 100), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_CNA_res_0.5", col.main="gray45")
text(clusters_propbarplot_CNA,clusters_proptable_CNA, paste(clusters_proptable_CNA), cex=0.8, pos=3, col="gray45")



