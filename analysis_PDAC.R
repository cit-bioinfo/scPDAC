#Loading dataset 

#pool
metadata_pool <- read.table("metadata_pool.tsv", sep="\t", header=F, fill=T)
metadata_pool <- metadata_pool[-1,]
colnames(metadata_pool) <- c("cell", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res_1", "seurat_clusters")
head(metadata_pool)

active_ident_pool <- read.table("active_ident_pool.tsv", sep="\t", header=F, fill=T)
active_ident_pool <- active_ident_pool[-1,]
colnames(active_ident_pool) <- c("cell", "cluster")
metadata_pool <- merge(active_ident_pool, metadata_pool, by="cell")
write.table(metadata_pool, file="~/Documents/projet_stage/data/analysis/pool_metadata.tsv", sep="\t", row.names=T, col.names=T, quote=F)

pool <- AddMetaData(pool, metadata = metadata, col.name="peng_all_annotations")
pool[['RNA_snn_res.1.5']] <- NULL
pool[['RNA_snn_res.0.5']] <- NULL

#pool2
metadata_pool2 <- read.table("metadata_pool2.tsv", sep="\t", header=F, fill=T)
metadata_pool2 <- metadata_pool2[-1,]
colnames(metadata_pool2) <- c("cell", "orig.ident", "nCount_RNA", "nFeature_RNA", "seurat_clusters")
head(metadata_pool2)
write.table(metadata_pool2, file="~/Documents/projet_stage/data/analysis/pool2_metadata.tsv", sep="\t", row.names=T, col.names=T, quote=F)

load("~/Documents/projet_stage/data/analysis/analysis.RData")
load("/home/rose/Documents//projet_stage/data/data_clusters.RData")
library(ggplot2)

#Analysis pool
head(metadata_pool)
str(metadata_pool)         
summary(metadata_pool)
metadata_pool[, 4:6] <- as.numeric(unlist(metadata_pool[,4:6]))
metadata_pool$cell <- as.factor(metadata_pool$cell)
metadata_pool$cluster <- as.factor(metadata_pool$cluster)
metadata_pool$orig.ident <- as.factor(metadata_pool$orig.ident)
str(metadata_pool)

#Analysis pool2
head(metadata_pool2)
str(metadata_pool2)         
summary(metadata_pool2)
metadata_pool2[, 3:4] <- as.numeric(unlist(metadata_pool2[,3:4]))
metadata_pool2$cell <- as.factor(metadata_pool2$cell)
metadata_pool2$seurat_clusters <- as.factor(metadata_pool2$seurat_clusters)
metadata_pool2$orig.ident <- as.factor(metadata_pool2$orig.ident)
str(metadata_pool2)   

hist(metadata_pool2[,4], xlab="number of counts")


