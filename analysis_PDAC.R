#Loading dataset 

# metadata_pool <- read.table("metadata_pool.tsv", sep="\t", header=F, fill=T)
# metadata_pool <- metadata_pool[-1,]
# colnames(metadata_pool) <- c("cell", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res_1", "seurat_clusters")
# head(metadata_pool)
# 
# metadata_pool2 <- read.table("metadata_pool2.tsv", sep="\t", header=F, fill=T)
# metadata_pool2 <- metadata_pool2[-1,]
# colnames(metadata_pool2) <- c("cell", "orig.ident", "nCount_RNA", "nFeature_RNA", "seurat_clusters")
# head(metadata_pool2)
# write.table(metadata_pool2, file="~/Documents/projet_stage/data/analysis/pool2_metadata.tsv", sep="\t", row.names=T, col.names=T, quote=F)
# 
# 
# active_ident_pool <- read.table("active_ident_pool.tsv", sep="\t", header=F, fill=T)
# active_ident_pool <- active_ident_pool[-1,]
# colnames(active_ident_pool) <- c("cell", "cluster")
# metadata_pool <- merge(active_ident_pool, metadata_pool, by="cell")
# write.table(metadata_pool, file="~/Documents/projet_stage/data/analysis/pool_metadata.tsv", sep="\t", row.names=T, col.names=T, quote=F)

load("~/Documents/projet_stage/data/analysis/analysis.RData")


         