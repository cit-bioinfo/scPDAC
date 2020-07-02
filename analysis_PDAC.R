#Loading dataset 
load("~/Documents/projet_stage/data/data/pool_object.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/metadata_pool_object.RData")

#Data treatment 
#Add metadat
#metadata_pool_object <- read.table("metadata_pool_object.tsv", sep="\t", header=F, fill=T)
#metadata_pool_object <- metadata_pool_object[-1,]
#colnames(metadata_pool_object) <- c("cell", "orig_ident", "nCount_RNA", "nFeature_RNA", "percent_mt", "RNA_snn_res_1", "seurat_clusters", "peng_all_annotations")
#head(metadata_pool)

##Add ident
#active_ident_pool_object <- read.table("active_ident_pool_object.tsv", sep="\t", header=F, fill=T)
#active_ident_pool_object <- active_ident_pool_object[-1,]
#colnames(active_ident_pool_object) <- c("cell", "my_annotations")
#metadata_pool_object <- merge(active_ident_pool_object, metadata_pool_object, by="cell")
#metadata_pool_object$RNA_snn_res_1 <- NULL
#metadata_pool_object <- metadata_pool_object[, c(1, 4, 5, 6, 3, 7, 2, 8)]
#write.table(metadata_pool_object, file="~/Documents/projet_stage/data/analysis_pool_object/pool_metadata_object.tsv", sep="\t", row.names=F, col.names=T, quote=F)

##Add N/T orig ident
#active_ident_pool_object_orig_ident <- read.table("active_ident_pool_object_orig_ident.tsv", sep="\t", header=F, fill=T)
#active_ident_pool_object_orig_ident <- active_ident_pool_object_orig_ident[-1,]
#colnames(active_ident_pool_object_orig_ident) <- c("cell", "orig_ident_N/T")
#metadata_pool_object <- merge(metadata_pool_object, active_ident_pool_object_orig_ident, by.group="cell")
#metadata_pool_object <- metadata_pool_object[, c(1, 2, 3, 4, 5, 9, 6, 7, 8)]
#write.table(metadata_pool_object, file="~/Documents/projet_stage/data/analysis_pool_object/pool_metadata_object.tsv", sep="\t", row.names=F, col.names=T, quote=F)


#Analysis 
head(metadata_pool)
str(metadata_pool)         
summary(metadata_pool)
metadata_pool[, 4:6] <- as.numeric(unlist(metadata_pool[,4:6]))
metadata_pool$cell <- as.factor(metadata_pool$cell)
metadata_pool$cluster <- as.factor(metadata_pool$cluster)
metadata_pool$orig.ident <- as.factor(metadata_pool$orig.ident)
str(metadata_pool)
