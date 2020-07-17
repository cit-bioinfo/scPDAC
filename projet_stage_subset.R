#Subset of epithelial cells.
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





pool <- AddMetaData(pool, metadata = metadata_pool, col.name="" )
