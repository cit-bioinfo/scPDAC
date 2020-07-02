#Loading dataset 
load("~/Documents/projet_stage/data/data/pool_object.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/metadata_pool_object.RData")
library(ggplot2)

#Data treatment 
#Add metadata
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
#colnames(metadata_pool_object)[which(names(metadata_pool_object)=="orig_ident_N/T")] <- "orig_ident_N_T"
#metadata_pool_object <- merge(metadata_pool_object, active_ident_pool_object_orig_ident, by.group="cell")
#metadata_pool_object <- metadata_pool_object[, c(1, 2, 3, 4, 5, 9, 6, 7, 8)]
#write.table(metadata_pool_object, file="~/Documents/projet_stage/data/analysis_pool_object/pool_metadata_object.tsv", sep="\t", row.names=F, col.names=T, quote=F)

##
# head(metadata_pool_object)
# str(metadata_pool_object) 
# summary(metadata_pool)
# metadata_pool_object[, 2:4] <- as.numeric(unlist(metadata_pool_object[,2:4]))
# metadata_pool_object$cell <- as.factor(metadata_pool_object$cell)
# metadata_pool_object$orig_ident <- as.factor(metadata_pool_object$orig_ident)
# metadata_pool_object$orig_ident_N_T <- as.factor(metadata_pool_object$orig_ident_N_T)
# metadata_pool_object$seurat_clusters <- as.factor(metadata_pool_object$seurat_clusters)
# metadata_pool_object$my_annotations <- as.factor(metadata_pool_object$my_annotations)
# metadata_pool_object$peng_all_annotations <- as.factor(metadata_pool_object$peng_all_annotations)
# str(metadata_pool_object)

#Analysis
##Myannotations
myannotations_table <- table(metadata_pool_object$orig_ident_N_T, metadata_pool_object$my_annotations)
myannotations_barplot <- barplot(myannotations_table, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, levels=c(3, 4, 1, 5, 6, 7, 9, 2, 10))
text(myannotations_barplot, myannotations_table, paste(myannotations_table), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

myannotations_proptable <- prop.table(myannotations_table)
myannotations_proptable
myannotations_proptable <- round(myannotations_proptable, 3)
myannotations_propbarplot <- barplot(myannotations_proptable, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 0.25), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8) 
text(myannotations_propbarplot,myannotations_proptable, paste(myannotations_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

#Pengall annotations
pengall_table <- table(metadata_pool_object$orig_ident_N_T, metadata_pool_object$peng_all_annotations)
pengall_barplot <- barplot(pengall_table, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(pengall_barplot,pengall_table, paste(pengall_table), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

pengall_proptable <- prop.table(pengall_table)
pengall_proptable
pengall_proptable <- round(pengall_proptable, 3)
pengall_propbarplot <- barplot(pengall_proptable, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 0.25), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(pengall_propbarplot,pengall_proptable, paste(pengall_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

