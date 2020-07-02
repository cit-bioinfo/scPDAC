#Loading dataset 
load("~/Documents/projet_stage/data/data/pool_object.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/metadata_pool_object.RData")

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
#metadata_pool_object$my_annotations <- NULL
#metadata_pool_object <- metadata_pool_object[, c(1, 3, 4, 5, 6, 7, 8, 2, 9)]
#write.table(metadata_pool_object, file="~/Documents/projet_stage/data/analysis_pool_object/pool_metadata_object.tsv", sep="\t", row.names=F, col.names=T, quote=F)

##Add N/T orig ident
#active_ident_pool_object_orig_ident <- read.table("active_ident_pool_object_orig_ident.tsv", sep="\t", header=F, fill=T)
#active_ident_pool_object_orig_ident <- active_ident_pool_object_orig_ident[-1,]
#colnames(active_ident_pool_object_orig_ident) <- c("cell", "orig_ident_N/T")
#colnames(metadata_pool_object)[which(names(metadata_pool_object)=="orig_ident_N/T")] <- "orig_ident_N_T"
#metadata_pool_object <- merge(metadata_pool_object, active_ident_pool_object_orig_ident, by.group="cell")
#metadata_pool_object <- metadata_pool_object[, c(1, 2, 3, 4, 5, 9, 6, 7, 8)]
#write.table(metadata_pool_object, file="~/Documents/projet_stage/data/analysis_pool_object/pool_metadata_object.tsv", sep="\t", row.names=F, col.names=T, quote=F)

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
###All cells
myannotations_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast", "Stellate cell", "Macrophage", "T cell", "B cell")))
myannotations_table 
myannotations_barplot <- barplot(myannotations_table, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(myannotations_barplot, myannotations_table, paste(myannotations_table), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F, inset=c(0.02, 0.02))

#myannotations_barplot <- barplot(myannotations_table, col=c("lavender", "lightblue"), width=.3, beside=FALSE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)

### Ductal cells
myannotations_table_ductal_cells <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2")))
myannotations_table_ductal_cells 
myannotations_barplot_ductal_cells <- barplot(myannotations_table_ductal_cells, col=c("lavender", "lightblue"), beside=TRUE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(myannotations_barplot_ductal_cells, myannotations_table_ductal_cells, paste(myannotations_table_ductal_cells), cex=0.8, pos=3, col="gray45")
legend("top", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

###Prop_table
myannotations_proptable <- prop.table(myannotations_table)
myannotations_proptable
myannotations_proptable <- round(myannotations_proptable, 3)
myannotations_propbarplot <- barplot(myannotations_proptable, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 0.25), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8) 
text(myannotations_propbarplot,myannotations_proptable, paste(myannotations_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

##Pengall annotations
###All cells
pengall_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$peng_all_annotations, levels=c("Ductal cell type 1", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast cell", "Stellate cell", "Macrophage cell", "T cell", "B cell")))
pengall_table
pengall_barplot <- barplot(pengall_table, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(pengall_barplot,pengall_table, paste(pengall_table), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

###Ductal cells
pengall_table_ductal_cells <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$peng_all_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2"))) 
pengall_table_ductal_cells 
pengall_barplot_ductal_cells <- barplot(pengall_table_ductal_cells, col=c("lavender", "lightblue"), beside=TRUE, ylim=c(0, 12000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(pengall_barplot_ductal_cells, pengall_table_ductal_cells, paste(pengall_table_ductal_cells), cex=0.8, pos=3, col="gray45")
legend("top", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)


###Prop_table
pengall_proptable <- prop.table(pengall_table)
pengall_proptable
pengall_proptable <- round(pengall_proptable, 3)
pengall_propbarplot <- barplot(pengall_proptable, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 0.25), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
text(pengall_propbarplot,pengall_proptable, paste(pengall_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

