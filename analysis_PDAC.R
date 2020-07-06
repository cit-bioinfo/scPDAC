#Loading dataset----------- 
load("~/Documents/projet_stage/data/data/pool_object.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/metadata_pool_object.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/barplot_data.RData")

library(plotrix)
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

#names(metadata_pool_object)[match("peng_all_annotations", names(metadata_pool_object))] <- "peng_al_annotations"
#Analysis-----------

##Myannotations-----------
###All cells
par(mfrow=c(1,2))
myannotations_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast", "Stellate cell", "Macrophage", "T cell", "B cell")))
myannotations_table 
myannotations_barplot <- barplot(myannotations_table, col=c("lavender", "lightblue"), width=.3, ylim=c(0, 12000),cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45", beside=T)
text(myannotations_barplot, myannotations_table, paste(myannotations_table), cex=0.8, pos=3, col="gray45")
par(xpd=TRUE)
legend("toprigh", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

myannotations_table
myannotations_barplot_2 <- barplot(myannotations_table, col=c("lavender", "lightblue"), width=.3, beside=FALSE, ylim=c(0, 12000),cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45")
legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)


### Ductal cells
myannotations_table_ductal_cells <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2")))
myannotations_table_ductal_cells 
myannotations_barplot_ductal_cells <- barplot(myannotations_table_ductal_cells, col=c("lavender", "lightblue"), beside=TRUE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45")
text(myannotations_barplot_ductal_cells, myannotations_table_ductal_cells, paste(myannotations_table_ductal_cells), cex=0.8, pos=3, col="gray45")
legend("top",  inset=c(-0.2, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

###Prop_table
myannotations_proptable <- prop.table(myannotations_table)
myannotations_proptable
myannotations_proptable <- round(myannotations_proptable, 3)
myannotations_propbarplot <- barplot(myannotations_proptable, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 0.25), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45") 
text(myannotations_propbarplot,myannotations_proptable, paste(myannotations_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

##Peng_al annotations-----------
###All cells
pengal_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$peng_al_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast cell", "Stellate cell", "Macrophage cell", "T cell", "B cell")))
pengal_table
pengal_barplot <- barplot(pengal_table, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
text(pengal_barplot,pengal_table, paste(pengal_table), cex=0.8, pos=3, col="gray45")
par(xpd=TRUE)
legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.8)

pengal_table
pengal_barplot_2 <- barplot(pengal_table, col=c("lavender", "lightblue"), width=.3, beside=FALSE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45",  cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

###Ductal cells
pengal_table_ductal_cells <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$peng_al_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2"))) 
pengal_table_ductal_cells 
pengal_barplot_ductal_cells <- barplot(pengal_table_ductal_cells, col=c("lavender", "lightblue"), beside=TRUE, ylim=c(0, 12000), cex.names = 0.6, border=F, font.axis=2, col.axis="gray45", cex.axis=0.6, main="peng_al_annotations", col.main="gray45")
text(pengal_barplot_ductal_cells, pengal_table_ductal_cells, paste(pengal_table_ductal_cells), cex=0.8, pos=3, col="gray45")
legend("top",  inset=c(-0.2, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

###Prop_table
pengal_proptable <- prop.table(pengal_table)
pengal_proptable
pengal_proptable <- round(pengal_proptable, 3)
pengal_propbarplot <- barplot(pengal_proptable, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 0.25), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
text(pengal_propbarplot,pengal_proptable, paste(pengal_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9)

##Clusters-----------
cluster_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$seurat_clusters, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))
cluster_table
cluster_barplot <- barplot(cluster_table, col=c("lavender", "lightblue"), width=.3, beside=TRUE, ylim=c(0, 5000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_N/T", col.main="gray45")
text(cluster_barplot, cluster_table, paste(cluster_table), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.1,-0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

# cluster_table
# cluster_barplot_2 <- barplot(cluster_table, col=c("lavender", "lightblue"), width=.3, beside=FALSE, ylim=c(0, 5000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
# legend("topright", inset=c(-0.1,-0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

# test_table <- table(metadata_pool_object$seurat_clusters, metadata_pool_object$my_annotations)
# test_table
# test_barplot <- barplot(test_table, width=.3, ylim=c(0, 15000), col=rainbow(34), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_celltype", col.main="gray45")
# legend("topright", inset=c(0.05,-0.2), legend= c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34"), col =rainbow(34), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.35, horiz=T)

# test2_table <- table(metadata_pool_object$seurat_clusters, metadata_pool_object$peng_all_annotations)
# test2_table
# barplot(test2_table)

# test3_table <- table(metadata_pool_object$my_annotations, metadata_pool_object$peng_all_annotations)
# test3_table
# barplot(test3_table)
