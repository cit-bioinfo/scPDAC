#Loading packages----------- 
library(plotrix)

library(devtools)
library(usethis)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

#Loading dataset----------- 
load("~/Documents/projet_stage/data/data/pool_object.RData")
load("~Documents/projet_stage/data/analysis_pool_subset/pool_subset.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/metadata_pool_object.RData")
load("~/Documents/projet_stage/data/analysis_pool_object/metadata_pool_subset.RData")

#POOL OBJECT - Analysis-----------

##My annotations-----------
###All cells
par(mfrow=c(1,1))
myannotations_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast", "Stellate cell", "Macrophage", "T cell", "B cell")))
myannotations_table 
myannotations_barplot <- barplot(myannotations_table, col=c("skyblue1", "red2"), width=.3, ylim=c(0, 12000),cex.names = 0.7, border=F,  font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45", beside=T)
text(myannotations_barplot, myannotations_table, paste(myannotations_table), cex=0.8, pos=3, col="gray45")
legend("topright", legend= c("Control samples", "Cancer samples"), col =c("skyblue1", "red2"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)
dev.off()

# myannotations_table
# myannotations_barplot_2 <- barplot(myannotations_table, col=c("lavender", "lightblue"), width=.3, beside=FALSE, ylim=c(0, 12000),cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45")
# legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("lavender", "lightblue"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

### Ductal cells
myannotations_table_ductal_cells <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2")))
myannotations_table_ductal_cells
myannotations_barplot_ductal_cells <- barplot(myannotations_table_ductal_cells, col=c("springgreen3", "red3"), beside=TRUE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45")
text(myannotations_barplot_ductal_cells, myannotations_table_ductal_cells, paste(myannotations_table_ductal_cells), cex=0.8, pos=3, col="gray45")
legend("top",  inset=c(-0.2, -0.2), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

###Prop_table
addmargins(myannotations_table)
myannotations_proptable <- prop.table(myannotations_table, 2)*100
myannotations_proptable <- round(myannotations_proptable, 3)
myannotations_proptable
myannotations_propbarplot <- barplot(myannotations_proptable, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 100), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45") 
text(myannotations_propbarplot,myannotations_proptable, paste(myannotations_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.01, -0.25), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9, horiz=T)

##Peng_al annotations-----------
###All cells
pengal_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$peng_al_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast cell", "Stellate cell", "Macrophage cell", "T cell", "B cell")))
pengal_table 
pengal_barplot <- barplot(pengal_table, col=c("skyblue1", "red2"), width=.3, beside=TRUE, ylim=c(0, 12000), cex.names = 0.7, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
text(pengal_barplot,pengal_table, paste(pengal_table), cex=0.8, pos=3, col="gray45")
par(xpd=TRUE)
legend("topright",legend= c("Control samples", "Cancer samples"), col =c("skyblue1", "red2"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.8)

# pengal_table
# pengal_barplot_2 <- barplot(pengal_table, col=c("springgreen3", "red3"), width=.3, beside=FALSE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45",  cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
# legend("topright", inset=c(-0.1, -0.2), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

###Ductal cells
pengal_table_ductal_cells <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$peng_al_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2"))) 
pengal_table_ductal_cells 
pengal_barplot_ductal_cells <- barplot(pengal_table_ductal_cells, col=c("springgreen3", "red3"), beside=TRUE, ylim=c(0, 12000), cex.names = 0.6, border=F, font.axis=2, col.axis="gray45", cex.axis=0.6, main="peng_al_annotations", col.main="gray45")
text(pengal_barplot_ductal_cells, pengal_table_ductal_cells, paste(pengal_table_ductal_cells), cex=0.8, pos=3, col="gray45")
legend("top",  inset=c(-0.2, -0.2), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

###Prop_table
addmargins(pengal_table)
pengal_proptable <- prop.table(pengal_table, 2)*100
pengal_proptable <- round(pengal_proptable, 3)
pengal_proptable
pengal_propbarplot <- barplot(pengal_proptable, col=c("springgreen3", "red3"), width=.3, las=2, beside=TRUE, ylim=c(0, 100), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
text(pengal_propbarplot,pengal_proptable, paste(pengal_proptable), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.01, -0.25), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", box.lty=0, bty="n", pch=20, pt.cex=2, cex=0.9, horiz=T)

##Clusters-----------
cluster_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$seurat_clusters, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))
cluster_table 
cluster_barplot <- barplot(cluster_table, col=c("skyblue1", "red2"), width=.3, beside=TRUE, ylim=c(0, 5000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_N/T", col.main="gray45")
text(cluster_barplot, cluster_table, paste(cluster_table), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.1,-0.2), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#prop table
cluster_table <- table(metadata_pool_object$orig_ident_N_T, factor(metadata_pool_object$seurat_clusters, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))
cluster_table <- prop.table(cluster_table, 2)*100
cluster_table <- round(cluster_table, 3)
cluster_barplot <- barplot(cluster_table, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 100), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_N/T", col.main="gray45")
text(cluster_barplot, cluster_table, paste(cluster_table), cex=0.8, pos=3, col="gray45")

# cluster_table
# cluster_barplot_2 <- barplot(cluster_table, col=c("springgreen3", "red3"), width=.3, beside=FALSE, ylim=c(0, 5000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8)
# legend("topright", inset=c(-0.1,-0.2), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)


#Complex HeatMap-----------
load("~/Documents/projet_stage/data/signatures/avg_module_score_object.RData")

heatmap <- Heatmap(matrix, name="mat", column_title="Clusters", column_title_side="top", column_title_gp = gpar(fontsize=15, fontface="bold"), column_dend_height=unit(2, "cm"), column_names_side="top",  column_names_rot=0, row_title ="Signatures", row_title_side="left", row_title_gp = gpar(fontsize=15, fontface="bold"), row_names_side="left", column_split=paste0(matrix2), row_order=rownames(matrix), top_annotation=ha)
heatmap
ha = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=1:10)))

#inferCNV-----------
par(mfrow=c(1,1))

#my_annotations / CNA
#all cells
myannotations_table_CNA <- table(metadata_pool_subset$CNA, factor(metadata_pool_subset$my_annotations, levels=c("Ductal cell 1", "Ductal cell 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast", "Stellate cell", "Macrophage", "T cell", "B cell")))
myannotations_table_CNA
myannotations_barplot_CNA <- barplot(myannotations_table_CNA, col=c("springgreen3", "red3"), beside=TRUE, ylim=c(0, 5000), las=2, cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations_CNA", col.main="gray45")
text(myannotations_barplot_CNA, myannotations_table_CNA, paste(myannotations_table_CNA), cex=0.8, pos=3, col="gray45")
legend("top",  inset=c(-0.2, -0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

#prop table
addmargins(myannotations_table_CNA)
myannotations_proptable_CNA <- prop.table(myannotations_table_CNA, 2)*100
myannotations_proptable_CNA <- round(myannotations_proptable_CNA, 3)
myannotations_proptable_CNA
myannotations_propbarplot_CNA <- barplot(myannotations_proptable_CNA, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 100), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations_CNA", col.main="gray45", horiz=F) 
text(myannotations_propbarplot_CNA,myannotations_proptable_CNA, paste(myannotations_proptable_CNA), cex=0.8, pos=3, col="gray45")

#pengal / CNA
#all cells
pengal_table_CNA <- table(metadata_pool_subset$CNA, factor(metadata_pool_subset$peng_al_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast cell", "Stellate cell", "Macrophage cell", "T cell", "B cell")))
pengal_table_CNA
pengal_barplot_CNA <- barplot(pengal_table_CNA, col=c("springgreen3", "red3"), beside=TRUE, ylim=c(0, 5000), cex.names = 0.8, las=2, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_CNA", col.main="gray45")
text(pengal_barplot_CNA, pengal_table_CNA, paste(pengal_table_CNA), cex=0.8, pos=3, col="gray45")
legend("top",  inset=c(-0.2, -0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=T)

#prop table
addmargins(pengal_table)
pengal_proptable_CNA <- prop.table(pengal_table_CNA, 2)*100
pengal_proptable_CNA <- round(pengal_proptable_CNA, 3)
pengal_proptable_CNA
pengal_propbarplot_CNA <- barplot(pengal_proptable_CNA, col=c("springgreen3", "red3"), width=.3, las=2, beside=TRUE, ylim=c(0, 100), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations_CNA", col.main="gray45")
text(pengal_propbarplot_CNA,pengal_proptable_CNA, paste(pengal_proptable_CNA), cex=0.8, pos=3, col="gray45")

#clusters 
cluster_table_CNA <- table(metadata_pool_subset$CNA, factor(metadata_pool_subset$seurat_clusters, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))
cluster_table_CNA
cluster_barplot_CNA <- barplot(cluster_table_CNA, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 1200), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_CNA", col.main="gray45")
text(cluster_barplot_CNA, cluster_table_CNA, paste(cluster_table_CNA), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.02,-0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

cluster_table_CNA <- table(metadata_pool_subset$CNA, factor(metadata_pool_subset$seurat_clusters, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34")))
cluster_table_CNA <- prop.table(cluster_table_CNA, 2)*100
cluster_table_CNA <- round(cluster_table_CNA, 3)
cluster_barplot_CNA <- barplot(cluster_table_CNA, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 100), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_CNA", col.main="gray45")
text(cluster_barplot_CNA, cluster_table_CNA, paste(cluster_table_CNA), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.02,-0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

cluster_table_CNA_positiv <- cluster_table_CNA[, -15]
cluster_barplot_CNA_positiv <- barplot(cluster_table_CNA_positiv, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 1500), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_CNA_positiv", col.main="gray45")
text(cluster_barplot_CNA_positiv, cluster_table_CNA_positiv, paste(cluster_table_CNA_positiv), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.02,-0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

cluster_table_CNA_positiv <- prop.table(cluster_table_CNA_positiv, 2)*100
cluster_table_CNA_positiv <- round(cluster_table_CNA_positiv, 3)
cluster_barplot_CNA_positiv <- barplot(cluster_table_CNA_positiv, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 120), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_CNA", col.main="gray45")
text(cluster_barplot_CNA_positiv, cluster_table_CNA_positiv, paste(cluster_table_CNA_positiv), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.02,-0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#ductal cell 2 + CNA
#my_annotations
cluster_ductal_cell_2_CNA_my_annotations <- metadata_pool_subset
cluster_ductal_cell_2_CNA_my_annotations[,3] <- NULL
cluster_ductal_cell_2_CNA_my_annotations <- cluster_ductal_cell_2_CNA[cluster_ductal_cell_2_CNA$my_annotations=="Ductal cell 2",] 

cluster_ductal_cell_2_CNA_table_my_annotations <- table(cluster_ductal_cell_2_CNA_my_annotations$CNA, factor(cluster_ductal_cell_2_CNA_my_annotations$seurat_clusters, levels=c("2", "9", "14", "18", "20", "24", "25")))
cluster_ductal_cell_2_CNA_table_my_annotations
cluster_ductal_cell_2_CNA_barplot_my_annotations <- barplot(cluster_ductal_cell_2_CNA_table_my_annotations, col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 1500), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_ductal_cell_2_CNA_my_annotations", col.main="gray45")
text(cluster_ductal_cell_2_CNA_barplot_my_annotations, cluster_ductal_cell_2_CNA_table_my_annotations, paste(cluster_ductal_cell_2_CNA_table_my_annotations), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.005,-0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#peng_al_annotations 
cluster_ductal_cell_2_CNA_peng_al <- cluster_ductal_cell_2_CNA_my_annotations
cluster_ductal_cell_2_CNA_peng_al[,2] <- NULL

cluster_ductal_cell_2_CNA_table_peng_al  <- table(cluster_ductal_cell_2_CNA_peng_al $CNA, factor(cluster_ductal_cell_2_CNA_peng_al $seurat_clusters, levels=c("2","3", "8", "9", "14", "18", "19", "20", "22", "23", "24", "25")))
cluster_ductal_cell_2_CNA_table_peng_al 
cluster_ductal_cell_2_CNA_barplot_peng_al  <- barplot(cluster_ductal_cell_2_CNA_table_peng_al , col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 1500), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="clusters_ductal_cell_2_CNA_peng_al", col.main="gray45")
text(cluster_ductal_cell_2_CNA_barplot_peng_al , cluster_ductal_cell_2_CNA_table_peng_al , paste(cluster_ductal_cell_2_CNA_table_peng_al ), cex=0.8, pos=3, col="gray45")
legend("topright", inset=c(-0.005,-0.2), legend= c("CNA-", "CNA+"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#pourcentage de cellules tumorales ne faisant pas parties des cellules ductal 2 
p <- prop.table(cluster_table_CNA_positiv)*100
p <- p[, c(1, 4, 5, 6, 9, 11, 12, 15, 19, 23)]
p <- prop.table(p)*100
rowSums(p)
#CNA+ : 1.0508
#CNA- : 98.9492

#POOL SUBSET - Analysis-----------
##Clusters / CNA-----------
clusters_table_CNA <- table(metadata_pool_subset_analysis$CNA, factor(metadata_pool_subset_analysis$RNA_snn_res_0.5, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21")))
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

##Peng_al annotations-----------
#annotations_table <- table(metadata_pool_subset_analysis$orig_ident_N_T, factor(metadata_pool_subset_analysis$peng_all_annotations, levels=c("Ductal cell type 1", "Ductal cell type 2", "Acinar cell", "Endocrine cell", "Endothelial cell", "Fibroblast cell", "Stellate cell", "Macrophage cell", "T cell", "B cell")))
#annotations_table 
#annotations_barplot <- barplot(annotations_table,col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
#text(annotations_barplot,annotations_table, paste(annotations_table), cex=0.8, pos=3, col="gray45")
#legend("topright", inset=c(-0.25,-0.05), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#addmargins(annotations_table)
#annotations_proptable <- prop.table(annotations_table)*100
#annotations_proptable <- round(annotations_proptable, 3)
#annotations_barplot_proptable <- barplot(annotations_proptable,col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 120), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
#text(annotations_barplot_proptable,annotations_proptable, paste(annotations_proptable), cex=0.8, pos=3, col="gray45")
#legend("topright", inset=c(-0.25,-0.05), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

##My annotations-----------
#annotations_table 
#annotations_barplot <- barplot(annotations_table,col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 12000), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="my_annotations", col.main="gray45")
#text(annotations_barplot,annotations_table, paste(annotations_table), cex=0.8, pos=3, col="gray45")
#legend("topright", inset=c(-0.25,-0.05), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

#addmargins(annotations_table)
#annotations_proptable <- prop.table(annotations_table)*100
#annotations_proptable <- round(annotations_proptable, 3)
#annotations_proptable
#annotations_barplot_proptable <- barplot(annotations_proptable,col=c("springgreen3", "red3"), width=.3, beside=TRUE, ylim=c(0, 120), cex.names = 0.8, border=F, font.axis=2, col.axis="gray45", cex.axis=0.8, main="peng_al_annotations", col.main="gray45")
#text(annotations_barplot_proptable,annotations_proptable, paste(annotations_proptable), cex=0.8, pos=3, col="gray45")
#legend("topright", inset=c(-0.25,-0.05), legend= c("Control samples", "Cancer samples"), col =c("springgreen3", "red3"), text.col="gray45", bty="n", pch=20, pt.cex=2, cex=0.8, horiz=F)

##Complex HeatMap-----------

load("~/Documents/projet_stage/data/signatures/avg_module_score_subset.RData")

col = colorRamp2(seq(min(matrix_subset), max(matrix_subset), length=3), c("blue", "white", "red"))
top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =2:6)))
heatmap_subset <- Heatmap(matrix_subset, name="mat", col=col, column_title="Clusters", column_title_side="top", column_title_gp = gpar(fontsize=15, fontface="bold"), column_dend_height=unit(2, "cm"), column_names_side="top",  column_names_rot=0, row_title ="Signatures", row_title_side="left", row_title_gp = gpar(fontsize=15, fontface="bold"), row_names_side="left", show_column_dend = TRUE, show_row_dend = FALSE, row_order=rownames(matrix_2_subset), row_split=paste0(matrix_2_subset))
heatmap_subset
