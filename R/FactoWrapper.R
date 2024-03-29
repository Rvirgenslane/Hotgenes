#' FactoMiner wrapper used by Shiny_DE_Viz for PCA/HCPC of samples and genes
#' @import FactoMineR
#' @import factoextra
#' @param biplot Boolean value. If TRUE (Default), both individuals and variables
#' will be shown on the plot. If FALSE, only individuals will be plotted.
#' @param Top_var a numeric value indicating the number of top contributing
#' genes to show on the biplot.
#' @param ellipse.level numeric value for ellipse size, from 0 to 1.0
#' @param ellipse.alpha numeric value for ellipse transparency, from 0 to 1.0.
#' @param habillage_selection string indicating the grouping to be shown.
#' Default is "clust", which represents the HCPC ward clusters.
#' @param label_sel String indicating if individuals "ind" and variables
#' "var" should be labeled on the biplot.
#' @param readouts either "rlog" (default) or "vsd" indicating the variance
#' stabilization method preferred for visualization.
#' @param pointsize numeric value for point size.
#' @param labelsize numeric value for label size.
#' @param selected_contrast string indicating the DE contrast to be explored.
#' @param Output_DEseq2 R object generated by DEseq2_export function
#' @export
#' @return FactoMiner plot showing HCPC of genes, samples, and conditions
#' @examples
#' Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' FactoWrapper(Example_Hotgenes)

FactoWrapper<-function(Output_DEseq2=NULL,
biplot=TRUE,
readouts=1,
Top_var=10,
ellipse.level= 0.5,
ellipse.alpha=0.5,
habillage_selection="clust",
label_sel=c("ind", "var"),
pointsize = 1, labelsize = 1,
selected_contrast= 1){

title=""
# Normalized expression

gene_levels<-Output_DEseq2$Normalized_Expression[[readouts]]


# getting ids
contrast_name<-names(Output_DEseq2$Output_DE[selected_contrast])
contrast_name_lists <- setNames(vector(length(contrast_name),
mode="list"), contrast_name)

for (i in contrast_name) {
sig_ids<-rownames(Output_DEseq2$Output_DE[[i]])
contrast_name_lists[[i]]<-sig_ids
}

gene_ids<-unique(unlist(contrast_name_lists, use.names = FALSE))
tdm <- data.frame(t(gene_levels[gene_ids,]),
check.names = FALSE,
stringsAsFactors = FALSE)
design_input<-Output_DEseq2$design_data

# merging with phenodata
dm_pheno<-merge(tdm, design_input,
by = "row.names")
rownames(dm_pheno)<-dm_pheno$Row.names
dm_pheno$Row.names<-NULL
quant_sup<-colnames(design_input)

# PCA start
res<-PCA(dm_pheno, scale.unit=TRUE, ncp=5,
quali.sup=match(quant_sup,names(dm_pheno)),
graph = FALSE)

res.hcpc<-HCPC(res ,nb.clust=-1, consol=FALSE,min=3,max=5,graph=FALSE)
summary_PCA<-summary(res,
nb.dec = 3, nbelements=100, nbind = 100, ncp = 3, file="")

if(biplot==TRUE){
res_PPI_pa_1<-fviz_pca_biplot(res,
axes = c(1,2), repel = TRUE,  label=label_sel,
habillage=res.hcpc$data.clust[,habillage_selection] ,
col.quanti.sup = "red",
col.var = c("black"),
pointsize = pointsize, labelsize = labelsize,
select.var = list(contrib = Top_var),
col.ind.sup = "black", ellipse.alpha = ellipse.alpha,
title = title,
addEllipses=TRUE, ellipse.level=ellipse.level) + theme_classic() +
scale_shape_manual(values = rep(20,
length(res.hcpc$data.clust[,habillage_selection] )))
} else if(biplot==FALSE) {
res_PPI_pa_1<-fviz_pca_ind(res, axes = c(1,2), repel = TRUE,  label=label_sel,
habillage=res.hcpc$data.clust[,habillage_selection] ,
col.quanti.sup = "red",
col.var = c("black"),
pointsize = pointsize, labelsize = labelsize,
select.var = list(contrib = Top_var),
col.ind.sup = "black", ellipse.alpha = ellipse.alpha,
title = title,
addEllipses=TRUE, ellipse.level=ellipse.level) + theme_classic() +
scale_shape_manual(values = rep(20,
length(res.hcpc$data.clust[,habillage_selection] )))
}
output_PCA<-list(
res=res,
summary_PCA=summary_PCA,
res.hcpc=res.hcpc,
res_PPI_pa_1=res_PPI_pa_1)
return(output_PCA)
}
