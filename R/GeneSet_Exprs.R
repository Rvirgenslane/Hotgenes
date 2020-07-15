#' Exports Expression data for a GeneSet of interest
#' @param GSEA_Hotgenes HotgenesObj appended by BatchGSEA function
#' @param contrastID contrast selection
#' @param pathwayID String or numeric index indicating the desired pathway
#' @param normalizedData indicates the normalized data in
#' HotgenesObj to return.
#' @return Expresion data
#' @export
#' @examples Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' library(msigdbr)
#' m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#' qbat<-BatchGSEA(m_df= m_df, HotgenesObj=Example_Hotgenes)
#' expsD<-GeneSet_Exprs(GSEA_Hotgenes=qbat)
#' library(pheatmap)
#' pheatmap(expsD$E,scale = "row",
#' main=expsD$Main)


# pheatmap_exporting ------------------------------------------------------

GeneSet_Exprs<-function(GSEA_Hotgenes=NULL,
normalizedData=1,
contrastID=1,
pathwayID=1){

if(is.null(GSEA_Hotgenes$OuputGSEA)) {
stop("No GSEA detected, make sure to run BatchGSEA first")
}

# Getting data    
data_heat<-GSEA_Hotgenes$Normalized_Expression[[normalizedData]]

# Getting GSEA
ht_batch<-GSEA_Hotgenes$OuputGSEA

# Getting contrast
c_id<-names(ht_batch[contrastID])

# Getting list of pathways
sel_pathway<-ht_batch[[contrastID]]$top

# Selecting pathway by index number or string
if(is.numeric(pathwayID)){
pID<-sel_pathway[pathwayID,]$pathway
ids_leadingEdge<-unlist(sel_pathway[pathwayID,]$leadingEdge)
}else if (is.character(pathwayID)){

pID<-sel_pathway[sel_pathway$pathway == pathwayID,]$pathway
ids_leadingEdge<-unlist(sel_pathway[sel_pathway$pathway == pathwayID,]$leadingEdge)
}  

# Subsetting expression data based on leading Edge genes
E<-data_heat[rownames(data_heat) %in% ids_leadingEdge,]

# making pheatmap main text label
Main<-paste(c_id, pID, sep = ": ")

# Final lists
ExpOut<-list(E=E,Main=Main)

return(ExpOut)}
