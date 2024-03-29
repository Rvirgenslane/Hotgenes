#' Generates a custom GeneSet from a Hotgenes object
#' @export
#' @param Hotgenes_obj R object generated by DEseq2_export or 
#' Limma_export functions
#' @return data.frame with "gene_symbol" and "gs_name".
#' GeneSet contains genes that are either Enriched or Depleted by the
#' condition (based on logfoldchange).
#' @examples Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' Custom_GeneSet<-Hotgenes2GeneSet(Example_Hotgenes)
#' head(Custom_GeneSet)

Hotgenes2GeneSet<-function(Hotgenes_obj=NULL){

# Getting ids
DE_ids<-names(Hotgenes_obj$Output_DE)
Enriched_map<-setNames(vector(length(DE_ids), mode = "list"),DE_ids)
Depleted_map<-setNames(vector(length(DE_ids), mode = "list"),DE_ids)

# Building lists
for (i in DE_ids) {
Tmp_df<-Hotgenes_obj$Output_DE[[i]]
Tmp_df$gene_symbol<-rownames(Tmp_df)
Tmp_df$gs_name<-c("")

# Checking for Enrichment/Depletion
Enriched_length<-length(Tmp_df[Tmp_df$log2FoldChange > 0,]$gs_name)
Depleted_length<-length(Tmp_df[Tmp_df$log2FoldChange < 0,]$gs_name)

# Preparing lists
if(Enriched_length != 0){
Tmp_df[Tmp_df$log2FoldChange > 0,]$gs_name<-paste("Enriched",i, sep = "_")
}

if(Depleted_length != 0){
Tmp_df[Tmp_df$log2FoldChange < 0,]$gs_name<-paste("Depleted",i, sep = "_")
}

Enriched_map[[i]]<-Tmp_df[Tmp_df$log2FoldChange > 0, c("gene_symbol","gs_name")]
Depleted_map[[i]]<-Tmp_df[Tmp_df$log2FoldChange < 0, c("gene_symbol","gs_name")]

}

# melting lists
GeneSet<-rbind(reshape2::melt(Enriched_map),
reshape2::melt(Depleted_map))

# returning GeneSet
return(GeneSet[c("gene_symbol","gs_name")])
}