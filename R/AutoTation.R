#' Builds color annotations for pheatmap
#' @param Design Hotgenes object or a data.frame containing experimental design data 
#' @param VarIds Vector of Strings to select design variables. Default (NULL), returns
#' @param Colors vector for desired colors
#' @param binary_var vector for desired binary colors.
#' @export
#' @return List of colors matching design variables
#' @examples
#' Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' AutoTation(Example_Hotgenes, 
#' Colors=c("cyan", "magenta", "blue", "pink"),
#' binary_var=c("white", "black"))


AutoTation<-function(Design=NULL, VarIds=NULL, Colors=NULL,
binary_var=c("white", "black")){

Max_levels<-length(Colors)

if(!is.null(Design[["design_data"]])){
metadata<- Design[["design_data"]]

} else if(is.null(Design[["design_data"]])){
metadata<- Design
}


if(!is.null(VarIds)){
metadata<- metadata[VarIds]

}


annotation_col_frame <-Filter(is.factor, metadata)
annotation_col_frame <- annotation_col_frame[, sapply(annotation_col_frame, nlevels) > 1]
an_col<-names(annotation_col_frame)

an_col.list <- setNames(vector("list", length(an_col)), 
an_col)

for(i in names(an_col.list)){

length_var<-length(levels(annotation_col_frame[,i]))


if(length_var ==2){
Var1<-binary_var

names(Var1) <- levels(annotation_col_frame[,i])
an_col.list[[i]]<-Var1
}else if(length_var > 2 & length_var <= Max_levels)
Var1        <- Colors[1:length_var]
names(Var1) <- levels(annotation_col_frame[,i])
an_col.list[[i]]<-Var1
}


return(an_col.list)
}


