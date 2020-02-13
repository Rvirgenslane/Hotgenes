#' Exports hotgenes object to csv files
#' @export
#' @importFrom utils write.table
#' @param hotgenes object created by DEseq2_export or Limma_export
#' functions.
#' @param Exps_Out logical, if TRUE (default), available expression data
#' is also exported to csv file.
#' @param sep field separator string
#' @param readouts integer indicated the desired normalized data to be 
#' written to csv file. If 1 (default), the first method is exported
#' @param dir_out string providing the desired directory for export. 
#' If not provided, the working directory will be used.
#' @return multiple csv files
#' @examples
#' Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' hotgenes_csv(Example_Hotgenes)


# export to csv 
hotgenes_csv<-function(hotgenes=NULL,
Exps_Out=TRUE,
readouts=1,sep="\t",
dir_out=NULL)  {

if(is.null(dir_out)){
    dir_out<-file.path(getwd(),"hotgenes_csv")
    dir.create(dir_out)
}
# output path
dir.create(file.path(dir_out, "DE_tables"))
DE_path<-file.path(dir_out, "DE_tables")    

# DE tables
for(i in names(hotgenes$Output_DE)){
table_out<-hotgenes$Output_DE[[i]]
table_out$Genes <- rownames(table_out)
rownames(table_out)<-NULL
table_out<-table_out[c(7,1:6)]

write.table(table_out, 
quote = FALSE,sep=sep, row.names = FALSE,
file = file.path(DE_path, paste(i,".csv", sep = "")))
}



# output path
dir.create(file.path(dir_out, "Enriched_Genes"))
en_path<-file.path(dir_out, "Enriched_Genes")

# Enriched genes
for(i in names(hotgenes$Enriched_by)){
Enriched_by_table_out<-data.frame(hotgenes$Enriched_by[[i]],
check.names = FALSE)

names(Enriched_by_table_out)<-paste("Enriched_by_",i,sep="")

write.table(Enriched_by_table_out, 
quote = FALSE,sep=sep,
file = file.path(en_path, paste("Enriched_by_",i,".csv", sep = "")))
}



# output path
dir.create(file.path(dir_out, "Depleted_Genes"))
depleted_path<-file.path(dir_out, "Depleted_Genes")
# Depletion tables
for(i in names(hotgenes$Depleted_by)){
Depleted_by_table_out<-data.frame(hotgenes$Depleted_by[[i]],
check.names = FALSE)
names(Depleted_by_table_out)<-paste("Depleted_by_",i,sep="")

write.table(Depleted_by_table_out, 
quote = FALSE,sep=sep,
file = file.path(depleted_path, paste("Depleted_by_",i,".csv", sep = "")))
}


if(isTRUE(Exps_Out)){
# output path
dir.create(file.path(dir_out, "Expression_table"))
Exps_path<-file.path(dir_out, "Expression_table")    

Exps_table_out<-hotgenes$Normalized_Expression[[readouts]]
write.table(Exps_table_out, 
quote = FALSE, sep=sep,
file = file.path(Exps_path, paste("Expression",".csv", sep = "")))
}
}  


