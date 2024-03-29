#' Generates a table of quantitative variables from FactoMineR data 
#' @export
#' @param res.hcpc R object generated by HCPC function
#' from FactoMineR
#' @param Round_To numeric value indicating digit for rounding
#' default = 2
#' @param Format Logical, if TRUE (default), output formated by
#' group.
#' @return data.frame with HCPC results
#' @examples Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' Results<-FactoWrapper(Example_Hotgenes)
#' Quanti_Table(Results$res.hcpc)

# FactoMiner Functions ----------------------------------------------------
Quanti_Table<-function(res.hcpc=NULL, Round_To = 2,
                       Format = TRUE){
    
# Table_Express 
    Table_Express<-function(df_cat){
        EqDatedf_cat <- as.data.frame(df_cat[1,])
        
        cat_cols<-colnames(df_cat)
        df_cat_lists <- setNames(vector(length(cat_cols), mode="list"), cat_cols)
        
        for (i in cat_cols){
            df_cat_lists[[i]]<-("")
        }
        EmptyLine <- data.frame(df_cat_lists, check.names = FALSE)
        
        
        for (i in 2:nrow(df_cat)){
            if(as.vector(df_cat$Groups[i])  ==  as.vector(df_cat$Groups[i-1])){
                EqDatedf_cat <- rbind(EqDatedf_cat, df_cat[i,])
            } else {
                EqDatedf_cat <- rbind(EqDatedf_cat, EmptyLine)
                EqDatedf_cat <- rbind(EqDatedf_cat, df_cat[i,])
            }
        }
        
        rownames(EqDatedf_cat)<-NULL
        return(EqDatedf_cat)
    }
    
    
    
    # Variable
    df_quanti <- plyr::ldply(res.hcpc$desc.var$quanti, data.frame)
    df_quanti$Variable<- names(do.call(rbind,res.hcpc$desc.var$quanti)[,1])
    df_quanti<-df_quanti[c(1,8,2:7)]
    colnames(df_quanti)<-gsub(".id", "Groups", fixed = TRUE, colnames(df_quanti))
    
    # round only numeric columns
    numeric_colnames<-colnames(dplyr::select_if(df_quanti, is.numeric))
    df_quanti[,numeric_colnames]<-round(df_quanti[,numeric_colnames], Round_To)
    
    df_quanti<-df_quanti[,seq(1,5)]
    df_quanti$Variable<-gsub(".", " ", df_quanti$Variable, fixed = TRUE)
    colnames(df_quanti)<-gsub(".", " ",colnames(df_quanti), fixed = TRUE)
    colnames(df_quanti)<-gsub("Variable", "Genes",colnames(df_quanti), fixed = TRUE)
    
    if(Format == TRUE){
        Table_df_quanti<-Table_Express(df_quanti)
        return(Table_df_quanti)
    }else if(Format != TRUE){
    return(df_quanti)
    }    
}



#' Generates a table of Categorical variables from FactoMineR data 
#' @export
#' @param res.hcpc R object generated by HCPC function
#' from FactoMineR
#' @param Round_To numeric value indicating digit for rounding
#' default = 2
#' @param Format Logical, if TRUE (default), output formated by
#' group.
#' @return data.frame with HCPC results
#' @examples Example_Hotgenes_dir<-system.file("extdata",
#' "Example_Hotgenes.Rdata",
#' package = "Hotgenes", mustWork = TRUE)
#' load(Example_Hotgenes_dir)
#' Results<-FactoWrapper(Example_Hotgenes)
#' Categorical_Table(Results$res.hcpc)

# Categorical table ------------------------------------------------------------
Categorical_Table<-function(res.hcpc = NULL, Round_To = 2,
                            Format=TRUE){
    
    # Table_Express 
    Table_Express<-function(df_cat){
        EqDatedf_cat <- as.data.frame(df_cat[1,])
        
        cat_cols<-colnames(df_cat)
        df_cat_lists <- setNames(vector(length(cat_cols), mode="list"), cat_cols)
        
        for (i in cat_cols){
            df_cat_lists[[i]]<-("")
        }
        EmptyLine <- data.frame(df_cat_lists, check.names = FALSE)
        
        
        for (i in 2:nrow(df_cat)){
            if(as.vector(df_cat$Groups[i])  ==  as.vector(df_cat$Groups[i-1])){
                EqDatedf_cat <- rbind(EqDatedf_cat, df_cat[i,])
            } else {
                EqDatedf_cat <- rbind(EqDatedf_cat, EmptyLine)
                EqDatedf_cat <- rbind(EqDatedf_cat, df_cat[i,])
            }
        }
        
        rownames(EqDatedf_cat)<-NULL
        return(EqDatedf_cat)
    }
    
    
    # Categories
    df_cat <- plyr::ldply(res.hcpc$desc.var$category, data.frame)
    colnames(df_cat)<-gsub(".id", "Groups", fixed = TRUE, colnames(df_cat))
    df_cat$Category<- names(do.call(rbind,res.hcpc$desc.var$category)[,1])
    
    df_cat$Category<-gsub("_"," ",fixed = TRUE,
                          df_cat$Category)
    
    df_cat$Category<-gsub("=",": ",fixed = TRUE,
                          df_cat$Category)
    df_cat<-df_cat[c(1,7,2:6)]
    
    # round numeric columns
    numeric_colnames<-colnames(dplyr::select_if(df_cat, is.numeric))
    df_cat[,numeric_colnames]<-round(df_cat[,numeric_colnames], Round_To)
    
    df_cat$p.value<-NULL
    colnames(df_cat)<-gsub(".", " ", colnames(df_cat), fixed = TRUE)
    col_percent<-colnames(df_cat[c(3:5)])
    
    for (i in col_percent) {
        df_cat[,c(i)]<-paste(df_cat[,c(i)], "%", sep = "")
    }
    
    if(Format == TRUE){
        Table_df_cat<-Table_Express(df_cat)
        return(Table_df_cat)
    }else if(Format != TRUE){
        return(df_cat)
    }   
    
}

