# Hotgenes (deprecated repo)

If you read about this package or saw a talk about it, this is not the official repo. 


This package reads the DESeq2 or Limma output files to provide 
users with a way to visualize the results in their 
browser. Currently only supports DE with "Wald" testing. 
"LRT" will still work, but contrast selection has not been optimized. 

```
# Install Hotgenes

# For R 3.6.1, install XML accordingly:
install.packages("XML", type = "binary")

install.packages("devtools")
devtools::install_github("Rvirgenslane/Hotgenes")

# Download and try with example data!
library(Hotgenes)

Example_Hotgenes_dir<-system.file("extdata",
"Example_Hotgenes.Rdata",
package = "Hotgenes", mustWork = TRUE)

load(Example_Hotgenes_dir)
if(interactive()){
Shiny_DE_viz(Example_Hotgenes)}

# GSEA support
Example_Hotgenes_dir<-system.file("extdata",
"Example_Hotgenes.Rdata",
package = "Hotgenes", mustWork = TRUE)
load(Example_Hotgenes_dir)
library(msigdbr)

# GO annotations
m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# Reactome annotations
# m_df = msigdbr(species = "Homo sapiens", category = "C2", 
# subcategory = "CP:REACTOME")  

qbat<-BatchGSEA(HotgenesObj=Example_Hotgenes,
m_df= m_df)

# View enriched pathways
lapply(qbat$OuputGSEA, function(x) head(x$fgRes$pathway))

# plot of enriched pathways
qbat$OuputGSEA$shEWS$g

# Prep for GSEA plot
m_list = qbat$m_list

# top pathways
qbat$OuputGSEA$shEWS$top$pathway
pyid<-qbat$OuputGSEA$shEWS$top$pathway[[2]]
pyid

Gene_Ranks <- qbat$OuputGSEA$shEWS$Gene_Ranks
fgsea::plotEnrichment(m_list[[pyid]], Gene_Ranks,
gseaParam = 1, ticksSize = 0.2) +
ggplot2::labs(title=stringr::str_wrap(pyid, 20))


# Explore your own DESeq2 analysis:
Input_Hotgenes<-DEseq2_export(DEseq2_object = dds_con,
padj_cut = 0.1)

Shiny_DE_viz(Input_Hotgenes) # Visualize your results!

# Explore Limma DE analysis:
Example_Hotgenes_dir<-system.file("extdata",
"Example_Hotgenes.Rdata",
package = "Hotgenes", mustWork = TRUE)
load(Example_Hotgenes_dir)

library(limma)

exp<-Example_Hotgenes$Normalized_Expression$rld
design_m<-Example_Hotgenes$design_data

design_matrix <- model.matrix(~sh*Hrs+Bio_Rep,   
data = design_m)

aw <- arrayWeights(exp, design_matrix)

fit <- lmFit(exp, design=design_matrix, weights = aw)
fit <- eBayes(fit, robust = TRUE) 

L_out<-Limma_export(Expression_dat = exp, design_data = design_m, 
limmafit = fit)

if(interactive()){
Shiny_DE_viz(L_out)}
```