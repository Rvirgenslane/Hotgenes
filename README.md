# Hotgenes
This package reads the DESeq2 output file to provide 
users with a way to visualize the results in their 
browser. Currently only supports DE with "Wald" testing. 
"LRT" will still work, but contrast selection has not been optimized. 

In the Shiny_DE_viz app:
-Tab 1 integrates FactoMineR's HCPC function to 
detect relationships between conditions and differentially expressed genes
--HCPC can help find the most representation DE genes for validation 
and export them to csv file.

-Tab 2 Reports the DE coefficients calculated by DESeq2
--DE coefficients can also be exported to csv file.

-Tab 3 integrates Pheatmap
--Representative Genes identified in tab 1 and saved to 
a csv file can be imported for visualization.
--Alternatively, Genes can be identified by the selecting the listed contrasts.
--Heatmap can be customized and exported to PDF.

# Download and try with example data!
library(Hotgenes)

Example_Hotgenes_dir<-system.file("extdata",
"Example_Hotgenes.Rdata",
package = "Hotgenes", mustWork = TRUE)

load(Example_Hotgenes_dir)
if(interactive()){
  Shiny_DE_viz(Example_Hotgenes)}

# have your own DESeq2 file to explore:
Input_Hotgenes<-DEseq2_export(DEseq2_object = dds_con,
  padj_cut = 0.1)
              
Shiny_DE_viz(Input_Hotgenes) # Visualize your results!
