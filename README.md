# Hotgenes
This package reads the DESeq2 output file to provide 
users with a way to visualize the results in their 
browser. Currently only supports DE with "Wald" testing. 
"LRT" will still work, but contrast selection has not been optimized. 

In the Shiny_DE_viz app:

Tab 1 Normalization QC:

-Visualize your data as normalized by different methods

Tab 2 Find Hotgenes:

-FactoMineR's HCPC function can find relationships between conditions 
and differentially expressed genes
-Export genes to csv file

Tab 3 DE coefficients:

-DE coefficients calculated by DESeq2 can be viewed

Tab 4 Heatmap:

-Representative Genes identified in tab 1 and saved to 
a csv file can be imported for visualization.
-Alternatively, Genes can be identified by the selecting the listed contrasts.
-Heatmap can be customized and exported to PDF.

# Install Hotgenes
        install.packages("devtools")
        library(devtools)
        install_github("Rvirgenslane/Hotgenes")

# Download and try with example data!
    library(Hotgenes)

    Example_Hotgenes_dir<-system.file("extdata",
    "Example_Hotgenes.Rdata",
    package = "Hotgenes", mustWork = TRUE)

    load(Example_Hotgenes_dir)
    if(interactive()){
      Shiny_DE_viz(Example_Hotgenes)}

# Explore your own DESeq2 analysis:
    Input_Hotgenes<-DEseq2_export(DEseq2_object = dds_con,
      padj_cut = 0.1)

    Shiny_DE_viz(Input_Hotgenes) # Visualize your results!
