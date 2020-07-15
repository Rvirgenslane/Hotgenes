## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Hotgenes)
Example_Hotgenes_dir<-system.file("extdata",
"Example_Hotgenes.Rdata",
package = "Hotgenes", mustWork = TRUE)
 load(Example_Hotgenes_dir)

 FactoWrapper(Example_Hotgenes)



