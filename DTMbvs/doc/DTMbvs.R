## ----setup, include=FALSE------------------------------------------------
library(RefManageR)
options(knitr.table.format = 'markdown')
file.name <- system.file("Bib", "microB.bib", package="RefManageR")
bib <- ReadBib("/Volumes/USB/mkoslovsky_mac/Rice Post Doc/Microbiome/Code/DTM/microB.bib")
BibOptions(style = "markdown", bib.style = "numeric",cite.style = "numeric")

## ----message=FALSE,  echo = FALSE----------------------------------------
library(mvtnorm)
library(MCMCpack)
library(Rcpp)
library(RcppArmadillo)
library(ape) # For handling phylogenetic tree data
library(GGMselect)
library(ggplot2)
sourceCpp("/Volumes/USB/mkoslovsky_mac/Rice Post Doc/Microbiome/Code/DTM/DTMbvs/src/DTMbvs_Final.cpp")
source("/Volumes/USB/mkoslovsky_mac/Rice Post Doc/Microbiome/Code/DTM/DTMbvs/R/DTMbvs_R.R")
source("/Volumes/USB/mkoslovsky_mac/Rice Post Doc/Microbiome/Code/DTM/DTMbvs/R/selected.R")
source("/Volumes/USB/mkoslovsky_mac/Rice Post Doc/Microbiome/Code/DTM/DTMbvs/R/bfdr.R")

## ----message=FALSE-------------------------------------------------------
DTMbvs_R( eval = FALSE )

## ----message = FALSE, echo = FALSE, fig.width = 7, fig.height= 4---------
tree.model$tip.label <- c(1,2,3,4,5)
plot( tree.model,use.edge.length = FALSE  )
nodelabels( )
edgelabels(c( "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"))

## ----message = FALSE-----------------------------------------------------
head( Y.model )

## ------------------------------------------------------------------------
DTM <- "micro( ( 1, 2 ), ( 3, ( 4, 5 ) ) );"
cat( DTM, file = "dtm.tre", sep = "\n" )
tree_DTM <- read.tree( "dtm.tre" )

## ------------------------------------------------------------------------
DM <- "micro( 1, 2, 3, 4, 5 );"
cat( DM, file = "dm.tre", sep = "\n" )
tree_DM <- read.tree( "dm.tre" )

## ----message=FALSE, results = "hide"-------------------------------------
test <- DTMbvs_R( tree = tree.model, Y = Y.model, X = X.model, aa = 0.1, bb = 0.9)

## ----message=FALSE-------------------------------------------------------
str(test)

## ----message = FALSE-----------------------------------------------------
 tree.model$edge

## ----message=FALSE, fig.width = 7, fig.height= 4-------------------------
edge_lab <- c( "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8")
selections <- selected( dtm_obj = test, threshold = c( 0.5 ), burnin = 2500, plotting = TRUE, edge_lab = edge_lab )

## ----echo = FALSE, results = 'asis'--------------------------------------
library(knitr)
kable(selections$selected_name, align = 'c')

## ----echo = FALSE, results = 'asis'--------------------------------------
covariate_list <- floor( ( truth - 1 )/8 ) + 1
branch_list <- ( truth - 1 ) %% 8 + 1
edge_lab <- c( "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8")
cov_lab <- seq(1:50)
R <- cbind( edge_lab[ branch_list ], cov_lab[ covariate_list ] )
colnames(R) <- c( "Branch", "Covariate" )
kable( R, align = "c") 

## ----message=FALSE, results = "hide"-------------------------------------
G_example <- matrix( 0, 50 ,50 )
G_example[ 1:10, 1:10 ] <- 1
diag( G_example ) <- 0 
test_G <- DTMbvs_R( iterations = 50000, tree = tree.model, Y = Y.model, X = X.model, prior = "MRF_fixed", G = G_example )

## ----message=FALSE, fig.width = 7, fig.height= 4-------------------------
selections_G <- selected( dtm_obj = test_G, threshold = c( 0.5 ), burnin = 2500, plotting = TRUE, edge_lab = edge_lab )

## ----echo = FALSE, results = 'asis'--------------------------------------
library(knitr)
kable(selections_G$selected_name, align = 'c')

## ----message=FALSE-------------------------------------------------------
bfdr_threshold <- bfdr( selections$mppi_zeta )$threshold
bfdr_threshold
bfdr_selected <- selected( test, burnin = 2500, threshold = bfdr_threshold, edge_lab = edge_lab)


## ----echo = FALSE, results = 'asis'--------------------------------------
library(knitr)
kable(bfdr_selected$selected_name, align = 'c')

## ----echo = FALSE, results = 'asis'--------------------------------------
PrintBibliography(bib)

