## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----style, echo = FALSE, results = 'asis'------------------------------------
library(BiocStyle)

## ----knitr, echo = FALSE------------------------------------------------------
library(knitr)

## ----Retrieve PPI network-----------------------------------------------------
library(PPIpriority)
#Set  extended timeout, not necessary while using the package
options(timeout=700)

network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
score_threshold=700, input_directory=getwd())

## ----Run prioritization-------------------------------------------------------
output_df<-prioritize_targets(network_ppi_obj, seeds = "ERP44")

head(output_df,15) # show the df

## ----Run prioritization with ENTREZ gene id-----------------------------------
output_df_erp44_entrez<-prioritize_targets(network_ppi_obj, seeds = "23071")

# is the output equal to out previous dataset obtained with HUGO symbols?
all.equal(output_df_erp44_entrez,output_df)

## ----Run prioritization, multiple seeds and candidates------------------------
seeds_genes<-c("ERP44","LMAN1", "MCFD2")
candidate_genes<-c("AREG", "IL1A", "IL1B", "MMP3", "PRL", "F5")

output_df_msmc<-prioritize_targets(network_ppi_obj, 
                              seeds = seeds_genes,
                              candidates=candidate_genes)

# show the df
output_df_msmc

## ----Access netowrk information-----------------------------------------------
get_network(network_ppi_obj)

## ----Access STRINGdb obj------------------------------------------------------
get_STRINGdb_connection(network_ppi_obj)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

