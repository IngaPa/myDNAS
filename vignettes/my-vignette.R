## ----global_options, include=FALSE---------------------------------------

library(knitr)
library(myDNA)
library(gwascat)
library(dplyr)
library(rmarkdown)
library(stringr)

opts_chunk$set(warning = FALSE,
               message= FALSE,
               fig.align='center',
               fig.path='Figures', 
               dev='png',
               fig.show='hold', 
               cache=FALSE)



## ---- eval = FALSE-------------------------------------------------------
#  
#  # example myHeritage
#  Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#  myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#  
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  myScreenDNA <- myDNAScreenDB(myDNA = myDNA,
#                               database="ebicat37" )
#  
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  
#  # IDEA is as follows: create SNPs objectby extracting rs and snp risk alleles, and extract those for coronary heart disease
#  
#  
#        # filter CHD SNPS
#        dd.heart <- filter(data.frame(ebicat37),
#                           MAPPED_TRAIT=="coronary heart disease")
#  
#        snpsData <-   dd.heart$STRONGEST.SNP.RISK.ALLELE
#        # greping risk allele
#        riskAlleles <- gsub("[^\\-]*-([ATCG?])", "\\1",
#                            snpsData)
#  
#        snpsData <- cbind(dd.heart$SNPS,riskAlleles)
#  
#  
#  # run myDNAScreenSNPS
#  myDNAScreenSNPS(myDNA = myDNA,
#                  snpsData = snpsData )
#  
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  # screen db for all traits
#  myAchille <- traitScreen(myScreenDNA, risk="high")
#  
#  # screen db for one trait for protective SNPA
#  myDepression <- traitScreen(myScreenDNA,
#                              risk="protective",
#                              trait="unipolar depression")
#  
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#   expandByLD(myDNA = myGenome[1:10],
#            pop="EUR",
#             R.squared=0.8,
#             chunks=1000,
#             out.dir="/data/akalin/Projects/AAkalin_myDNA/Results/")
#  

