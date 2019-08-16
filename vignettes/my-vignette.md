---
title: "myDNA"
author: "Inga Patarcic"
date: "2019-08-16"
output:  
  html_document:
        toc: true
        toc_float: true
        number_sections: true
        toc_depth: 3
vignette: >
  %\VignetteIndexEntry{myDNA}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




# Introduction

This package is used to analyze personal genotype which is obtained by using one of the direct-to-consumer (DTC) DNA tests from DNA testing companies such as myHeritageDNA, 23andMe, AncestryDNA, etc.
Currently, user can:
1) import genotypes, 
2a) overlap personal genotypes with SNP-triat database results 
(such as GWASCatalog) or,
2b) a list od SNP unique identifiers and corresponding risk alleles,
3) identify personal risk alleles, 
4) expand SNP list by identifying SNPs that are in LD with reported SNPs. 




-----------------------------------------------------------------------------


# Import genotypes using importDNA()

One can import personal genotypes obtained by using one of the
direct-to-consumer (DTC) DNA tests from DNA testing companies: myHeritageDNA, 23andMe, AncestryDNA


```r
# example myHeritage
Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
```




# Overlap genotypes with SNP-trait database with myDNAScreenDB()

Once can overlap between my genotypes (myDNA object) and database which stores 
info about SNP-trait associations. Currently available only GWASCatalog
(hg19,hg38 as ebicat37,ebicat38, and makeCurrentGWASC). However,
in the future it will use info from SNPedia


```r
myScreenDNA <- myDNAScreenDB(myDNA = myDNA,
                             database="ebicat37" )
```






#Overlaps myDNA genotypes with SNPs of interest myDNAScreenSNPS()

Overlap my genotypes (myDNA object) and data.frame which stores info about SNPs (their unique ID -rsID, and risk allele -in that order. This function is useful when user already has a list of SNPs associated with a trait of interest.


```r
# IDEA is as follows: create SNPs objectby extracting rs and snp risk alleles, and extract those for coronary heart disease

      # filter CHD SNPS
      dd.heart <- filter(data.frame(ebicat37),
                         MAPPED_TRAIT=="coronary heart disease")
      
      snpsData <-   dd.heart$STRONGEST.SNP.RISK.ALLELE
      # greping risk allele
      riskAlleles <- gsub("[^\\-]*-([ATCG?])", "\\1",
                          snpsData)
      
      snpsData <- cbind(dd.heart$SNPS,riskAlleles)


# run myDNAScreenSNPS
myDNAScreenSNPS(myDNA = myDNA,
                snpsData = snpsData )
```






# Filter prescreened myDNA object for trait traitScreen()

This function allows you to filter your imported myDNA object +
overlapped with database of interest (GWASCatalog) by selecting trait in which you are interested in, and selecting category of SNPs. For example,high risk SNPs are those which have reported OR or beta >1.5.


```r
# screen db for all traits
myAchille <- traitScreen(myScreenDNA, risk="high")

# screen db for one trait for protective SNPA
myDepression <- traitScreen(myScreenDNA,
                            risk="protective",
                            trait="unipolar depression")
```



#Identify proxy SNPs for all SNPs using expandByLD()


For analyzed genotypes (myDNA) this function identifies all SNPs
 in LD above predefined threshold. 
 
 

```r
 expandByLD(myDNA = myGenome[1:10],
          pop="EUR",
           R.squared=0.8,
           chunks=1000,
           out.dir="/data/akalin/Projects/AAkalin_myDNA/Results/")
```
