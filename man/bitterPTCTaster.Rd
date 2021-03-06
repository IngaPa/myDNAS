% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpleTests.R
\name{bitterPTCTaster}
\alias{bitterPTCTaster}
\title{Function to screen for PTC bitter tasting ability}
\usage{
bitterPTCTaster(myDNA)
}
\arguments{
\item{myDNA}{(character) path to the genotype file}
}
\description{
Function to screen for PTC bitter tasting ability
}
\details{
SOURCE: 23andMe report

INFO:
  rs713598 - a genetic marker that affects your chances of being able to detect
a certain bitter chemical called "PTC." Some vegetables like raw broccoli
and brussels sprouts, contain bitter chemicals similar to PTC.
CC - unable to detect, GG/CG - able to detect.
The TAS2R38 gene contains instructions for a protein, or taste receptor,
that can detect the bitter chemical called "PTC."

NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
can be imputed and accessed by DNALand imputation tool

!!!! myAllele column reports whether risk allele was identified in my genome
}
\examples{

\dontrun{
# example myHeritage
 library(myDNA)
 library(plyr)
 library(stringr)
  Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
   myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
 bitterPTCTaster(myDNA)
}
}
\author{
Inga Patarcic
}
