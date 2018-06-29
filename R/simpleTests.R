

##########################################
# --- screen for lactose intolerance ----
##########################################



#' Function to screen for lactose intolerance
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details SOURCE: 23andMe report
#' INFO:
#' Lactose intolerance is influenced by a genetic marker near the LCT gene.
#' GG - likely lactose intolerant. CG/CC - likely not lactose intolerant
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
# library(myDNA)
# library(dplyr)
# library(stringr)
# Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
# myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
# lactoseIntolerance(myDNA)
#' }
#' @author Inga Patarcic
#' @export
lactoseIntolerance <- function(myDNA) {

      SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                snpsData = cbind("rs4988235","G")))



  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}

}




###########################################
# --- screen for celiac disease ---
###########################################

#' Function to screen for celiac disease based on 23andMe tested SNPs
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SCREEN 1: rs2187668 (T) and rs7454108 (G)
#' SOURCE: https://blog.23andme.com/health-traits/new-23andme-report-celiac-disease/
#'
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'   celiacDisease23AndMe(myDNA)
#' }
#' @author Inga Patarcic
#' @export
celiacDisease23AndMe <- function(myDNA) {

    SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                         snpsData = cbind(c("rs2187668",
                                            "rs7454108"),
                                          c("T","G"))))



  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}

}


#' Function to screen for celiac disease based on 6 CD associated SNPs
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#'  SOURCE: https://www.mygenefood.com/genetics-celiac-disease-need-know/
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2386975/
#'   DQ2.5	rs2187668	T (0.09)
#' DQ8	rs7454108	G (0.18)
#' DQ2.2	rs2395182	T (0.71)
#' rs7775228	G (0.10)
#' rs4713586	G (0.025)
#' DQ7	rs4639334	A (0.09)
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'   celiacDiseaseExtended(myDNA)
#' }
#' @author Inga Patarcic
#' @export
celiacDiseaseExtended <- function(myDNA){


  SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                         snpsData = cbind(c("rs2395182",
                                            "rs7775228",
                                            "rs4713586",
                                            "rs2187668",
                                            "rs4639334" ,
                                            "rs7454108"),
                                          c("T","G","G","T","A","G"))))


  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}

}











###########################################
# --- screen for coffee consumption   ---
###########################################



#' Function to screen for higher coffee consumption
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE: 23andMe report
#' INFO:
#'   Report is based on genetic variants near two genes that play a role in how
#' your body handles caffeine. The first gene, CYP1A2, contains instructions for
#' an enzyme that breaks down 95% of the caffeine you consume.
#' The second gene, AHR (rs4410790), contains instructions for a protein that ramps
#' up production of the CYP1A2 enzyme. Variants in these genes may affect
#' how quickly the body breaks down and clears away caffeine.
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#' The CYP1A2 gene (rs2472297) contains instructions for an enzyme that breaks
#' down many substances, including caffeine.
#' This enzyme is a member of a large family of enzymes called cytochrome P450.
#' Presence of variants: likely higher coffee consumption.
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'   cofeeConsumption(myDNA)
#' }
#' @author Inga Patarcic
#' @export
cofeeConsumption <- function(myDNA){

                  SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                                        snpsData = cbind(c("rs2472297",
                                                           "rs4410790"),
                                                         c("T","C"))))


  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}

}

###########################################
# --- screen for deep sleep  ---
###########################################


#' Function to screen for deep sleep
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE: 23andMe report
#' INFO from 23andMe:
#' The genetic marker in this report is in the ADA gene,
#' which contains instructions for an enzyme that helps control adenosine levels.
#' Scientists think that adenosine builds up more quickly in people with one or
#' two copies of the T variant at this marker. This extra adenosine increases
#' sleepiness, leading to stronger delta waves. Because of this stronger sleep
#' pressure, people with the T variant also report feeling sleepier than other
#' people after a night of missed sleep.
#'
#'
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'   deepSleeper(myDNA)
#' }
#' @author Inga Patarcic
#' @export
deepSleeper <- function(myDNA){

    SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                             snpsData = cbind("rs73598374","T")))

  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}
}



###########################################
# --- screen for Muscle Composition  ---
###########################################


#' Function to screen for production of ACTN3 gene
#'
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE: 23andMe report
#' INFO:
#' A genetic marker in the ACTN3 gene.
#' This marker controls whether muscle cells produce a protein
#' (called alpha-actinin-3) that's found in fast-twitch muscle fibers.
#' While some people don't produce this protein at all, almost all of the elite
#' power athletes who have been studied have a genetic variant that allows them to
#' produce the protein. This suggests that the protein may be beneficial at least
#' at the highest levels of power-based athletic competition.
#'
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'  muscleComposition(myDNA)
#' }
#' @author Inga Patarcic
#' @export
muscleComposition <- function(myDNA){

  SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                             snpsData = cbind("rs1815739","C")))

  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}

}




###########################################
# --- screen for PTC bitter taste  ---
###########################################



#' Function to screen for PTC bitter tasting ability
#'
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE: 23andMe report
#'
#' INFO:
#'   rs713598 - a genetic marker that affects your chances of being able to detect
#' a certain bitter chemical called "PTC." Some vegetables like raw broccoli
#' and brussels sprouts, contain bitter chemicals similar to PTC.
#' CC - unable to detect, GG/CG - able to detect.
#' The TAS2R38 gene contains instructions for a protein, or taste receptor,
#' that can detect the bitter chemical called "PTC."
#'
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'  bitterPTCTaster(myDNA)
#' }
#' @author Inga Patarcic
#' @export
bitterPTCTaster <- function(myDNA) {

  SCREEN <- (myDNAScreenSNPS(myDNA = myDNA,
                                snpsData = cbind("rs713598","G")))

  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}
  }






###########################################
# --- screen for blue eyes  ---
###########################################


#' Function to screen for production of ACTN3 gene
#'
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE 1: 23andMe report
#' INFO:
#'   The genetic marker in this report is located near a gene called OCA2 that
#' affects how much brown pigment your cells produce. People with 1 or 2 copies
#' of the A variant of this marker tend to have more brown pigment in their eyes,
#' so they are likely to have darker eyes.
#' AA/AG - likely brown/hazel eyes, GG - likely blue/green eyes
#' NOTE! MyHeritage genotypes do not report this SNP, however missing genotypes
#' can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'  blueEyes(myDNA)
#' }
#' @author Inga Patarcic
#' @export
blueEyes <- function(myDNA) {

  SCREEN <- myDNAScreenSNPS(myDNA = myDNA,
                            snpsData = cbind("rs12913832","G"))

  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}
}







#' Function to screen for testing eye color
#'
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' Source 2: SNPedia for rs12913832
#' resuts source 2: rs12913832 is also part of a haplotype spanning
#' 166kB on chromosome 15, defined by 13 SNPs listed below,
#' that is found in 97% of all Caucasians with blue eyes.
#' The "h-1" haplotype found in homozygous state in 97% of individuals
#' with blue eye color is composed as follows [PMID 18172690]:
#'
#' rs4778241(C)
#' rs1129038(A)
#' rs12593929(A)
#' rs12913832(G)  OCA2/HERC2 has a score of 51.5 and an estimated allelic OR of 8.43
#' rs7183877(C)
#' rs3935591(G)
#' rs7170852(A)
#' rs2238289(T)
#' rs3940272(C)
#' rs8028689(T)
#' rs2240203(A)
#' rs11631797(G)
#' rs916977(G)
#' ADDED from SNPedia: rs1667394 increases susceptibility to Blond rather than
#' brown hair 4.94 times for carriers of the A allele
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'  blueEyesExtended(myDNA)
#' }
#' @author Inga Patarcic
#' @export
blueEyesExtended <- function(myDNA) {

  SCREEN <- myDNAScreenSNPS(myDNA = myDNA,
                            snpsData = cbind(c("rs1667394","rs4778241",
                                             "rs1129038",
                                             "rs12593929",
                                             "rs12913832",
                                             "rs7183877",
                                             "rs3935591",
                                             "rs7170852",
                                             "rs2238289",
                                             "rs3940272",
                                             "rs8028689",
                                             "rs2240203",
                                             "rs11631797",
                                             "rs916977"),
                            c("A","C","A","A","G","C","G","A","T","C","T","A","G","G")))


  if (nrow(SCREEN)==0) {
    stop("SNP not present in the genotype data, try with imputed dataset")}
  if (nrow(SCREEN)!=0) {return(SCREEN)}

}


###########################################
# --- screen for Neanderthal SNPs  ---
###########################################



#' Function to screen for certain Neanderthal variants
#'
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE 1: 23andMe report
#' If you have certain Neanderthal variants, it means that some of your physical
#' traits may trace back to your Neanderthal ancestors.
#'
#' Tested variants/gene/trait:
#' rs7544462/Gene: MEAF6/Trait: Height
#' rs1877547/Gene: LPP/Trait: Height
#' rs4849721/Gene: Near the EN1 gene/Trait: Less back hair
#' rs12458349/Gene: Near the PHLPP1 gene/Trait: Straight hair
#' rs11213819/Gene: Near the C11orf53 gene/Trait: Less likely to
#' sneeze after eating dark chocolate
#'
#' NOTE! MyHeritage genotypes do not report some if the SNPs, however missing
#' genotypes can be imputed and accessed by DNALand imputation tool
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'  neanderthalMe(myDNA)
#' }
#' @author Inga Patarcic
#' @export
neanderthalMe <- function(myDNA) {

  SCREEN <- myDNAScreenSNPS(myDNA = myDNA,
                               snpsData = cbind(c("rs7544462",
                                                  "rs1877547",
                                                  "rs4849721",
                                                  "rs12458349",
                                                  "rs11213819"),
                                                c("C","A","T","G","T")))



    if (nrow(SCREEN)==0) {
        stop("SNP not present in the genotype data, try with imputed dataset")}
    if (nrow(SCREEN)!=0) {return(SCREEN)}

}




###########################################
# --- screen for red hair  ---
###########################################


#' Function to screen for  red hairiness
#'
#'
#' @param myDNA (character) path to the genotype file
#'
#' @details
#' SOURCE 1: 23andMe report
#' NOTE! MyHeritage genotypes do not report some if the SNPs, however missing
#' genotypes can be imputed and accessed by DNALand imputation tool
#' INFO:
#' Several variants in a single gene, MC1R, can cause red hair by increasing
#' the amount of pheomelanin in your hair.
#' Hair color is determined by not just how much pigment you have,
#' but also what kind. Whether you have red hair depends on your levels of
#' red/yellow pigment (pheomelanin). But the lightness or darkness of your hair
#' depends on your levels of a brown/black pigment called eumelanin
#' REPORTED SNPS: rs1805007,rs1805008,i3002507
#'
#' SOURCE 2: https://www.snpedia.com/index.php/I3002507
#' i3002507 alias 	rs1805009 identified in SNPedia
#' proxy calculations needed
#'
#' !!!! myAllele column reports whether risk allele was identified in my genome
#'
#' @import stringr
#'
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#'  library(myDNA)
#'  library(dplyr)
#'  library(stringr)
#'   Genome="/data/akalin/Projects/AAkalin_myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#'    myDNA <- importDNA(myGenotypes = Genome,type = "myHeritage" )
#'  redHair(myDNA)
#' }
#' @author Inga Patarcic
#' @export
redHair <-  function(myDNA) {


  SCREEN <- suppressWarnings(myDNAScreenSNPS(myDNA = myDNA,
                            snpsData = cbind(c("rs1805007",
                                               "rs1805008",
                                               "rs1805009"),
                                             c("T","T","C"))))


        if (nrow(SCREEN)==0) {
          stop("SNP not present in the genotype data, try with imputed dataset")}
        if (nrow(SCREEN)!=0) {return(SCREEN)}

}

