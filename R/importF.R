#' Function to import myHeritage, 23andMe, AncestryDNA, Illumina personal genotypes
#'
#' @param myGenotypes (character) path to the genotype file
#'
#' @param type (character) info about the company used to genotype data: OPTIONS
#' "myHeritage", "23andMe", "AncestryDNA", "Illumina"
#'
#' @param asGRanges (default FALSE) Result of the function is GRanges or
#' data.frame object
#'
#' @details Import genotypes into R
#'
#' @import GenomicRanges
#' @import stringr
#' @examples
#'
#' \dontrun{
#' # example myHeritage
#' library(GenomicRanges)
#' library(stringr)
#' myGenotypes="D:/GoogleDrive/Projects/myDNA/Data/MyHeritage/MyHeritage_raw_dna_dataInga/MyHeritage_raw_dna_data.csv"
#' myDNA <- importDNA(myGenotypes)
#' }
#' @author Inga Patarcic
#' @export
importDNA <- function(myGenotypes,
                      type="myHeritage",
                      asGRanges=FALSE) {

  if (type=="myHeritage"){

        myDNA <- read.table(myGenotypes,he=F, skip=7,sep=",",
                              col.names=c("rsid", "chrom",
                                          "position", "genotype"))


        }


  if (type=="23andMe"){

        myDNA <- read.table(myGenotypes,he=T, skip=15,sep="\t",
                                    col.names=c("rsid", "chrom",
                                        "position", "genotype"))

  }

  if (type=="AncestryDNA"){

        myDNA <- read.table(myGenotypes,he=T, skip=17,sep="\t",
                            col.names=c("rsid", "chrom",
                                        "position", "genotype","genotype2"))

        # Ancestry reports genotypes separately so i need to merge info
        myDNA$genotype <- paste0(myDNA$genotype,myDNA$genotype2)
        myDNA <- myDNA[,-5]

        }

 # if (type=="FamilyTree"){ myDNA <-  read.table(myGenotypes,he=T, skip=15,sep="\t")}


   # returns GRanges object

  if (asGRanges==TRUE){ return(asGR(myDNA))}
  if (asGRanges==FALSE){ return(myDNA)}



}







#' Creates Granges object from genotyping results
#'
#' @param myDNA data.frame. Results of read.data for genotyping results
#'
#' @param type (character) info about the company used to genotype data:
#' myHeritage, 23andMe, AncestryDNA, Illumina
#'
#' @details outputs granges
#'
#' @import GenomicRanges
#'
#' @author Inga Patarcic
#' @keywords internal
#' @export
asGR <- function(myDNA){

  GR <- GRanges(paste0("chr",myDNA$chrom),
                IRanges(myDNA$position,
                        myDNA$position))

  GR$rsid <- myDNA$rsid
  GR$genotype <- myDNA$genotype

  return(GR)

}




