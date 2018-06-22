#' Overlaps myDNA genotypes with SNP-trait database of interest
#'
#' @param myDNA (data.frame) analyzed genotypes, an output of importDNA function
#'
#' @param database (character) info database which stores info about genotyping
#' OPTIONS: "ebicat37","ebicat38", "CurrentGwascat", "taSNP"
#'
#' NOTE! if currentGWAS option used then input myDNA needs to be liftovered
#' to hg38
#'
#'
#' @details This function identifies overlap between my genotypes (myDNA object)
#' and database which stores info about SNP-trait associations. Currently,
#' GWASCatalog (hg19,hg38 as ebicat37,ebicat38, and makeCurrentGWASC). However,
#' in the future it will use info from SNPedia as well.#'
#'
#' @import stringr
#' @import dplyr
#' @import plyr
#' @examples
#' \dontrun{
#' myScreenDNA <- myDNAScreenDB(myDNA = myDNA,
#' database="ebicat37" )
#' }
#' @author Inga Patarcic
#' @export
myDNAScreenDB <- function(myDNA,
                         database="ebicat37"){

  # step 1: import database: GWASCatalog
  require(gwascat)
  require(plyr)
      if (database=="ebicat37"){

          # STEP 1 read DB
                    data(ebicat37)
                    db <- as.data.frame(ebicat37)

          # step 2: overlap myDNA and db
                    myDNA$SNPS <- myDNA$rsid
                    myDNAAdded <- join(myDNA,
                                     db,
                                     by="SNPS",
                                     type="inner")
        }


      if (database=="ebicat38"){

        # STEP 1 read DB
                  data(ebicat37)
                  db <- as.data.frame(ebicat37)

        # step 2: overlap myDNA and db
                  myDNA$SNPS <- myDNA$rsid
                  myDNAAdded <- join(myDNA,
                                   db,
                                   by="SNPS",
                                   type="inner")
      }


      if (database=="CurrentGwascat"){

        # STEP 1 read DB
                  db <- as.data.frame(makeCurrentGwascat())

        # step 2: overlap myDNA and db
                  myDNA$SNPS <- myDNA$rsid
                  myDNAAdded <- join(myDNA,
                                   db,
                                   by="SNPS",
                                   type="inner")
      }

      #if (database=="taSNP"){ data(taSNP); db <- as.data.frame(taSNP)


  # assessing if risk allele is present or not
  # maybe adjust for each database

         myDNAAdded <- identifyRiskAllele(myDNAAdded = myDNAAdded)

      return(myDNAAdded)


}


#' Overlaps myDNA genotypes with SNPs of interest
#'
#'
#' @param myDNA (data.frame) analyzed genotypes, an output of importDNA function
#'
#' @param snpsData (data.frame) with min 2 columns: rsID (SNPid) and risk allele
#' In that order!
#'
#' @import stringr
#' @import dplyr
#'
#' @details This function identifies overlap between my genotypes (myDNA object)
#' and data.frame which stores info SNPs (their unique ID -rsID, and risk allele
#' -in that order. This function is useful when user already has a list of SNPs
#' associated with a trait of interest.
#'
#' @examples
#' \dontrun{
#' # example - creating example for myDNAScreenSNPS()
#' # IDEA is as follows: extract rs and snp risk alleles, and extract those for
#'  coronary heart disease
#'  # filter CHD SNPS
#' dd.heart <- filter(db,
#'                    MAPPED_TRAIT=="coronary heart disease")
#'
#' snpsData <-   dd.heart$STRONGEST.SNP.RISK.ALLELE
#' # greping risk allele
#' riskAlleles <- gsub("[^\\-]*-([ATCG?])", "\\1",
#'                     snpsData)
#'
#' snpsData <- cbind(dd.heart$SNPS,riskAlleles)
#'
#'
#' # check my genome for given list of SNPs
#' myDNAScreenSNPS(myDNA, snpsData )
#'
#'  }
#' @export
myDNAScreenSNPS <- function(myDNA,
                            snpsData){


      # step 1. Adjusting column names
          colnames(snpsData) <- c("SNPS","riskAlleles")
          snpsData <- data.frame(snpsData,stringsAsFactors = F)

      # step 2: overlap myDNA and db
          myDNA$SNPS <- myDNA$rsid
          myDNAAdded <- join(myDNA,
                           snpsData,
                           by="SNPS",
                           type="inner")

      # step 3. identify if I have a risk or not
          myDNAAdded <- identifyRiskAllele(myDNAAdded = myDNAAdded,
                                            riskSNP = "riskAlleles")

          return(myDNAAdded)

}



#' f which adds info if i have a risk or not
#' @param riskSNP column name where risk SNP info is stores
#' @param myDNAAdded data.table with genotypes and added info from snp-diease DB
#'
#' @author Inga Patarcic
#' @keywords internal
identifyRiskAllele <- function(myDNAAdded,
                               riskSNP="STRONGEST.SNP.RISK.ALLELE"){

      # greping risk allele
      riskAlleles <- gsub("[^\\-]*-([ATCG?])", "\\1",
            myDNAAdded[,riskSNP])

    # assessing if analyzed genome has a risk allele
        myRisk <- mapply(
              function(risk,
                       mine) {
                risk %in% unlist(str_split(mine, ""))
                        },
                      riskAlleles,
                      myDNAAdded$genotype)

    # adding info if i have risk or not
        myDNAAdded$myAllele <- myRisk

                return(myDNAAdded)

}






#' Filter prescreened myDNA object filter for trait of interest
#'
#' @param myScreenDNA GRanges object, output of importDNA function
#'
#' @param trait (character) which trait are you interested in
#'
#' @param risk (character; default "high")
#' "high" -> OR>=1.5,
#'  "elavated" -> OR >1&<1.5
#'  "protective" -> OR >0.75&<1
#'  "high.protective" -> OR <0.75,
#'
#' @details This function allows you to filter your imported myDNA object +
#' overlapped with database of interest (GWASCatalog) by selecting trait in
#' which you are interested in, and selecting category of SNPs. For example,
#' high risk SNPs are those which have reported OR or beta >1.5.
#'
#' @import stringr
#' @import dplyr
#' @examples
#' \dontrun{
#'myAchille <- traitScreen(myScreenDNA, risk="high")
#'myDepression <- traitScreen(myScreenDNA, risk="all",
#'trait="unipolar depression")
#'myDepression.high <- traitScreen(myScreenDNA, risk="high",
#' trait="unipolar depression")
#' }
#' @author Inga Patarcic
#' @export
traitScreen <- function(myScreenDNA,
                        trait="all",
                        risk="high.protective"){

  require(dplyr)

  # select only those with risk
        myGenomeAtRisk <- myScreenDNA[myScreenDNA$myAllele,]

  # select trait

        if (trait!="all"){myGenomeAtRisk <- filter(myGenomeAtRisk,
                                                  MAPPED_TRAIT==trait)}

        if (risk=="high.protective"){ return(filter(myGenomeAtRisk,
                                                    OR.or.BETA<=0.75))}

        if (risk=="protective"){ return(filter(myGenomeAtRisk,
                                             (OR.or.BETA<1.0&OR.or.BETA>0.75)))}


        if (risk=="elavated"){ return(filter(myGenomeAtRisk,
                                             (OR.or.BETA>1.0&OR.or.BETA<=1.5)))}

        if (risk=="high"){ return(filter(myGenomeAtRisk,OR.or.BETA>1.5))}


        if (risk=="all"){ return(myGenomeAtRisk)}



}







