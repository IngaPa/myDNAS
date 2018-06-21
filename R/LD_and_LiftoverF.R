#' For genotypes identify proxy SNPs (with R2 above threshold)
#'
#' @param myDNA GRanges object, output of importDNA function
#'
#' @param pop (character, default EUR) the name of a 1000 Genomes population
#' (AMR,AFR,ASN,EUR,...). Set this to NA to use all populations.
#'
#' @param  R.squared (numeric) threshold for R.squared statistics used to
#'  assess LD
#'
#' @param  chunks (numeric, default=1). Enables running function for N chunks
#'  of SNPs. It is to intensive to run it on the full set of genotypes (700k)
#'  so it is better to split it into cca 1000 chunks
#'
#' @param out.dir (character) path to the directory where data will be stored
#'
#' @details For analyzed genotypes (myDNA) this function identifies all SNPs
#' that in LD above predefined threshold. Since number of SNPs is large (>700k)
#' analysis is separated into chunks, and results are pooled together in the end
#' In addition, chunks and final results are saved in the out.dir
#'
#' @import proxysnps
#' @import parallel
#' @examples
#' \dontrun{}
#' @author Inga Patarcic
#' @export
expandByLD <- function(myDNA,
                       pop="CEU",
                       R.squared=0.8,
                       chunks=1,
                       out.dir){

  require(parallel)
  # split data into chunks
  Split.factor <- split(1:nrow(myDNA),
                        sort(1:nrow(myDNA)%%chunks))

  # Run LD in chunks
  tmp <- mclapply(Split.factor, function(chunk) {


    SNPs_LD <- mclapply(chunk,function(x,
                                       pop,
                                       R.squared,
                                       myDNA){

      print(x)
      require(proxysnps)


      LD_snps <-  get_proxies(chrom = myDNA$chrom[x],
                              pos = myDNA$position[x],
                              pop = pop)

      # filtering by LD
      LD_snps <- LD_snps[which(LD_snps$R.squared>=R.squared),]
      if (nrow(LD_snps)>0){ LD_snps$OrginalSNP <- myDNA$rsid[x]}

      return(LD_snps)

    },
    pop=pop,
    R.squared=R.squared,
    myDNA=myDNA,
    mc.cores=5)

    # binding data together
    SNP_LD_all <- do.call("rbind.data.frame",SNPs_LD)
    # removing duplicated entries
    dupl <- duplicated(SNP_LD_all$ID)
    SNP_LD_allSNP_LD_all <- SNP_LD_all[!dupl,]


    saveRDS(object = SNP_LD_allSNP_LD_all,
            file = paste0(out.dir,"/myDNA_LDchucks_",chunk[1],".rds"))


  },mc.cores = 5)


  # read-in back all LD data
  AllLDchunks <- list.files(out.dir,
                            pattern =  "myDNA_LDchucks_",
                            full.names = T)


  AllLDdata <- lapply(AllLDchunks,readRDS)
  AllLDdata <- do.call("rbind.data.frame",AllLDdata) #bind
  AllLDdata <- AllLDdata[!duplicated(AllLDdata$ID),] #!dupl

  saveRDS(object = AllLDdata,
          file = paste0(out.dir,"/myDNA_LD_full.rds"))

  return(AllLDdata)



}











# SNPs_LD <- lapply(1:10,
#                   function(x,
#                            pop,
#                            R.squared){
#
#
#    require(proxysnps)
#
#           LD_snps <-  get_proxies(chrom = myDNA$chrom[x],
#                           pos = myDNA$position[x],
#                           pop = pop)
#
#           # filtering by LD
#           LD_snps <- LD_snps[which(LD_snps$R.squared>=R.squared),]
#           LD_snps$OrginalSNP <- myDNA$rsid[x]
#
#     return(LD_snps)
#
# },pop=pop,R.squared=R.squared)
# install.packages("rsnps")
# #library(rsnps)
# #rsnps::LDSearch("rs420358")
#
# install.packages("devtools")
# devtools::install_github("slowkow/proxysnps",force=T)
# library(proxysnps)
