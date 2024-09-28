
#' Run up to 4 different copy number segmentation algorithms.
#'
#' \code{QDNAseqSegtoCNAclinic} Converts segmented copy number data from QDNAseq
#' into a CNAclinicData object, & calls gains & losses for downstream plotting 
#' and visualisation. Assumes CBS algorithm was used for segmentation and 
#' imports segmentation data into both {segCBS} and {segSummary} dataframes 
#' within the CNAclinicData object
#'
#' @param x A \code{\link{QDNAseqCopyNumbers}} object with copy number and 
#' segmentation data.
#' @param callThreshLog2R Thresholds used to call segments as a "loss" or "gain".
#' Defaults to c(-0.3, 0.3)
#'
#' @return Returns an object of class \code{\link{CNAclinicData}} with
#' segmentation results from QDNAseq
#' @importFrom QDNAseq normalizeSegmentedBins
#' @importClassesFrom QDNAseq QDNAseqReadCounts QDNAseqCopyNumbers
#' @export QDNAseqSegtoCNAclinic
#'
#' @author Dineika Chandrananda & Allan Lui
#' @examples
#'      \dontrun{
#'       vignette("CNAclinic")
#'      }
#'
QDNAseqSegtoCNAclinic=function(x,
                               callThreshLog2R=c(-0.3, 0.3)){

#############################################################################

# Check arguments & initialise dataframes

#############################################################################
  
  binsToUse <- NULL
  sampleNames <- NULL
  totalBins <- NULL
  
  segCBS <- segLACBS <- segHMM <- segPLS <- segSummary <- calls <- data.frame()
  
  if(missing(x))
    stop("x is missing")
  
  if(class(x) == "QDNAseqCopyNumbers"){
  
  sampleNames <- Biobase::sampleNames(x)
  
  copyNumber<-as.data.frame(
    QDNAseq:::log2adhoc(Biobase::assayDataElement(x, "copynumber")))
  
  names(copyNumber) <- sampleNames
  
  segCBS <- as.data.frame(
    QDNAseq:::log2adhoc(Biobase::assayDataElement(x, "segmented")),stringsAsFactors=FALSE)
  
  tempDF <- data.frame(
    chromosome=QDNAseq::chromosomes(x),
    start=QDNAseq::bpstart(x),
    end=QDNAseq::bpend(x),
    usebin=QDNAseq:::binsToUse(x),
    copyNumber,
    stringsAsFactors=FALSE)
  
  x <- tempDF
  
  rm(list=c("tempDF")); gc(FALSE)
  
  }else{
    stop(paste("\nx should either be a CNAclinicData object,",
               "a data.frame with specific columns or a QDNAseqCopyNumbers object",
               sep="\n"))
  }
  
#############################################################################

# Transform input data to correct format

#############################################################################
  
  # At this point x is a data.frame with the genomic regions & log2R values
  # with chromosome, start, end, usebin and sample columns
  
  # Check chromosome names, remove 'chr' suffix
  # convert to character 1:22, X, Y, MT and
  # set usebin=FALSE for other types of chromosomes
  
  binsToUse <- x$usebin
  tempDF <- .convertChrom(x$chromosome)
  x$chromosome <- tempDF$chromosome
  x$usebin <- binsToUse & tempDF$valid
  
  rm(list=c("binsToUse"))
  rm(list=c("tempDF")); gc(FALSE)
  
  
  binsToUse <- x$usebin
  tempDF <- .convertChrom(x$chromosome)
  x$chromosome <- tempDF$chromosome
  x$usebin <- binsToUse & tempDF$valid
  
  rm(list=c("binsToUse"))
  rm(list=c("tempDF")); gc(FALSE)
  
  x[!x$usebin, sampleNames] <- NA_real_
  
  
  # Get  totalBins
  totalBins <- nrow(x)
  
  x <- .reorderByChrom(x)
  binsToUse <- x$usebin
  
  listSegDF = list()
  listSegDF[["CBS"]] <- segCBS
  listVec <- lapply(listSegDF, c, recursive=TRUE)
  listMatrix <- do.call(cbind, listVec); rm(list=c('listVec')); gc(FALSE)
  
  segSummary <- apply(listMatrix, 1, mean, na.rm=TRUE)
  segSummary <- matrix(segSummary, ncol=length(sampleNames))
  
  bins <- data.frame(chromosome = x$chromosome,
                     start = x$start,
                     end = x$end,
                     usebin = binsToUse,
                     stringsAsFactors = FALSE)
  
  segSummary <- as.data.frame(segSummary, col.names=sampleNames)
  segSummary[!binsToUse, ] <- NA
  names(segSummary) <- sampleNames
  
  segCBS <- as.data.frame(segCBS, col.names=sampleNames)
  segCBS[!binsToUse, ] <- NA
  names(segCBS) <- sampleNames
  
  copyNumber <- as.data.frame(copyNumber, col.names=sampleNames)
  copyNumber[!binsToUse, ] <- NA
  names(copyNumber) <- sampleNames
  
  CNAdata <- CNAclinicData(
    bins=bins,
    copyNumber=copyNumber,
    segCBS=segCBS,
    segLACBS=segLACBS,
    segHMM=segHMM,
    segPLS=segPLS,
    segSummary=segSummary)
  
#############################################################################

# Call bins

#############################################################################
  CNAdata <- callData(CNAdata,
                      callTypeLog2R="summary",
                      callThreshLog2R=c(-0.3, 0.3))
  
  return(CNAdata)
  
}


#############################################################################

# Helper functions needed for chromosome formatting

#############################################################################

.reorderByChrom <- function(x){
  
  # x is a data.frame(chromosome, start, end, usebin, samp1, samp2 ...)
  # The chromosomes are assumed to be character and only inluclude X, Y, MT
  # Order chromosomes & the data from 1:22, 23, 24, 25
  # Re-convert the chromosomes back to 1:22, X, Y, MT
  
  
  chromosome <- as.character(x$chromosome)
  chromosome[which(chromosome == "X")] <- "23"
  chromosome[which(chromosome == "Y")] <- "24"
  chromosome[which(chromosome == "MT")] <- "25"
  
  x$chromosome <- as.numeric(chromosome)
  
  x <- dplyr::arrange(x,chromosome)
  
  x$chromosome <- as.character(x$chromosome)
  
  # Replace 23 by X:
  x$chromosome[which(x$chromosome == "23")] <- "X"
  
  # Replace 24 by Y
  x$chromosome[which(x$chromosome == "24")] <- "Y"
  
  # Replace 25 by MT
  x$chromosome[which(x$chromosome == "25")] <- "MT"
  
  return(x)
  
}
