#' create a GenomicFiles from the results of sesamizeToBeds (or similar) 
#' 
#' see sesamizeToBeds for an easy way to go from IDATs to tabixed BED files. 
#' 
#' @param beds    a character vector of files or data.frame with the column $bed
#' @param rr      rowRanges for GenomicFiles object (default: union BED ranges)
#' @param cd      colData for GenomicFiles object (default: skeletal info)
#' @param verbose squawk? (TRUE) 
#'
#' @return  a GenomicFiles object
#' 
#' @import rtracklayer
#' @import GenomicFiles
#' 
#' @examples
#'
#' # median methylation across a range
#' mm <- function(range, file, ...) { 
#'   res <- rtracklayer::import(file, which=range, ...)$score
#'   if (!is.null(res)) res <- median(res, na.rm=TRUE)
#'   if (is.null(res)) is.na(res) <- TRUE 
#'   return(res)
#' }
#' 
#' # EPIC data, hg38 mappings, parallel
#' if (require("minfiDataEPIC")) {
#'
#'   baseDir <- system.file("extdata", package = "minfiDataEPIC")
#'   sheet <- minfi::read.metharray.sheet(baseDir) 
#'
#'   # don't panic, feedback from sesamizeToBeds will appear eventually 
#'   beds <- sesamizeToBeds(sheet, refversion="hg38", BPPARAM=MulticoreParam(3))
#'
#'   # package up the resulting beds 
#'   gf <- genomicFilesFromBeds(beds)
#'   show(gf)
#'
#'   # get one probe's worth of data
#'   cg21870274 <- GRanges("chr1", IRanges(69591), "*", name="cg21870274")
#'   names(cg21870274) <- "cg21870274"
#'   genome(cg21870274) <- "hg38"
#'   unlist(reduceByRange(ranges=cg21870274, files=files(gf), MAP=mm)[[1]])
#'
#'   # summarize across some CpG islands
#'   data("HMM_CpG_islands.hg38")
#'   set.seed(1234)
#'   someIndices <- sample(seq_along(HMM_CpG_islands.hg38), 20)
#'   cpgi <- HMM_CpG_islands.hg38[someIndices]
#'   se <- reduceByRange(ranges=cpgi, files=files(gf), MAP=mm, summarize=TRUE)
#'   show(se)
#'   head(assays(se)$data)
#' 
#'   # information extracted from each BED
#'   colData(gf)
#' 
#' }
#' 
#' @export
genomicFilesFromBeds <- function(beds, rr=NULL, cd=NULL, verbose=TRUE) { 
 
  if (is(beds, "data.frame") | is(beds, "data.frame")) {
    stopifnot("bed" %in% names(beds))
    beds <- beds$bed
  } 

  refversion <- unique(.getGenome(beds))
  if (length(refversion) > 1) stop("Multiple genome assemblies detected!")
  
  if (is.null(rr)) {
    rr <- tabixFileListToGr(beds)
    genome(rr) <- refversion
  }
  
  n <- length(beds)
  if (is.null(cd)) {
    platform <- .getPlatform(beds)
    cd <- DataFrame(format=rep("bed", n),
                    platform=.getPlatform(beds),
                    genome=rep(refversion, n))
    rownames(cd) <- .getStubs(beds)
  }

  GenomicFiles(files=beds, rowRanges=rr, colData=cd)

} 


# helper fn
.strPop <- function(x, y=".", z=5) { 
  xx <- strsplit(x, y, fixed=TRUE)[[1]]
  if (z > length(xx)) return(c())
  else return(xx[seq(length(xx) - z + 1, length(xx))])
}


# helper fn
.strUnshift <- function(x, y=".", z=5) { 
  xx <- strsplit(x, y, fixed=TRUE)[[1]]
  if (z > length(xx)) return('')
  else return(paste(xx[seq(1, length(xx) - z)], collapse=y))
}


# helper fn 
.getPlatform <- function(beds) vapply(beds, .strPop, character(5))[1,]


# helper fn
.getGenome <- function(beds) vapply(beds, .strPop, character(5))[3,]


# helper fn
.getStubs <- function(beds) vapply(beds, .strUnshift, character(1))
