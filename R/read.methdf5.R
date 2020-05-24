#' process [an] IDAT[s] directly into an RGChannelSet backed by HDF5
#' 
#' For processing entire sample sheets, see read.methdf5.sheet
#' 
#' @param basenames basename(s) (i.e., path(s) to pair(s) of IDAT files) 
#' @param filepath  a directory in which to put the HDF5 data ("my_h5_se")
#' @param ...       other arguments, passed to minfi:::read.metharray2
#' 
#' @return an HDF5-backed RGChannelSet
#'
#' @import minfi
#' @import HDF5Array 
#' 
#' @seealso read.methdf5.sheet
#' 
#' @examples
#'
#' # cheezy 450k example
#' if (require("minfiData")) {
#'   baseDir <- system.file("extdata", package = "minfiData")
#'   sheet <- read.metharray.sheet(baseDir)
#'   inCore <- read.metharray(basenames=sheet[1, "Basename", drop=FALSE])
#'   outOfCore <- read.methdf5(basenames=sheet[1, "Basename", drop=FALSE])
#' }
#'
#' # set TRUE for real test
#' RUN_TARGET_EPIC <- FALSE
#' 
#' # the real thing, with TARGET AML data
#' # yes all 500 will run on your laptop.
#' # no it will not be particularly fast.
#' # be certain that you have disk space!
#' if (RUN_TARGET_EPIC) {
#'
#'   if (!file.exists("GSE124413/GSE124413_RAW.tar")) {
#'     message("Warning: you are about to download 7GB of compressed IDATs.") 
#'     message("Exit *NOW* if you do not have the space or desire to do so.")
#'     getGEOSuppFiles("GSE124413") 
#'   } 
#'
#'   # ...time passes...
#'   if (length(list.files(patt="GSM3532678_200989060236_R08C01_")) < 2) {
#'     untar("GSE124413/GSE124413_RAW.tar")
#'   }
#' 
#' }
#' 
#' # how many arrays to test on? 10? 50? All 500? 
#' RUN_TARGET_HOWMANY <- 10 
#' 
#' if (RUN_TARGET_EPIC) { 
#'
#'   set.seed(123) 
#'   files <- list.files(patt="idat.gz")))
#'   stubs <- unique(gsub("_(Grn|Red).idat.gz", "", files))
#'   basenames <- sample(stubs, RUN_TARGET_HOWMANY)
#'
#'   outOfCore <- read.methdf5(basenames=stubs[indices])
#'   saveHDF5SummarizedExperiment(outOfCore, dir="TARGET_RGset", replace=TRUE)
#'   outOfCore <- loadHDF5SummarizedExperiment(dir="TARGET_RGset")
#' 
#'   inCore <- read.metharray(basenames=stubs[indices])
#'   stopifnot(identical(dim(inCore), dim(outOfCore)))
#' 
#'   # FIXME: refactor this as a unit test 
#'   verifyChannel <- function(ch, ram, hdf5) {
#'     idx <- seq(1, ncol(ram))
#'     funs <- c(Red=getRed, Green=getGreen)
#'     identical(as.matrix(funs[channel](hdf5[chunk, chunk])),
#'               as.matrix(funs[channel](ram[chunk, chunk])))
#'   }
#'   
#'   # FIXME: refactor this as a unit test 
#'   verifyRGsets <- function(ram, hdf5) {
#'     stopifnot(identical(ncol(ram), ncol(hdf5)))
#'     stopifnot(identical(nrow(ram), nrow(hdf5)))
#'     for (ch in assayNames(hdf5)) { 
#'       message("Testing ", ch, " channel...", appendLF=FALSE)
#'       if (verifyChannel(ch, ram, hdf5)) {
#'         message("OK.")
#'       } else { 
#'         message("FAILED.")
#'       }
#'     } 
#'   }
#'
#'   verifyRGsets(ram=inCore, hdf5=outOfCore)
#' 
#' } 
#' 
#' @export 
read.methdf5 <- function(basenames, filepath="my_h5_se", ...) {

  if (!dir.exists(filepath)) dir.create(filepath, recursive=TRUE)
  RGset <- minfi:::read.metharray2(basenames=basenames, dir=filepath, ...)
  annotation(RGset) <- annotation(read.metharray(basenames[1])) # FIXME!
  return(RGset)

}
