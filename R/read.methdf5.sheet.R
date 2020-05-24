#' process as sample sheet directly into an RGChannelSet backed by HDF5
#' 
#' For processing lists of IDATs, see read.methdf5
#' 
#' @param sheet     a sample sheet, as from read.metharray.sheet
#' @param filepath  a directory in which to put the HDF5 data ("my_h5_se")
#' @param ...       other arguments, passed to minfi:::read.metharray2
#' 
#' @return an HDF5-backed RGChannelSet
#'
#' @import utils
#' @import minfi
#' @import HDF5Array 
#' 
#' @seealso read.methdf5
#' 
#' @examples
#'
#' # cheezy 450k example
#' if (require("minfiData")) {
#'   baseDir <- system.file("extdata", package = "minfiData")
#'   sheet <- read.metharray.sheet(baseDir)
#'   inCore <- read.metharray(basenames=sheet$Basename)
#'   outOfCore <- read.methdf5.sheet(sheet)
#'   stopifnot(identical(dim(inCore), dim(outOfCore)))
#'   verifyRGsets(ram=inCore, hdf5=outOfCore)
#' }
#'
#' @export 
read.methdf5.sheet <- function(sheet, filepath="my_h5_se", ...) {

  if (!dir.exists(filepath)) dir.create(filepath, recursive=TRUE)
  RGset <- minfi:::read.metharray2(basenames=sheet$Basename, dir=filepath, ...)
  annotation(RGset) <- annotation(read.metharray(sheet[1, "Basename"]))
  return(RGset)

}
