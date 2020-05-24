#' Verify that an RGChannelSet backed by HDF5 files matches an in-core twin.
#' 
#' @param ram       the in-core RGChannelSet
#' @param hdf5      the hdf5-backed RGChannelSet
#' @param rows      how many rows to test? (default: ncol(ram))
#' 
#' @return a logical of length one 
#'
#' @seealso read.methdf5
#' 
#' @export
verifyRGsets <- function(ram, hdf5, rows=NULL) {

  stopifnot(identical(ncol(ram), ncol(hdf5)))
  stopifnot(identical(nrow(ram), nrow(hdf5)))

  passed <- TRUE 
  for (ch in assayNames(hdf5)) { 
    message("Testing ", ch, " channel...", appendLF=FALSE)
    res <- verifyChannel(ch=ch, ram=ram, hdf5=hdf5, rows=rows)
    if (res == TRUE) message("OK.")
    else message("FAILED.")
    passed <- passed & res 
  } 

  return(passed)

}
