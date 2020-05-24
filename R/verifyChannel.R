#' Verify that a channel of an HDF5 RGChannelSet matches its in-core twin.
#' 
#' @param ch        the channel name ("Red" or "Green") 
#' @param ram       the in-core RGChannelSet
#' @param hdf5      the hdf5-backed RGChannelSet
#' @param rows      how many rows to test? (default: ncol(ram))
#' 
#' @return a logical of length one 
#'
#' @seealso read.methdf5
#' 
#' @export
verifyChannel <- function(ch=c("Red","Green"), ram, hdf5, rows=NULL) {

  j <- seq(1, ncol(ram))
  if (is.null(rows)) {
    i <- j
  } else { 
    i <- seq(1, rows)
  }
  ch <- match.arg(ch)
  funs <- c(Red=minfi::getRed, 
            Green=minfi::getGreen)
  identical(as.matrix(funs[ch](hdf5[i, j])),
            as.matrix(funs[ch](ram[i, j])))

}
