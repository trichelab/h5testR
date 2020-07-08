#' Scan a single Tabix'ed file's ranges (only!) into a GRanges object.
#'
#' Modified from a prototype originally by Wanding Zhou (github.com/zwdzwd).
#' One useful ... argument is yieldSize, e.g. yieldSize=1e4, to avoid wedging.
#' Another useful ... argument is index, if not the usual tf.tbi scheme. 
#'
#' @param  tf         a tabix'ed file (either a character(1) or a TabixFile)
#' @param  verbose    Provide feedback on which file is being read? (TRUE) 
#' @param  is0based   Are the intervals in the tabixed file 0-based? (TRUE)
#' @param  ...        other arguments to pass to Rsamtools::TabixFile()
#'
#' @details
#'
#' Beware: the parameter 'is0based' may be extremely important for your usage. 
#' It is a good idea to look at the returned GRanges object to ensure this.
#' 
#' @return a GRanges
#'
#' @import Rsamtools
#'
#' @export
tabixToGr <- function(tf, verbose=TRUE, is0based=TRUE, ...) {

  if (!is(tf, "TabixFile")) tf <- Rsamtools::TabixFile(tf, ...)
  if (!file.exists(Rsamtools::index(tf))) stop("Your file lacks an index!")
  if (verbose) message("Reading GRanges from ", basename(path(tf)))
  
  chrs <- Rsamtools::headerTabix(tf)$seqnames
  cols <- Rsamtools::headerTabix(tf)$indexColumns
  what <- as.list(rep("character", max(cols)))
  what[cols[2:3]] <- "integer" # start, end 
  names(what) <- rep("ignore", length(what))
  names(what)[cols] <- c("seqnames", "start", "end")
 
  tfdf <- as.data.frame(scan(path(tf), what=what, sep="\t", flush=TRUE))[, cols]
  makeGRangesFromDataFrame(tfdf, starts.in.df.are.0based=is0based)

}


#' Scan a TabixFileList into a joint GRanges by union()'ing each tabixed file's
#'
#' @param tfl       a TabixFileList
#' @param verbose   squawk? (TRUE) 
#' 
#' @return a GRanges
#'
#' @import Rsamtools
#' @import GenomicRanges
#'
#' @examples
#' \dontrun{
#'   tfl <- TabixFileList(lapply(list.files(patt=".hg19.bed.gz$"), TabixFile))
#'   gr <- tabixFileListToGr(tfl) 
#' }
#'
#' @export
tabixFileListToGr <- function(tfl, verbose=TRUE) {

  if (!is(tfl, "TabixFileList")) tfl <- TabixFileList(lapply(tfl, TabixFile))
  Reduce(GenomicRanges::union, lapply(tfl, tabixToGr, verbose=verbose))

}


#' Scan Tabix'ed files into a merged GRanges[List] (for RaggedExperiments). 
#'
#' Modified from a prototype originally by Wanding Zhou (github.com/zwdzwd).
#'
#' @param ...   TABIXed files from which to generate a joint GenomicRanges
#'
#' @return a GRanges
#'
#' @import Rsamtools
#'
#' @examples
#' \dontrun{
#'   gr <- tabixedToGr(list.files(patt=".hg19.bed.gz$"))
#' }
#'
#' @export
tabixedToGr <- function(...) {

  tfs <- as.list(...) 
  if (length(as.list(...)) == 1) tfs <- tfs[[1]]
  tfl <- TabixFileList(lapply(tfs, TabixFile))
  tabixFileListToGr(tfl) 

}
