#' process IDAT pair[s] directly into tabix'ed BED files using Sesame
#' 
#' Given either a character vector or a list-like object (e.g. a data.frame), 
#' with elements $Basename and (optionally) $Sample_Name if using sample names,
#' write out beta values into a BED file per sample (see Details, below). 
#' If a data.frame is provided, the IDATs may be from multiple platforms; 
#' sesame::readIDATpair() will try to determine which platform they came from.
#'
#' @param target      an IDAT stub (e.g. "5723646052_R02C02") or a data.frame
#' @param refversion  a reference genome assembly: "hg19" (default) or "hg38"
#' @param renameBeds  if target has an element named Sample_Name, use it? (TRUE)
#' @param rschProbes  retain SNP and CpH probes in the output? (FALSE)
#' @param tabix       compress and tabix the file(s) generated? (TRUE) 
#' @param verbose     be verbose while processing? (TRUE) 
#' @param BPPARAM     a BiocParallelParam object of some sort (SerialParam())
#' @param ...         additional arguments for sesame::getBetas() 
#' 
#' @return the filename(s) of the generated BED file(s) and tabix index(es)
#'
#' @import Rsamtools
#' @import rtracklayer
#' @import GenomicFiles
#' @import BiocParallel
#' @import sesame
#' 
#' @details
#' 
#' If target is a stub (which can be a path), read the pair of IDATs for it, 
#' process typically (noob/nonlinearDyeBias/pOOBAH), mask typically, and write
#' out a BED file for the beta value (M+15/(M+15+U+15)) at each probe. 
#' 
#' The BED file(s) will be named target.platform.meth.refversion.bed.gz, e.g.
#'
#'   5723646052_R02C02.HM450.meth.hg19.bed.gz
#'
#' in the single-sample hm450 example below. This is tabixed if tabix = TRUE,
#' which by default it is, so the corresponding tabix index will be 
#' 
#'   5723646052_R02C02.HM450.meth.hg19.bed.gz.tbi
#' 
#' If target is a data.frame with columns Basename and, perhaps, Sample_Name, 
#' process a number of such files in parallel, using Basename as the stub and,
#' if column Sample_Name is present and renameBeds is true, substitute this in
#' as the prefix for the BED files. In the hm450 example below, this yields
#'
#'   GroupA_3.HM450.meth.hg19.bed.gz # and GroupA_3.HM450.meth.hg19.bed.gz.tbi
#'   ...
#'   GroupB_2.HM450.meth.hg19.bed.gz # and GroupB_2.HM450.meth.hg19.bed.gz.tbi
#' 
#' If renameBeds is FALSE, the files will be named as for single-sample runs.
#' The resulting BED files, regardless of arguments, are always of the format 
#'
#' chrom  start end probeName value * 
#' 
#' If BED files exist for the samples being processed, they may be overwritten.
#' 
#' The platform is included in the BED filename for use in annotating across
#' multiple BED files or platforms as necessary.
#' 
#' @examples
#'
#' # 450k data, hg19 mappings
#' if (require("minfiData")) {
#'   hm450BaseDir <- system.file("extdata", package = "minfiData")
#'   hm450Sheet <- minfi::read.metharray.sheet(hm450BaseDir) 
#'   hm450Files <- sesamizeToBeds(hm450Sheet[1, ]) # single-sample list
#'   unlink(hm450Files)
#'   hm450Files <- sesamizeToBeds(hm450Sheet$Basename[1]) # single-sample string
#'   unlink(hm450Files)
#'   with(hm450Sheet, sesamizeToBeds(Basename)) # multi-sample character vector
#' }
#' 
#' # EPIC data, hg38 mappings
#' if (require("minfiDataEPIC")) {
#'   epicBaseDir <- system.file("extdata", package = "minfiDataEPIC")
#'   epicSheet <- minfi::read.metharray.sheet(epicBaseDir) 
#'   epicFiles <- sesamizeToBeds(epicSheet$Basename[1], refversion="hg38")
#'   unlink(epicFiles)
#'   sesamizeToBeds(epicSheet, refversion="hg38") # multi-sample df 
#' }
#' 
#' @export 
sesamizeToBeds <- function(target, refversion=c("hg19","hg38"), renameBeds=TRUE, rschProbes=FALSE, tabix=TRUE, verbose=TRUE, BPPARAM=SerialParam(), ...) {

  # if processing multiple samples, bplapply here:
  if (is(target, "character") & length(target) > 1) {
    target <- data.frame(Basename=target)
  }
  if (.isDfOk(target)) {
    res <- as.data.frame(do.call(rbind,
                                 bplapply(.byBarcode(target), sesamizeToBeds, 
                                          ..., refversion=refversion, 
                                          renameBeds=renameBeds, 
                                          rschProbes=rschProbes,
                                          tabix=tabix, 
                                          verbose=verbose,
                                          BPPARAM=BPPARAM)))
  } else { 
   
    # processing a single sample 
    refversion <- match.arg(refversion) 
    mappings <- paste("mapped", "probes", refversion, sep=".")
    if (is(target, "character") & length(target) == 1) { 
      target <- list(Basename=target)
    }
    target <- .checkElt(target)
    stub <- basename(target$Basename)
    if ("Sample_Name" %in% names(target) & renameBeds) {
      stub <- target$Sample_Name
    }
    if (verbose) message("Reading IDATs for ", stub, "... ", appendLF=FALSE)
    sset <- dyeBiasCorrTypeINorm(noob(readIDATpair(target$Basename)))
    platform <- sset@platform
    if (verbose) message("done.")
    if (verbose) message("Processing with ", mappings, "... ", appendLF=FALSE)
    anno <- sesameDataGet(paste0(platform, ".probeInfo"))
    betas <- getBetas(sset, ...)
    retain <- names(which(!is.na(betas)))
    probes <- anno[[mappings]]
    retain <- intersect(retain, names(probes))
    probes <- probes[retain]
    probes$name <- retain
    probes$score <- betas[retain]
    if (verbose) message("done.")
    if (!rschProbes) probes <- subset(probes, grepl("^cg", probes$name))
    bedfile <- paste(stub, platform, "meth", refversion, "bed", sep=".")
    if (file.exists(bedfile)) unlink(bedfile)
    probes <- sort(probes)
    if (verbose) message("Writing betas to ", bedfile, "... ", appendLF=FALSE)
    export(probes, bedfile) 
    if (verbose) message("done.")
    gzipped <- paste0(bedfile, ".gz")
    if (verbose) message("Compressing into ", gzipped, "... ", appendLF=FALSE)
    if (file.exists(gzipped)) unlink(gzipped)
    bgzip(bedfile, dest=paste0(bedfile, ".gz"))
    if (verbose) message("done.")
    if (verbose) message("Removing ", bedfile, "... ", appendLF=FALSE)
    unlink(bedfile)
    if (verbose) message("done.")
    res <- c(bed=gzipped)
    if (tabix) {
      tabixed <- paste0(gzipped, ".tbi") 
      if (verbose) message("Indexing ", tabixed, "... ", appendLF=FALSE)
      if (file.exists(tabixed)) unlink(tabixed)
      tabixed <- indexTabix(gzipped, format="bed")
      if (verbose) message("done.")
      res <- c(bed=gzipped, tabix=tabixed)
    }
  }

  return(res)

}


# helper fn
.isDfOk <- function(x) {
  (is(x, "data.frame") | is(x, "DataFrame")) & "Basename" %in% names(x)
}


# helper fn
.byBarcode <- function(x) { 
  stopifnot(.isDfOk(x))
  columns <- intersect(c("Basename", "Sample_Name"), names(x))
  lapply(split(x[, columns], basename(x$Basename)), as.list)
}


# helper fn
.checkElt <- function(elt) {
  stopifnot("Basename" %in% names(elt))
  stopifnot(is(elt, "List") | is(elt, "list"))
  stopifnot(length(elt) %in% c(1, 2))
  items <- intersect(c("Basename", "Sample_Name"), names(elt))
  return(elt[items])
}
