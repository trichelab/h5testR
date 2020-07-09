#' Hidden Markov Model (HMM) based CpG island boundaries for various genomes.
#' 
#' In inst/extdata/HMM_CpG_Islands.R you can see how these were constructed. 
#' All of the island models are provided by Dr. Hao Wu of Emory University.
#' The makeCGI software provided by Dr. Wu can be used to generate new models.
#' 
#' hgXX datasets are for human genome assemblies (hg38 and hg19).
#' (Note that hg19 is equivalent to GRCh37 and hg38 is equivalent to GRCh38.)
#' 
#' mmXX datasets are for mouse genome assemblies (mm10 and mm9). 
#' (Note that mm10 is equivalent to GRCm37 and mm9 is equivalent to GRCm37.)
#'
#' The HMM island paper: \url{https://doi.org/10.1093/biostatistics/kxq005}
#' The citation: Wu et al., Redefining CpG islands using hidden Markov models,
#'               Biostatistics, Volume 11, Issue 3, July 2010, pages 499-514.
#'
#' @format A GRanges object with one metadata column (obsExp, observed/expected)
#' \describe{
#'   \item{seqnames}{the chromosome the island is on}
#'   \item{ranges}{an IRanges object with the start and end coordinates}
#'   \item{strand}{the strand the island is on (all islands are unstranded)}
#'   \item{obsExp}{observed over expected CpG density across this island}
#' }
#'
#' @aliases HMM_CpG_islands.hg38 HMM_CpG_islands.hg19
#' @aliases HMM_CpG_islands.mm10 HMM_CpG_islands.mm9
#' 
#' @seealso \url{http://www.haowulab.org/software/makeCGI/}
#' 
#' @source \url{http://www.haowulab.org/software/makeCGI/}
"HMM_CpG_islands.hg38"
