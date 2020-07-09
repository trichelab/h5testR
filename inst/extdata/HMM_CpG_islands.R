library(GenomicRanges)
library(GenomeInfoDb)
#
# Reconstruct HMM-based CpG island boundaries from Hao Wu's model 
# Thanks to Dr. Wu for kindly making mm10 islands available too
#
getHMMislands <- function(refversion) { 
  model_url <- paste0("http://www.haowulab.org/software/makeCGI/",
                     "model-based-cpg-islands-", refversion, ".txt")
  message("Using ", model_url)
  islands <- makeGRangesFromDataFrame(read.table(url(model_url), head=TRUE),
                                      keep.extra.columns=TRUE)
  names(islands) <- paste0(refversion, "_HMM_CGI_", seq_along(islands))
  seqinf <- with(getChromInfoFromUCSC(refversion),
                 Seqinfo(seqnames=chrom, 
                         seqlengths=size, 
                         isCircular=circular,
                         genome=refversion))
  seqinfo(islands) <- seqinf[seqlevels(islands)]
  mcols(islands) <- mcols(islands)[, "obsExp", drop=FALSE] 
  sort(islands)
}

# human hg38
HMM_CpG_islands.hg38 <- getHMMislands("hg38")
save(HMM_CpG_islands.hg38, file="../../data/HMM_CpG_islands.hg38.rda")
     
# human hg19
HMM_CpG_islands.hg19 <- getHMMislands("hg19")
save(HMM_CpG_islands.hg19, file="../../data/HMM_CpG_islands.hg19.rda")

# mouse mm10
HMM_CpG_islands.mm10 <- getHMMislands("mm10")
save(HMM_CpG_islands.mm10, file="../../data/HMM_CpG_islands.mm10.rda")

# mouse mm9
HMM_CpG_islands.mm9 <- getHMMislands("mm9")
save(HMM_CpG_islands.mm9, file="../../data/HMM_CpG_islands.mm9.rda")


