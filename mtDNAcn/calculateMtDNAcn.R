library(GenomicAlignments)
library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19.masked)

# Calculates the autosomal coverage, Mt coverage, and the mtDNAcn
calculateMtDNAcn <- function (sample) {
  
  bamfile <- paste("/mnt/mfs/ctcn/datasets/rosmap/wgs/ampad/bam/", sample, ".final.bam", sep="")
  stopifnot(file.exists(bamfile))

  gIndex <- read.table("/mnt/mfs/ctcn/resources/GRCh37/human_g1k_v37/human_g1k_v37.fasta.fai",
                       sep="\t", header=FALSE, stringsAsFactors=FALSE)
  gIndex <- gIndex[1:25, ]
  chrCovs <- data.frame(chr=gIndex$V1, chrSize=gIndex$V2, chrSizeUN=NA,
                        medianCov=NA, madCov=NA, meanCov=NA, sdCov=NA, numReads=NA)
  
  for (i in 1:nrow(chrCovs)) {
    print(i)
    param <- ScanBamParam(what="mapq", which=GRanges(seqname=chrCovs$chr[i], IRanges(1, chrCovs$chrSize[i])),
                          flag=scanBamFlag(isUnmappedQuery=FALSE))
    reads <- readGAlignments(bamfile, index=gsub(".bam", "", bamfile, fixed=TRUE), param=param)
    reads <- reads[elementMetadata(reads)$mapq != 0, ]
    cov <- coverage(reads)[[gIndex$V1[i]]]
    
    # Two default masks are activated: "assembly gaps" and "intra-contig ambiguities"
    m <- collapse(masks(BSgenome.Hsapiens.UCSC.hg19.masked[[paste("chr", gsub("MT", "M", chrCovs$chr[i]), sep="")]]))
    r <- nir_list(gaps(m))[[1]]
    
    # chr MT has no masked regions. In the BSgenome package, the MT genome has 16571 instead of 16569 bp,
    # indicating that the package comes with the older (originally used with hg19) MT version (NC_001807.4).
    # We work with the newer widely accepted rCRS/Mitomap sequence (NC_012920) and used this version at the
    # sequence alignment step.
    if (chrCovs$chr[i] == "MT") {
      end(r) <- 16569
    }
    
    chrCovs$chrSizeUN[i] <- sum(width(r))
    chrCovs$medianCov[i] <- median(cov[r])
    chrCovs$madCov[i] <- mad(cov[r])
    chrCovs$meanCov[i] <- mean(cov[r])
    chrCovs$sdCov[i] <- sd(cov[r])
    chrCovs$numReads[i] <- length(reads)
  }
  
  AutosomalCov <- median(chrCovs$medianCov[chrCovs$chr %in% as.character(1:22)])
  MtCov <- chrCovs$medianCov[chrCovs$chr %in% "MT"]
  mtDNAcn <- (chrCovs$medianCov[chrCovs$chr %in% "MT"] / AutosomalCov) * 2
  Xcn <- (chrCovs$medianCov[chrCovs$chr %in% "X"] / AutosomalCov) * 2
  Ycn <- (chrCovs$medianCov[chrCovs$chr %in% "Y"] / AutosomalCov) * 2
  
  return(c(AutosomalCov=AutosomalCov, MtCov=MtCov, mtDNAcnRaw=mtDNAcn, Xcn=Xcn, Ycn=Ycn))
}
