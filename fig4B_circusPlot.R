library(VariantAnnotation)
library(RColorBrewer)
library(circlize)
library(EnsDb.Hsapiens.v86)
library(ComplexHeatmap)

# Read in the vcf files containing the heteroplasmic mtDNA mutations.
# The vcf files are available at the AD Knowledge Portal: https://doi.org/10.7303/syn25618990
getMtVars <- function(region="DLPFC") {
  
  if (region == "DLPFC") {
    file <- "/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz"
    vcf <- readVcf(file)
    mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$BrainRegion == "DLPFC", ]
    vcf <- vcf[, mt$SampleID]
  } else if (region == "CB") {
    file <- "/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz"
    vcf <- readVcf(file)
    mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$BrainRegion == "CB", ]
    vcf <- vcf[, mt$SampleID]
  } else if (region == "PCC") {
    file <- "/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz"
    vcf <- readVcf(file)
    mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$BrainRegion == "PCC", ]
    vcf <- vcf[, mt$SampleID]
  } else if (region == "TCX") {
    file <- "/mnt/mfs/ctcn/team/hklein/analyses/mayo_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz"
    vcf <- readVcf(file)
    mt <- read.csv("mayo.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
    vcf <- vcf[, mt$SampleID]
  } else if (region == "FP") {
    file <- "/mnt/mfs/ctcn/team/hklein/analyses/msbb_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz"
    vcf <- readVcf(file)
    mt <- read.csv("msbb.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$RaceInferred == "W", ]
    vcf <- vcf[, mt$SampleID]
  } else {
    stop("Unknown region!")
  }
  
  relAlFreq <- apply(geno(vcf)$AF, 2, function(x) {
    freq <- rep(0, length(x))
    ind <- sapply(lapply(x, is.na), all)
    freq[!ind] <- sapply(x[!ind], function (a) {
      sum(a[!is.na(a)])
    })
    return(freq)
  })
  ind <- !apply(relAlFreq == 0, 1, all)
  relAlFreq <- relAlFreq[ind, ]
  vcf <- vcf[ind, ]
  
  # all(apply(relAlFreq, 2, function(x) {return(sum(x >= 0.03 & x <= 0.9))}) == mt$NumHeteroplasmy)
  
  return(list(relAlFreq=relAlFreq, vcf=vcf, mt=mt))
}


# Figure 4B - circus plot showing heteroplasmic mtDNA variants
plotCircus <- function (maxTh=0.9, minTh=0.03) {
  
  biotypeCols <- brewer.pal(n=4, name="Pastel1")
  names(biotypeCols) <- c("tRNA", "Non-coding", "Protein coding", "rRNA")
  
  mtGenes <- genes(EnsDb.Hsapiens.v86, filter=AnnotationFilterList(SeqNameFilter("MT")))
  mtLength <- seqlengths(seqinfo((mtGenes))) # 16569
  strand(mtGenes) <- "*"
  mt <- gaps(mtGenes)
  mt <- mt[strand(mt) == "*"]
  mt$Name <- ""
  mt$Type <- "Non-coding"
  mtGenes$Type <- c(Mt_rRNA="rRNA", Mt_tRNA="tRNA", protein_coding="Protein coding")[mtGenes$gene_biotype]
  mtGenes$Name <- gsub("MT-", "", mtGenes$symbol)
  mtAnnot <- c(mt, mtGenes[, c("Name", "Type")])
  mtAnnot$Color <- biotypeCols[mtAnnot$Type]
  mtAnnot$NameG <- ""
  mtAnnot$NameG[mtAnnot$Type %in% c("rRNA", "Protein coding")] <- mtAnnot$Name[mtAnnot$Type %in% c("rRNA", "Protein coding")]
  
  # MT genome is usually plotted anti-clockwise
  rmtAnnot <- mtAnnot
  end(rmtAnnot) <- mtLength
  start(rmtAnnot) <- mtLength - end(mtAnnot) + 1
  end(rmtAnnot) <- mtLength - start(mtAnnot) + 1
  rmtAnnot <- sort(rmtAnnot)
  
  circos.clear()
  circos.par(start.degree=90, track.height=0.15, gap.degree=0, gap.after=0,
             cell.padding=c(0.02, 0, 0.02, 0))
  circos.initialize(factors="MT", xlim=c(1, mtLength))
  circos.track(ylim=c(0,1), bg.border="white", panel.fun=function(x,y) {
    circos.rect(xleft=start(rmtAnnot), xright=end(rmtAnnot), ybottom=0, ytop=1,
                col=rmtAnnot$Color)
    ind <- rmtAnnot$NameG != ""
    circos.text(x=(start(rmtAnnot)[ind] + end(rmtAnnot)[ind])/2, y=CELL_META$cell.ylim[2], # + mm_y(2),
                labels=rmtAnnot$NameG[ind], facing=NULL, niceFacing=TRUE,
                adj = c(0.5, -0.5), cex = 0.6)
  })
  
  regions <- c("DLPFC", "TCX", "FP", "CB")
  tissueCols <- c("DLPFC/ROSMAP"="royalblue", "TCX/Mayo"="salmon3", "FP/MSBB"="aquamarine4", "CB/ROSMAP"="forestgreen")
  for (i in 1:4) {
    mtVars <- getMtVars(region=regions[i])
    circos.track(ylim=c(0,1), panel.fun=function(x, y) {
      x <- mtLength - start(mtVars$vcf) + 1  # start(mtVars$vcf)
      coords <- apply(mtVars$relAlFreq, 2, function (y) {
        y <- as.numeric(y)
        ind <- y >= minTh & y <= maxTh
        if (sum(ind) > 0) {
          return(data.frame(x=x[ind], y=y[ind]))
        } else {
          return(NULL)
        }
      })
      coords <- do.call(rbind, coords)
      circos.points(x=coords$x, y=coords$y, col=tissueCols[i], cex=0.5)
    })
    rm(mtVars)
  }
  
  lgd.tracks <- Legend(labels=names(tissueCols), title="Region/Study", type="points",
                       legend_gp=gpar(col=tissueCols),
                       title_gp = gpar(fontsize = 8, fontface = "bold"),
                       labels_gp = gpar(fontsize = 8))
  draw(lgd.tracks, x = unit(1, "npc") - unit(5, "mm"), y = unit(1, "npc") - unit(5, "mm"), 
       just = c("right", "top"))
  
  lgd.annot <- Legend(labels=names(biotypeCols)[c(3, 4, 1, 2)], title="MT annotation",
                      legend_gp=gpar(fill=biotypeCols[c(3, 4, 1, 2)]),
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8))
  draw(lgd.annot, x = unit(1, "npc") - unit(6, "mm"), y = unit(5, "mm"), 
       just = c("right", "bottom"))
  
}

pdf(file="figures/figure4B.pdf", width=6, height=6)
plotCircus()
dev.off()
