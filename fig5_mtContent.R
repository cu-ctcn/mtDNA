library(SummarizedExperiment)
library(ggcorrplot)
library(huge)
library(bootnet)
library(qgraph)
library(cowplot)
library(ggraph)
library(dplyr)

# Read the normalized TMT proteomics data available at Synapse:
#   data file: syn21266454; md5sum=fc68fceafbd81cb77b311f6171f92be6
#   phenotypes: syn21266449; md5sum=58eafa88bae2c14e20cc6df87590567b
readTmtData <- function() {

  if (file.exists("pSet.Rdata")) {load("pSet.Rdata"); return(pSet)}
  
  file <- "/mnt/mfs/ctcn/datasets/rosmap/tmt/dlpfcTissue/synapseRel1/output/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"
  tmt <- read.table(file, sep=",", quote = "\"", header=TRUE, row.names=1, stringsAsFactors=FALSE)
  rn <- strsplit(rownames(tmt), "|", fixed=TRUE)
  gene <- sapply(rn, "[", 1)
  protein <- sapply(rn, "[", 2)
  rownames(tmt) <- protein
    
  pheno <- read.table("/mnt/mfs/ctcn/datasets/rosmap/tmt/dlpfcTissue/synapseRel1/input/rosmap_50batch_specimen_metadata_for_batch_correction.csv",
                      sep=",", header=TRUE, stringsAsFactors=FALSE)
  pheno$projid[!is.na(pheno$projid)] <- sprintf("%08d", pheno$projid[!is.na(pheno$projid)])
  stopifnot(colnames(tmt) %in% pheno$SampleID)
  pheno <- pheno[match(colnames(tmt), pheno$SampleID), ]
  stopifnot(pheno$SampleID == colnames(tmt))
  stopifnot(!any(duplicated(pheno$projid)))
  rownames(pheno) <- pheno$projid
  colnames(tmt) <- pheno$projid
  pheno <- pheno[, c(1, 2, 4, 5, 6, 7)]

  tmt <- as.matrix(tmt)
  pSet <- SummarizedExperiment(assays=list(tmt=tmt), colData=pheno)
  stopifnot(all(rownames(pSet) == protein))
  rowData(pSet)$UniProt <- protein
  rowData(pSet)$Symbol <- gene
    
  return(pSet)
}


# Proteins used to calculate mt content score
# https://www.proteinatlas.org/humanproteome/cell/mitochondria - Table 1 "Selection of proteins suitable as markers for mitochondria."
# Manual: Outer membrane proteins specific to MT based on Human Protein Atlas.
getMtContentProteins <- function () {
  
  mtContentProt <- read.csv(text="protein, symbol, source, note
O75390,CS,HumanProteinAtlas,mitochondrial matrix
P42704,LRPPRC,HumanProteinAtlas,RNA metabolism/regulation
Q6NUK1,SLC25A24,HumanProteinAtlas,carrier protein
O43615,TIMM44,HumanProteinAtlas,inner mitochondrial membrane
Q92947,GCDH,HumanProteinAtlas,enzyme in mitochondrial matrix
Q12931,TRAP1,HumanProteinAtlas,mitochondrial chaperone protein/HSP90 family
P10809,HSPD1,Manual,mitochondrial chaperone protein/HSP60 family
P45880,VDAC2,Manual,outer mitochondrial membrane
Q9Y277,VDAC3,Manual,outer mitochondrial membrane
Q15388,TOMM20,Manual,outer mitochondrial membrane
Q9NS69,TOMM22,Manual,outer mitochondrial membrane",
                            stringsAsFactors = FALSE)
  return(mtContentProt)
}


# Calculate Mt content score
calculateMtContentScore <- function(pSet) {
  mtContentProt <- getMtContentProteins()
  tmt <- assay(pSet)[mtContentProt$protein, ]
  stopifnot(!any(is.na(tmt)))
  med <- apply(tmt, 2, median)
  return(med)
}


# Calculate Mt complex scores based on the annotation given in
# Supplementary Excel File 3. Also calculates a score for
# mtDNA-encoded proteins.
calculateMtComplexScores <- function(pSet, excel3="tables/SupExcelFile3.csv") {
  
  # Read Supplementary Exel File 3 with gene annotation
  S3 <- read.csv(file=excel3, stringsAsFactors=FALSE)
  mtGenes <- split(S3$UniprotId, f=S3$Complex)
  mtGenes$mtEncoded <- S3$UniprotId[S3$DNA == "mtDNA"]
  
#  sapply(mtGenes, length)
#  Complex I  Complex II Complex III  Complex IV   Complex V   mtEncoded
#         55           6          15          40          20          11
  
  med <- numeric()
  for (i in 1:length(mtGenes)) {
    tmt <- assay(pSet)[mtGenes[[i]], ]
    med <- cbind(med, apply(tmt, 2, median, na.rm=TRUE))
  }
  colnames(med) <- names(mtGenes)
  
  return(med)
}


# Fig S5A - correlation between proteins used for Mt content score
plotMtContentScore <- function()  {
  
  mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- mt[mt$BrainRegion == "DLPFC", ]
  pSet <- readTmtData()
  pSet <- pSet[, colnames(pSet) %in% mt$ProjID]
  mtcProts <- getMtContentProteins()
  
  tmt <- assay(pSet)[mtcProts$protein, ]
  rownames(tmt) <- mtcProts$symbol
  
  cmat <- cor(t(tmt))
  p <- ggcorrplot(cmat, type="lower", show.diag=FALSE, lab=TRUE, lab_size=2.6, tl.cex=10)

  return(p)
}


# Fig 5A - correlation matrix of the mt measures and phenotypes
# Fig S5B - samples sizes for the pairwise correlations calculated for Fig 5A
plotCormat <- function() {
  
  mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- mt[mt$BrainRegion == "DLPFC", ]
  mt$NumHeteroplasmy <- as.vector(huge.npn(x=matrix(mt$NumHeteroplasmy_03, ncol=1)))
  rownames(mt) <- mt$ProjID
  pSet <- readTmtData()
  pSet <- pSet[, colnames(pSet) %in% mt$ProjID]
  # ncol(pSet) # 156
  
  mt$MtContent <- NA
  mt[colnames(pSet), ]$MtContent <- calculateMtContentScore(pSet)
  mtComplex <- calculateMtComplexScores(pSet)
  mt$MtComplexI <- NA
  mt$MtComplexII <- NA
  mt$MtComplexIII <- NA
  mt$MtComplexIV <- NA
  mt$MtComplexV <- NA
  mt$mtDNAencoded <- NA
  mt[rownames(mtComplex), ]$MtComplexI <- mtComplex[, "Complex I"]
  mt[rownames(mtComplex), ]$MtComplexII <- mtComplex[, "Complex II"]
  mt[rownames(mtComplex), ]$MtComplexIII <- mtComplex[, "Complex III"]
  mt[rownames(mtComplex), ]$MtComplexIV <- mtComplex[, "Complex IV"]
  mt[rownames(mtComplex), ]$MtComplexV <- mtComplex[, "Complex V"]
  mt[rownames(mtComplex), ]$mtDNAencoded <- mtComplex[, "mtEncoded"]
  
  vars <- c("MtComplexI", "MtComplexII", "MtComplexIII", "MtComplexIV", "MtComplexV",
            "MtContent", "mtDNAencoded", "amyloid_sqrt", "tangles_sqrt",  "Neu",
            "cogn_global_lv", "age_death", "NumHeteroplasmy", "mtDNAcn")
  rename <- c("Mt complex I", "Mt complex II", "Mt complex III", "Mt complex IV", "Mt complex V",
              "Mt content", "mtDNA-encoded Prot.", "Amyloid", "Tau", "Proportion neurons",
              "Cognition", "Age", "mtDNA heteroplasmy", "mtDNAcn")
  mt <- mt[, vars]
  colnames(mt) <- rename
  notNa <- as.matrix(!is.na(mt))
  nmat <- t(notNa) %*% notNa # number of samples available for pairwise correlations
  cmat <- cor(mt, use="pairwise.complete")
  pmat <- cor_pmat(mt)
  #  p <- ggcorrplot(cmat, p.mat=pmat, type="lower", show.diag=FALSE, lab=TRUE, lab_size=2.6, tl.cex=10)
  p <- ggcorrplot(cmat, type="lower", show.diag=FALSE, lab=TRUE, lab_size=2.6, tl.cex=10)
  pn <- ggcorrplot(nmat, type="lower", show.diag=FALSE, lab=TRUE, lab_size=2.6, tl.cex=10, show.legend=FALSE)
  
  return(list(p=p, pn=pn))
}


# Fig 5B - sparse graph of partial correlation between phenotypes and mt variables
# Fig S5C - assessment of the stability of the graph based on bootstrapping
plotPartialCor <- function() {
  
  mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- mt[mt$BrainRegion == "DLPFC", ]
  mt$NumHeteroplasmy <- as.vector(huge.npn(x=matrix(mt$NumHeteroplasmy_03, ncol=1)))
  rownames(mt) <- mt$ProjID
  pSet <- readTmtData()
  pSet <- pSet[, colnames(pSet) %in% mt$ProjID]
  # ncol(pSet) # 156
  
  mt$MtContent <- NA
  mt[colnames(pSet), ]$MtContent <- calculateMtContentScore(pSet)
  
  # Removed "educ" and "msex" for simplicity since they don't alter
  # the network structure of the other variables 
  vars <- c("MtContent", "mtDNAcn", "NumHeteroplasmy", "amyloid_sqrt", "tangles_sqrt",  "Neu",
            "cogn_global_lv", "age_death")
  rename <- c("Mt\ncontent", "mtDNAcn", "mtDNA\nhetero-\nplasmy", "Amyloid", "Tau", "Proportion\nneurons",
              "Cognition", "Age")
  renameLong <- c("Mt content", "mtDNAcn", "mtDNA heteroplasmy", "Amyloid", "Tau", "Proportion neurons",
                  "Cognition", "Age")
  
  mt <- mt[, vars]
  colnames(mt) <- rename
  net <- estimateNetwork(data=mt, default="EBICglasso", corMethod="cov", tuning=0.5)
  colnames(mt) <- renameLong
  bn <- bootnet(data=mt, default="EBICglasso", corMethod="cov", nBoots=1000, tuning=0.5)
  
  net$graph <- net$graph * -1 # invert colors
  qnet <- plot(net)
  qnet$graphAttributes$Nodes$width <- qnet$graphAttributes$Nodes$width * 1.6
  ggnet <- as.ggraph(qnet)
  
  return(list(graph=ggnet, bootstrap=plot(bn)))
}


# Copied from qgraph 1.6.5 (was removed in later versions)
as.ggraph <- function(object){
  # To graph object using tidygraph:
  df <- as.data.frame(object$Edgelist)
  
  if (any(df$directed)){
    stop("Not yet supported for directed graphs")
  }
  
  tidyG <- tidygraph::tbl_graph(edges = df[,1:2], directed = all(df$directed))
  
  edges <- NULL
  nodes <- NULL
  order <- NULL
  n <- NULL
  
  # Add ID factor:
  tidyG <- tidyG %>% tidygraph::activate(edges) %>% dplyr::mutate(id = as.factor(seq_len(dplyr::n())))
  tidyG <- tidyG %>% tidygraph::activate(nodes) %>% dplyr::mutate(id = as.factor(seq_len(dplyr::n())), label = object$graphAttributes$Nodes$labels)
  
  # Put the edges in the right order:
  tidyG <- tidyG %>% tidygraph::activate(edges) %>% dplyr::mutate(order = order(object$graphAttributes$Graph$edgesort)) %>% dplyr::arrange(order)
  
  # Create plot:
  ggraph::ggraph(tidyG, layout = object$layout) + 
    
    # Edges:
    ggraph::geom_edge_fan(aes(
      color = id,
      edge_width = id,
      edge_linetype = id
    ), show.legend = FALSE) + 
    ggraph::scale_edge_color_manual(values = object$graphAttributes$Edges$color) +
    ggraph::scale_edge_width_manual(values = object$graphAttributes$Edges$width/2) +
    ggraph::scale_edge_linetype_manual(values = sapply(object$graphAttributes$Edges$lty,switch,
                                                       '1' = "solid", '2' = "dashed", '3'= "dotted", '4' = "dotdash", '5' = "longdash", '6' = "twodash"
    )) + 
    
    # Nodes:
    ggraph::geom_node_point(
      aes_string(
        size = "id",
        colour = "id",
        fill = "id"
      ), shape = 21, show.legend = FALSE
    ) + 
    ggplot2::scale_size_manual(values = object$graphAttributes$Nodes$width*2)  + 
    ggplot2::scale_colour_manual(values = object$graphAttributes$Nodes$border.color) +
    ggplot2::scale_fill_manual(values = object$graphAttributes$Nodes$color) +
    
    # Labels:
    ggraph::geom_node_text(
      aes_string(
        label = "label"
      )
    ) +
    
    # Limits:
    ggplot2::xlim(-1-object$plotOptions$mar[2],1+object$plotOptions$mar[3]) + 
    ggplot2::ylim(-1-object$plotOptions$mar[1],1+object$plotOptions$mar[3]) +
    
    # Theme
    ggraph::theme_graph(plot_margin=margin(0,0,0,0))
}


# Figure 5 panel
generateFigure5 <- function() {
  
  cormat <- plotCormat()
  net <- plotPartialCor()
  
  p <- plot_grid(cormat$p, net$graph, nrow=1, rel_widths=c(1.4, 1), labels=c("A", "B"))
  ggsave("figures/figure5.pdf", plot=p, width=11, height=5.5)
}


# Figure S5 panel
generateFigureS5 <- function() {
  
  mtcScore <- plotMtContentScore()
  cormat <- plotCormat()
  net <- plotPartialCor()
  
  left <- plot_grid(mtcScore, cormat$pn, nrow=2, labels=c("A", "B"))
  right <- plot_grid(net$bootstrap, nrow=1, labels="C")
  p <- plot_grid(left, net$bootstrap, nrow=1, labels=c("", "C"))
  ggsave("figures/figureS5.pdf", plot=p, width=12, height=8)
}


generateFigure5()
generateFigureS5()