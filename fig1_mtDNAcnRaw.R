library(ggplot2)
library(ggpubr)
library(gridExtra)
library(gtable)
library(grid)
library(cowplot)


# Fig 1A-D - Boxplot raw mtDNAcn for ROSMAP
plotMtDNAcnRosmap <- function(csv="rosmap.csv") {
  
  tissueCols <- c("DLPFC"="royalblue", "PCC"="deepskyblue", "CB"="forestgreen")
  relEffectSize <- rep(NaN, 4)
  names(relEffectSize) <- c("dlpfcQ", "dlpfcA", "pcc", "cb")
  
  mt <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  
  mtDlpfcQ <- mt[mt$BrainRegion == "DLPFC" & mt$ExtractionKit == "QIAamp", ]
  mtDlpfcA <- mt[mt$BrainRegion == "DLPFC" & mt$ExtractionKit == "AllPrepUniversal", ]
  mtPCC <- mt[mt$BrainRegion == "PCC", ]
  mtCB <- mt[mt$BrainRegion == "CB", ]
  
  # DLPFC QIAamp
  n <- table(mtDlpfcQ$pathoAD)
  labels <- c("0"=paste("Non AD\n(n=", n["0"] , ")", sep=""),
              "1"=paste("AD\n(n=", n["1"] , ")", sep=""))
  mtDlpfcQ$Dx <- factor(labels[as.character(mtDlpfcQ$pathoAD)], levels=labels)
  compare <- list(c(labels["0"], labels["1"]))
  
  pDQ <- ggplot(mtDlpfcQ, aes(x=Dx, y=mtDNAcnRaw)) +
    geom_boxplot(fill=tissueCols["DLPFC"], outlier.shape=NA) +
    geom_jitter(width=0.1) +
    stat_compare_means(comparisons=compare, method="wilcox.test") +
    xlab("") +
    ylab("mtDNAcn in DLPFC\n(QIAamp)") +
    ylim(c(floor(min(mtDlpfcQ$mtDNAcnRaw)), max(mtDlpfcQ$mtDNAcnRaw) + 600)) + 
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  relEffectSize["dlpfcQ"] <- median(mtDlpfcQ$mtDNAcnRaw[mtDlpfcQ$pathoAD == 1]) / median(mtDlpfcQ$mtDNAcnRaw[mtDlpfcQ$pathoAD == 0])
  
  
  # DLPFC AllPrepUniversal
  n <- table(mtDlpfcA$pathoAD)
  labels <- c("0"=paste("Non AD\n(n=", n["0"] , ")", sep=""),
              "1"=paste("AD\n(n=", n["1"] , ")", sep=""))
  mtDlpfcA$Dx <- factor(labels[as.character(mtDlpfcA$pathoAD)], levels=labels)
  compare <- list(c(labels["0"], labels["1"]))
  
  pDA <- ggplot(mtDlpfcA, aes(x=Dx, y=mtDNAcnRaw)) +
    geom_boxplot(fill=tissueCols["DLPFC"], outlier.shape=NA) +
    geom_jitter(width=0.1) +
    stat_compare_means(comparisons=compare, method="wilcox.test") +
    xlab("") +
    ylab("mtDNAcn in DLPFC\n(AllPrep Universal)") +
    ylim(c(floor(min(mtDlpfcA$mtDNAcnRaw)), max(mtDlpfcA$mtDNAcnRaw) + 250)) + 
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  relEffectSize["dlpfcA"] <- median(mtDlpfcA$mtDNAcnRaw[mtDlpfcA$pathoAD == 1]) / median(mtDlpfcA$mtDNAcnRaw[mtDlpfcA$pathoAD == 0])

  
  # PCC
  n <- table(mtPCC$pathoAD)
  labels <- c("0"=paste("Non AD\n(n=", n["0"] , ")", sep=""),
              "1"=paste("AD\n(n=", n["1"] , ")", sep=""))
  mtPCC$Dx <- factor(labels[as.character(mtPCC$pathoAD)], levels=labels)
  compare <- list(c(labels["0"], labels["1"]))
  
  pPCC <- ggplot(mtPCC, aes(x=Dx, y=mtDNAcnRaw)) +
    geom_boxplot(fill=tissueCols["PCC"], outlier.shape=NA) +
    geom_jitter(width=0.1) +
    stat_compare_means(comparisons=compare, method="wilcox.test") +
    xlab("") +
    ylab("mtDNAcn in PCC\n(AllPrep Universal)") +
    ylim(c(floor(min(mtPCC$mtDNAcnRaw)), max(mtPCC$mtDNAcnRaw) + 200)) + 
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  relEffectSize["pcc"] <- median(mtPCC$mtDNAcnRaw[mtPCC$pathoAD == 1]) / median(mtPCC$mtDNAcnRaw[mtPCC$pathoAD == 0])
  
  
  # CB
  n <- table(mtCB$pathoAD)
  labels <- c("0"=paste("Non AD\n(n=", n["0"] , ")", sep=""),
              "1"=paste("AD\n(n=", n["1"] , ")", sep=""))
  mtCB$Dx <- factor(labels[as.character(mtCB$pathoAD)], levels=labels)
  compare <- list(c(labels["0"], labels["1"]))
  
  pCB <- ggplot(mtCB, aes(x=Dx, y=mtDNAcnRaw)) +
    geom_boxplot(fill=tissueCols["CB"], outlier.shape=NA) +
    geom_jitter(width=0.1) +
    stat_compare_means(comparisons=compare, method="wilcox.test") +
    xlab("") +
    ylab("mtDNAcn in CB\n(Gentra Puregene)") +
    ylim(c(floor(min(mtCB$mtDNAcnRaw)), max(mtCB$mtDNAcnRaw) + 500)) + 
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  relEffectSize["cb"] <- median(mtCB$mtDNAcnRaw[mtCB$pathoAD == 1]) / median(mtCB$mtDNAcnRaw[mtCB$pathoAD == 0])
  
  return(list(dlpfcQ=pDQ, dlpfcA=pDA, pcc=pPCC, cb=pCB, relChange=relEffectSize))
}


# Fig 1E - Boxplot raw mtDNAcn for Mayo
plotMtDNAcnMayo <- function(csv="mayo.csv") {
  
  mt <- read.csv(csv, colClasses=c(SampleID="character", apoe_genotype="character"),
                   stringsAsFactors=FALSE)
  
  relEffectSize <- rep(NaN, 3)
  names(relEffectSize) <- c("AD", "PSP", "PathoAging")
  
  n <- table(mt$Diagnosis)
  labels <- c("Control"=paste("Control\n(n=", n["Control"] , ")", sep=""),
              "AD"=paste("AD\n(n=", n["AD"] , ")", sep=""),
              "PSP"=paste("PSP\n(n=", n["PSP"] , ")", sep=""),
              "PathologicAging"=paste("Path. Aging\n(n=", n["PathologicAging"] , ")", sep=""))
  mt$Dx <- factor(labels[as.character(mt$Diagnosis)], levels=labels)
  
  compare <- list(c(labels["Control"], labels["AD"]), c(labels["Control"], labels["PSP"]), c(labels["Control"], labels["PathologicAging"]))
  
  p <- ggplot(mt, aes(x=Dx, y=mtDNAcnRaw)) +
    geom_boxplot(fill="salmon", outlier.shape=NA) +
    geom_jitter(width=0.1) +
    stat_compare_means(comparisons=compare, method="wilcox.test") +
    xlab("") +
    ylab("mtDNAcn in TCX") +
    ylim(c(floor(min(mt$mtDNAcnRaw)), max(mt$mtDNAcnRaw) + 1200)) +
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  
  relEffectSize["AD"] <- median(mt$mtDNAcnRaw[mt$Diagnosis == "AD"]) / median(mt$mtDNAcnRaw[mt$Diagnosis == "Control"])
  relEffectSize["PSP"] <- median(mt$mtDNAcnRaw[mt$Diagnosis == "PSP"]) / median(mt$mtDNAcnRaw[mt$Diagnosis == "Control"])
  relEffectSize["PathoAging"] <- median(mt$mtDNAcnRaw[mt$Diagnosis == "PathologicAging"]) / median(mt$mtDNAcnRaw[mt$Diagnosis == "Control"])
  
  return(list(p=p, relChange=relEffectSize))
}


# Fig 1F - Boxplot raw mtDNAcn for MSBB
plotMtDNAcnMsbb <- function(csv="msbb.csv") {
  
  mt <- read.csv(csv, colClasses=c(SampleID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  
  relEffectSize <- c("AD"=NaN)
  
  n <- table(mt$pathoAD)
  labels <- c("0"=paste("Control\n(n=", n["0"] , ")", sep=""),
              "1"=paste("AD\n(n=", n["1"] , ")", sep=""))
  mt$Dx <- factor(labels[as.character(mt$pathoAD)], levels=labels)
  
  compare <- list(c(labels["0"], labels["1"]))
  
  p <- ggplot(mt, aes(x=Dx, y=mtDNAcnRaw)) +
    geom_boxplot(fill="aquamarine3", outlier.shape=NA) +
    geom_jitter(width=0.1) +
    stat_compare_means(comparisons=compare, method="wilcox.test") +
    xlab("") +
    ylab("mtDNAcn in FP") +
    ylim(c(floor(min(mt$mtDNAcnRaw)), max(mt$mtDNAcnRaw) + 400)) +
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  
  relEffectSize["AD"] <- median(mt$mtDNAcnRaw[mt$pathoAD == 1]) / median(mt$mtDNAcnRaw[mt$pathoAD == 0])
  
  return(list(p=p, relChange=relEffectSize))
}


# Fig S1A-E - QQ plots of normalized mtDNAcn for each region/study
generateQqPlots <- function(rosmap="rosmap.csv", mayo="mayo.csv", msbb="msbb.csv") {
  
  tissueCols <- c(DLPFC="royalblue", PCC="deepskyblue", CB="forestgreen",
                  TCX="salmon", FP="aquamarine3")
  
  
  # ROSMAP 
  mt <- read.csv(rosmap, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- split(mt, mt$BrainRegion)
  
  mtDNAcn <- data.frame(mtDNAcn=mt$DLPFC$mtDNAcn)
  stopifnot(all(is.na(mtDNAcn$mtDNAcn) == FALSE))
  ks <- ks.test(x=mtDNAcn$mtDNAcn, y="pnorm")$p.value # 0.49
  label <- paste("p=", round(ks, digits=2), ", n=", nrow(mtDNAcn), sep="")
  p1 <- ggplot(data=mtDNAcn, aes(sample=mtDNAcn)) +
    geom_qq(distribution = "qnorm") +
    geom_abline(intercept=0, slope=1, color=tissueCols["DLPFC"]) +
    annotate(geom="text", x=-1.6, y=2.7, label=label) +
    ylab("DLPFC mtDNAcn quantiles") +
    xlab("Theoretical quantiles") +
    theme_light()
  
  mtDNAcn <- data.frame(mtDNAcn=mt$PCC$mtDNAcn)
  stopifnot(all(is.na(mtDNAcn$mtDNAcn) == FALSE))
  ks <- ks.test(x=mtDNAcn$mtDNAcn, y="pnorm")$p.value # 0.89
  label <- paste("p=", round(ks, digits=2), ", n=", nrow(mtDNAcn), sep="")
  p2 <- ggplot(data=mtDNAcn, aes(sample=mtDNAcn)) +
    geom_qq(distribution = "qnorm") +
    geom_abline(intercept=0, slope=1, color=tissueCols["PCC"]) +
    annotate(geom="text", x=-1.3, y=1.65, label=label) +
    ylab("PCC mtDNAcn quantiles") +
    xlab("Theoretical quantiles") +
    theme_light()
  
  mtDNAcn <- data.frame(mtDNAcn=mt$CB$mtDNAcn)
  stopifnot(all(is.na(mtDNAcn$mtDNAcn) == FALSE))
  ks <- ks.test(x=mtDNAcn$mtDNAcn, y="pnorm")$p.value # 0.56
  label <- paste("p=", round(ks, digits=2), ", n=", nrow(mtDNAcn), sep="")
  p3 <- ggplot(data=mtDNAcn, aes(sample=mtDNAcn)) +
    geom_qq(distribution = "qnorm") +
    geom_abline(intercept=0, slope=1, color=tissueCols["CB"]) +
    annotate(geom="text", x=-1.5, y=3.05, label=label) +
    ylab("CB mtDNAcn quantiles") +
    xlab("Theoretical quantiles") +
    theme_light()
  
  
  # Mayo 
  mt <- mt <- read.csv(mayo, colClasses=c(SampleID="character", apoe_genotype="character"),
                       stringsAsFactors=FALSE)
  mtDNAcn <- data.frame(mtDNAcn=mt$mtDNAcn)
  stopifnot(all(is.na(mtDNAcn$mtDNAcn) == FALSE))
  ks <- ks.test(x=mtDNAcn$mtDNAcn, y="pnorm")$p.value  # 0.91
  label <- paste("p=", round(ks, digits=2), ", n=", nrow(mtDNAcn), sep="")
  p4 <- ggplot(data=mtDNAcn, aes(sample=mtDNAcn)) +
    geom_qq(distribution = "qnorm") +
    geom_abline(intercept=0, slope=1, color=tissueCols["TCX"]) +
    annotate(geom="text", x=-1.55, y=2.2, label=label) +
    ylab("TCX mtDNAcn quantiles") +
    xlab("Theoretical quantiles") +
    theme_light()
  
  
  # MSBB
  mt <- read.csv(msbb, colClasses=c(SampleID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mtDNAcn <- data.frame(mtDNAcn=mt$mtDNAcn)
  stopifnot(all(is.na(mtDNAcn$mtDNAcn) == FALSE))
  ks <- ks.test(x=mtDNAcn$mtDNAcn, y="pnorm")$p.value  # 0.23
  label <- paste("p=", round(ks, digits=2), ", n=", nrow(mtDNAcn), sep="")
  p5 <- ggplot(data=mtDNAcn, aes(sample=mtDNAcn)) +
    geom_qq(distribution = "qnorm") +
    geom_abline(intercept=0, slope=1, color=tissueCols["FP"]) +
    annotate(geom="text", x=-1.55, y=3.05, label=label) +
    ylab("FP mtDNAcn quantiles") +
    xlab("Theoretical quantiles") +
    theme_light()
  
  return(list(p1, p2, p3, p4, p5))
}


# Figure 1 panel
generateFigure1 <- function() {
  rosmap <- plotMtDNAcnRosmap("rosmap.csv")
  mayo <- plotMtDNAcnMayo("mayo.csv")
  msbb <- plotMtDNAcnMsbb("msbb.csv")

  effSizes <- data.frame("Region/Study"=c("DLPFC/ROSMAP\n(QIAamp)", "DLPFC/ROSMAP\n(AllPrep Universal)", "PCC/ROSMAP", "CB/ROSMAP", "TCX/Mayo", "FP/MSBB"),
                         "mtDNAcn in AD\nrelative to controls"=paste(sprintf("%3.1f", round(c(rosmap$relChange, mayo$relChange["AD"], msbb$relChange), digits=3)*100), "%", sep=""),
                         stringsAsFactors=FALSE, check.names=FALSE)

  # Generate table as last figure in panel
  table1g <- tableGrob(effSizes, rows=NULL,
                       theme=ttheme_minimal(core=list(fg_params=list(hjust=c(rep(0, 6), rep(1, 6)), x=c(rep(0.05, 6), rep(0.95, 6)))),
                                            colhead=list(fg_params=list(hjust=c(0, 1), x=c(0.05, 0.95))),
                                            padding=unit(c(5,5), "mm"), base_size=9))
  table1g <- gtable_add_grob(table1g,
                             grobs=segmentsGrob(x0=c(unit(0, "npc"), x0=unit(0, "npc")) , y0=c(unit(1, "npc"), unit(0, "npc")),
                                                x1=c(unit(1, "npc"), unit(1, "npc")), y1=c(unit(1, "npc"), unit(0, "npc")), gp=gpar(lwd=2.0)),
                             t=1, b=1, l=1, r=ncol(table1g))
  table1g <- gtable_add_grob(table1g,
                             grobs=segmentsGrob(x0=unit(0, "npc"), y0=unit(0, "npc"), x1=unit(1, "npc"), y1=unit(0, "npc"), gp=gpar(lwd=2.0)),
                             t=nrow(table1g), b=1, l=1, r=ncol(table1g))
  
  # Generate panel
  upperRow <- plot_grid(rosmap$dlpfcQ, rosmap$dlpfcA, rosmap$pcc, rosmap$cb, nrow=1,
                        labels=c("A", "B", "C", "D"))
  lowerRow <- plot_grid(mayo$p, msbb$p, table1g, nrow=1, rel_widths=c(0.45, 0.25, 0.3),
                        labels=c("E", "F", "G"))
  p <- plot_grid(upperRow, lowerRow, nrow=2)

  ggsave("figures/figure1.pdf", plot=p, width=12, height=7)
}


# Figure S1 panel
generateFigureS1 <- function() {
  
  ps <- generateQqPlots(rosmap="rosmap.csv", mayo="mayo.csv", msbb="msbb.csv")
  
  p <- plot_grid(ps[[1]], ps[[2]], ps[[3]],
                 ps[[4]], ps[[5]], NULL,
                 labels=c("A", "B", "C", "D", "E", "", nrow=2, ncol=3))
  
  ggsave("figures/figureS1.pdf", plot=p, width=12, height=7)
}


generateFigure1()
generateFigureS1()
