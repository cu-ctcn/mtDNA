library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)

# Figure 4A - histogram of heteroplasmic mutations in DLPFC, TCX, FP, and CB.
plotHistogramHetPlasmy <- function () {
  tissueCols <- c("DLPFC/ROSMAP"="royalblue", "CB/ROSMAP"="forestgreen",
                  "TCX/Mayo"="salmon3", "FP/MSBB"="aquamarine3")
  counts <- list()
  
  mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
  counts$DLPFC <- data.frame("Region/Study"="DLPFC/ROSMAP",
                             Heteroplasmy=mt$NumHeteroplasmy_03[mt$BrainRegion == "DLPFC"],
                             stringsAsFactors=FALSE)
  
  mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
  counts$CB <- data.frame("Region/Study"="CB/ROSMAP",
                          Heteroplasmy=mt$NumHeteroplasmy_03[mt$BrainRegion == "CB"],
                          stringsAsFactors=FALSE)
  
  mt <- read.csv("mayo.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
  counts$TCX <- data.frame("Region/Study"="TCX/Mayo",
                           Heteroplasmy=mt$NumHeteroplasmy_03,
                           stringsAsFactors=FALSE)
  
  mt <- read.csv("msbb.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
  counts$FP <- data.frame("Region/Study"="FP/MSBB",
                           Heteroplasmy=mt$NumHeteroplasmy_03[mt$RaceInferred == "W"],
                           stringsAsFactors=FALSE)
  
  # sapply(counts, function (x) {return(mean(x$Heteroplasmy))})
  #    DLPFC       CB      TCX       FP 
  # 2.821586 1.016529 2.652672 2.585185
  
  maxN <- max(counts$DLPFC$Heteroplasmy, counts$CB$Heteroplasmy, counts$TCX$Heteroplasmy, counts$FP$Heteroplasmy)
  relFreqs <- lapply(counts, function (c) {
    freq <- table(c$Heteroplasmy)
    df <- data.frame(Region.Study=rep(c$Region.Study[1], maxN+1),
                     NumHet=0:maxN,
                     freq=rep(0, maxN+1),
                     stringsAsFactors=FALSE)
    ind <- match(names(freq), as.character(df$NumHet))
    df$freq[ind] <- freq
    df$relFreq <- df$freq/sum(df$freq)
    N=sum(df$freq)
    df$Region.Study.N <- paste(df$Region.Study, "\n(n=", N, ")", sep="")
    return(df)
  })
  relFreqs <- do.call(rbind, relFreqs)
  relFreqs$Region.Study.N <- factor(relFreqs$Region.Study.N, levels=sort(unique(relFreqs$Region.Study.N))[c(2, 4, 3, 1)])
  tissueCols <- tissueCols[c(1, 3, 4, 2)]
  names(tissueCols) <- levels(relFreqs$Region.Study.N)
  
  # ggplot(counts, aes(x=Heteroplasmy, fill=Region.Study)) + geom_histogram(bins=max(counts$Heteroplasmies) + 1)
  p <- ggplot(relFreqs, aes(x=NumHet, y=relFreq, fill=Region.Study.N)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8) +
    scale_fill_manual(name="Region/Study", values=tissueCols) +
    xlab("Number of heteroplasmic mtDNA mutations") + ylab("Relative frequency") +
    theme_light() +
    theme(legend.key.height=unit(1.6, "line"), legend.position=c(0.87,0.61))
  
  return(p)
}


# Figures 4 C, D, and E - Number of heteroplasmic mutations versus age in DLPFC, TCX, and FP.
plotHetPlasmyAge <- function(tissue="DLPFC") {
  
  if (tissue %in% c("DLPFC", "PCC", "CB")) {
    mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$BrainRegion == tissue, ]
    tissueCols <- c("DLPFC"="royalblue", "PCC"="deepskyblue", "CB"="forestgreen")
    
    fit <- glm(NumHeteroplasmy_03 ~ age_death, family=quasipoisson, data=mt)
    regline <- data.frame(age_death=round(range(mt$age_death)[1]):round(range(mt$age_death)[2]))
    regline$NumHeteroplasmy <- predict(fit, newdata=regline, type="response")
    pVal <- coefficients(summary(fit))["age_death", "Pr(>|t|)"]
    coef <- coefficients(summary(fit))["age_death", "Estimate"]
    cint <- confint(fit, level=0.95)["age_death", ]
    n <- nrow(mt)
    
    p <- ggplot(mt, aes(x=age_death, y=NumHeteroplasmy_03)) +
      geom_jitter(width=0, height=0.3) +
      geom_line(data=regline, aes(x=age_death, y=NumHeteroplasmy), col=tissueCols[tissue], size=1.5) +
      xlab("Age (years)") +
      ylab(paste("Number heteroplasmic mutations in", tissue)) +
      theme_light()
    
  } else if (tissue %in% c("TCX", "FP")) {
    if (tissue == "TCX") {
      mt <- read.csv("mayo.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
    } else if (tissue == "FP") {
      mt <- read.csv("msbb.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
      mt <- mt[mt$RaceInferred == "W", ]
    }
    tissueCols <- c("TCX"="salmon3", "FP"="aquamarine4")
    
    fit <- glm(NumHeteroplasmy_03 ~ AgeStrataNumeric, family=quasipoisson, data=mt)
    regline <- data.frame(AgeStrataNumeric=round(range(mt$AgeStrataNumeric)[1]):round(range(mt$AgeStrataNumeric)[2]))
    regline$NumHeteroplasmy <- predict(fit, newdata=regline, type="response")
    pVal <- coefficients(summary(fit))["AgeStrataNumeric", "Pr(>|t|)"]
    coef <- coefficients(summary(fit))["AgeStrataNumeric", "Estimate"]
    cint <- confint(fit, level=0.95)["AgeStrataNumeric", ]
    n <- nrow(mt)
    
    xb <- c(65, 75, 85, 95)
    xl <- c("[65, 70[", "[75, 80[", "[85, 90[", "95+")
    
    p <- ggplot(mt, aes(x=AgeStrataNumeric, y=NumHeteroplasmy_03)) +
      geom_jitter(width=0.5, height=0.3) +
      geom_line(data=regline, aes(x=AgeStrataNumeric, y=NumHeteroplasmy), col=tissueCols[tissue], size=1.5) +
      scale_x_continuous(breaks=xb, labels=xl) +
      xlab("Age (years)") +
      ylab(paste("Number heteroplasmic mutations in", tissue)) +
      theme_light()
  }  
  
  return(list(plot=p, pVal=pVal, coef=coef, cint=cint, n=n))
}


# Figures S4A-C and Table S5 - Generate plots and tables for multivariable models (heteroplasmies ~ AD + age + sex + mtDNAcn)
getMultiVariableModels <- function(tissue="DLPFC") {
  
  tissueCols <- c("DLPFC/ROSMAP"="royalblue", "CB/ROSMAP"="forestgreen",
                  "TCX/Mayo"="salmon3", "FP/MSBB"="aquamarine4")
  
  if (tissue == "DLPFC") {
    mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$BrainRegion == "DLPFC", ]
    mt$pathoAD <- factor(mt$pathoAD)
    mt$msex <- factor(mt$msex)
    sdAge <- sd(mt$age_death)
    mt$age <- (mt$age_death - mean(mt$age_death)) / sdAge
    fit <- glm(NumHeteroplasmy_03 ~ pathoAD + age + msex + mtDNAcn, family=quasipoisson, data=mt)
    s <- summary(fit)
    coefs <- c("pathoAD1", "age", "msex1", "mtDNAcn")
    results <- data.frame(Tissue="DLPFC/ROSMAP",
                          Variable=c("Pathologic AD", "Age", "Male sex", "mtDNAcn"),
                          Coef=coefficients(s)[coefs, "Estimate"],
                          CIlow=confint(fit, level=0.95)[coefs, "2.5 %"],
                          CIhigh=confint(fit, level=0.95)[coefs, "97.5 %"],
                          se=coefficients(s)[coefs, "Std. Error"],
                          p=coefficients(s)[coefs, "Pr(>|t|)"],
                          n=nrow(mt),
                          stringsAsFactors=FALSE)
    results <- results[c(2,3,1,4), ]
    
    pd <- results
    pd$Color <- "black"
    pd$Color[pd$p <= 0.05] <- tissueCols[pd$Tissue[pd$p <= 0.05]]
    pd$Variable <- factor(pd$Variable, levels=pd$Variable[4:1])
    p <- ggplot(pd, aes(x=Coef, xmin=CIlow, xmax=CIhigh, y=Variable)) +
      geom_pointrange(col=pd$Color) + 
      geom_vline(xintercept=0, lty=2) +
      xlab("Effect size") +
      ylab("") +
      theme_light()
    
    return(list(plot=p, results=results, sdAge=sdAge, color=tissueCols["DLPFC/ROSMAP"]))
    
  } else if (tissue == "CB") {
    mt <- read.csv("rosmap.csv", colClasses=c(ProjID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$BrainRegion == "CB", ]
    mt$pathoAD <- factor(mt$pathoAD)
    mt$msex <- factor(mt$msex)
    sdAge <- sd(mt$age_death)
    mt$age <- (mt$age_death - mean(mt$age_death)) / sdAge
    fit <- glm(NumHeteroplasmy_03 ~ pathoAD + age + msex + mtDNAcn, family=quasipoisson, data=mt)
    s <- summary(fit)
    coefs <- c("pathoAD1", "age", "msex1", "mtDNAcn")
    results <- data.frame(Tissue="CB/ROSMAP",
                          Variable=c("Pathologic AD", "Age", "Male sex", "mtDNAcn"),
                          Coef=coefficients(s)[coefs, "Estimate"],
                          CIlow=confint(fit, level=0.95)[coefs, "2.5 %"],
                          CIhigh=confint(fit, level=0.95)[coefs, "97.5 %"],
                          se=coefficients(s)[coefs, "Std. Error"],
                          p=coefficients(s)[coefs, "Pr(>|t|)"],
                          n=nrow(mt),
                          stringsAsFactors=FALSE)
    results <- results[c(2,3,1,4), ]
    
    pd <- results
    pd$Color <- "black"
    pd$Color[pd$p <= 0.05] <- tissueCols[pd$Tissue[pd$p <= 0.05]]
    pd$Variable <- factor(pd$Variable, levels=pd$Variable[4:1])
    p <- ggplot(pd, aes(x=Coef, xmin=CIlow, xmax=CIhigh, y=Variable)) +
      geom_pointrange(col=pd$Color) + 
      geom_vline(xintercept=0, lty=2) +
      xlab("Effect size") +
      ylab("") +
      theme_light()
    
    return(list(plot=p, results=results, sdAge=sdAge, color=tissueCols["DLPFC/CB"]))
    
  } else if (tissue == "TCX") {
    mt <- read.csv("mayo.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
    mt$Diagnosis <- factor(mt$Diagnosis, levels=c("Control", "AD", "PathologicAging", "PSP"))
    mt$Sex <- factor(mt$Sex)
    sdAge <- sd(mt$AgeStrataNumeric)
    mt$age <- (mt$AgeStrataNumeric - mean(mt$AgeStrataNumeric)) / sdAge
    fit <- glm(NumHeteroplasmy_03 ~ Diagnosis + age + Sex + mtDNAcn, family=quasipoisson, data=mt)
    s <- summary(fit)
    coefs <- c("DiagnosisAD", "DiagnosisPathologicAging", "DiagnosisPSP", "age", "SexM", "mtDNAcn")
    results <- data.frame(Tissue="TCX/Mayo",
                          Variable=c("AD", "Path. aging", "PSP", "Age", "Male sex", "mtDNAcn"),
                          Coef=coefficients(s)[coefs, "Estimate"],
                          CIlow=confint(fit, level=0.95)[coefs, "2.5 %"],
                          CIhigh=confint(fit, level=0.95)[coefs, "97.5 %"],
                          se=coefficients(s)[coefs, "Std. Error"],
                          p=coefficients(s)[coefs, "Pr(>|t|)"],
                          n=nrow(mt),
                          stringsAsFactors=FALSE)
    results <- results[c(4,5,1,2,3,6), ]
    
    pd <- results
    pd$Color <- "black"
    pd$Color[pd$p <= 0.05] <- tissueCols[pd$Tissue[pd$p <= 0.05]]
    pd$Variable <- factor(pd$Variable, levels=pd$Variable[6:1])
    p <- ggplot(pd, aes(x=Coef, xmin=CIlow, xmax=CIhigh, y=Variable)) +
      geom_pointrange(col=pd$Color) + 
      geom_vline(xintercept=0, lty=2) +
      xlab("Effect size") +
      ylab("") +
      theme_light()
    
    return(list(plot=p, results=results, sdAge=sdAge, color=tissueCols["TCX/Mayo"]))
    
    
  } else if (tissue == "FP") {
    mt <- read.csv("msbb.csv", colClasses=c(SampleID="character"), stringsAsFactors=FALSE)
    mt <- mt[mt$RaceInferred == "W", ]
    mt$pathoAD <- factor(mt$pathoAD)
    mt$Sex <- factor(mt$Sex)
    sdAge <- sd(mt$AgeStrataNumeric)
    mt$age <- (mt$AgeStrataNumeric - mean(mt$AgeStrataNumeric)) / sdAge
    fit <- glm(NumHeteroplasmy_03 ~ pathoAD + age + Sex + mtDNAcn, family=quasipoisson, data=mt)
    s <- summary(fit)
    coefs <- c("pathoAD1", "age", "SexM", "mtDNAcn")
    results <- data.frame(Tissue="FP/MSBB",
                          Variable=c("AD", "Age", "Male sex", "mtDNAcn"),
                          Coef=coefficients(s)[coefs, "Estimate"],
                          CIlow=confint(fit, level=0.95)[coefs, "2.5 %"],
                          CIhigh=confint(fit, level=0.95)[coefs, "97.5 %"],
                          se=coefficients(s)[coefs, "Std. Error"],
                          p=coefficients(s)[coefs, "Pr(>|t|)"],
                          n=nrow(mt),
                          stringsAsFactors=FALSE)
    results <- results[c(2,3,1,4), ]
    
    pd <- results
    pd$Color <- "black"
    pd$Color[pd$p <= 0.05] <- tissueCols[pd$Tissue[pd$p <= 0.05]]
    pd$Variable <- factor(pd$Variable, levels=pd$Variable[6:1])
    p <- ggplot(pd, aes(x=Coef, xmin=CIlow, xmax=CIhigh, y=Variable)) +
      geom_pointrange(col=pd$Color) + 
      geom_vline(xintercept=0, lty=2) +
      xlab("Effect size") +
      ylab("") +
      theme_light()
    
    return(list(plot=p, results=results, sdAge=sdAge, color=tissueCols["FP/MSBB"]))
  }
}


# Figure 4 panel (circus plot to be added manually as panel B)
generateFigure4 <- function() {

  a <- plotHistogramHetPlasmy()
  c <- plotHetPlasmyAge(tissue="DLPFC")
  d <- plotHetPlasmyAge(tissue="TCX")
  e <- plotHetPlasmyAge(tissue="FP")
  
  c$plot <- c$plot + annotate("text", x=78, y=8, parse=TRUE,
                              label=paste("p == ", gsub("e-0", "%*%10^-", formatC(c$pVal, format="e", digits=1)), sep=""))
  d$plot <- d$plot + annotate("text", x=68, y=8, parse=TRUE,
                              label=paste("p == ", gsub("e-0", "%*%10^-", formatC(d$pVal, format="e", digits=1)), sep=""))
  e$plot <- e$plot + annotate("text", x=71.5, y=6, parse=TRUE,
                              label=paste("p == ", gsub("e-0", "%*%10^-", formatC(e$pVal, format="e", digits=1)), sep=""))
  
  # Generate table
  effSizes <- data.frame("Region/\nStudy"=c("DLPFC/\nROSMAP", "TCX/\nMayo", "FP/\nMSBB"),
                         "Increase of heteroplas-\nmic mutations per year"=c(paste(sprintf("%3.2f", round(c$coef * 100, digits=2)), "% [", sprintf("%3.2f", round(c$cint[1] * 100, digits=2)), "%; ", sprintf("%3.2f", round(c$cint[2] * 100, digits=2)), "%]", sep=""),
                                                                             paste(sprintf("%3.2f", round(d$coef * 100, digits=2)), "% [", sprintf("%3.2f", round(d$cint[1] * 100, digits=2)), "%; ", sprintf("%3.2f", round(d$cint[2] * 100, digits=2)), "%]", sep=""),
                                                                             paste(sprintf("%3.2f", round(e$coef * 100, digits=2)), "% [", sprintf("%3.2f", round(e$cint[1] * 100, digits=2)), "%; ", sprintf("%3.2f", round(e$cint[2] * 100, digits=2)), "%]", sep="")),
                         stringsAsFactors=FALSE, check.names=FALSE)
  table4 <- tableGrob(effSizes, rows=NULL,
                      theme=ttheme_minimal(core=list(fg_params=list(hjust=c(rep(0, 3), rep(1, 3)), x=c(rep(0.05, 3), rep(0.95, 3)))),
                                           colhead=list(fg_params=list(hjust=c(0, 1), x=c(0.05, 0.95))),
                                           padding=unit(c(5,5), "mm"), base_size=9))
  table4 <- gtable_add_grob(table4,
                            grobs=segmentsGrob(x0=c(unit(0, "npc"), x0=unit(0, "npc")) , y0=c(unit(1, "npc"), unit(0, "npc")),
                                                x1=c(unit(1, "npc"), unit(1, "npc")), y1=c(unit(1, "npc"), unit(0, "npc")), gp=gpar(lwd=2.0)),
                             t=1, b=1, l=1, r=ncol(table4))
  table4 <- gtable_add_grob(table4,
                            grobs=segmentsGrob(x0=unit(0, "npc"), y0=unit(0, "npc"), x1=unit(1, "npc"), y1=unit(0, "npc"), gp=gpar(lwd=2.0)),
                            t=nrow(table4), b=1, l=1, r=ncol(table4))
  
  
  leftCol <- plot_grid(a, nullGrob(), ncol=1, labels=c("A", "B"), rel_heights=c(1, 2))
  rightCol <- plot_grid(c$plot, d$plot, e$plot, table4, ncol=2, labels=c("C", "D", "E", "F"))
  p <- plot_grid(leftCol, rightCol, nrow=1, rel_widths = c(2, 1.8))
  
  ggsave("figures/figure4.pdf", plot=p, width=12, height=8)
}


# Figure S4 panel
generateFigure43 <- function() {
  
  a <- getMultiVariableModels(tissue="DLPFC")
  b <- getMultiVariableModels(tissue="TCX")
  c <- getMultiVariableModels(tissue="FP")
  
  a$plot <- a$plot + annotate("text", x=0.1, y=4.25, parse=TRUE, color=a$color,
                              label=paste("p == ", gsub("e-0", "%*%10^-", formatC(a$results["age", "p"], format="e", digits=1)), sep=""))
  b$plot <- b$plot + annotate("text", x=0.22, y=6.3, parse=TRUE, color=b$color,
                              label=paste("p == ", gsub("e-0", "%*%10^-", formatC(b$results["age", "p"], format="e", digits=1)), sep="")) +
    annotate("text", x=0.23, y=4.3, parse=TRUE, color=b$color,
             label=paste("p == ", formatC(b$results["DiagnosisAD", "p"], digits=2), sep="")) +
    annotate("text", x=0.15, y=1.3, parse=TRUE, color=b$color,
             label=paste("p == ", formatC(b$results["mtDNAcn", "p"], digits=2), sep=""))
  c$plot <- c$plot + annotate("text", x=0.16, y=4.25, parse=TRUE, color=c$color,
                              label=paste("p == ", gsub("e-0", "%*%10^-", formatC(c$results["age", "p"], format="e", digits=1)), sep=""))
  
  p <- plot_grid(a$plot, b$plot, c$plot, nrow=1, labels=c("A", "B", "C"))
  ggsave("figures/figureS4.pdf", plot=p, width=12, height=4)
}


# Table S5 - Table with details for the multivariable models
# Effect size for age is rescaled from per standard deviation age to per year (corresponding Figure S4 shows effect size in sd of age)
generateTableS5 <- function() {
  
  a <- getMultiVariableModels(tissue="DLPFC")
  b <- getMultiVariableModels(tissue="TCX")
  c <- getMultiVariableModels(tissue="FP")
  
  table <- data.frame(Variable=c("Age", "Sex", "AD", "Path. aging", "PSP", "mtDNAcn"),
                      DLPFC.b=c(a$results["age", "Coef"] / a$sdAge, a$results["msex1", "Coef"], a$results["pathoAD1", "Coef"], NA, NA, a$results["mtDNAcn", "Coef"]),
                      DLPFC.se=c(a$results["age", "se"] / a$sdAge,  a$results["msex1", "se"],   a$results["pathoAD1", "se"],   NA, NA, a$results["mtDNAcn", "se"]),
                      DLPFC.p=c(a$results["age", "p"],              a$results["msex1", "p"],    a$results["pathoAD1", "p"],    NA, NA, a$results["mtDNAcn", "p"]),
                      
                      TCX.b= c(b$results["age", "Coef"] / b$sdAge, b$results["SexM", "Coef"], b$results["DiagnosisAD", "Coef"], b$results["DiagnosisPathologicAging", "Coef"], b$results["DiagnosisPSP", "Coef"], b$results["mtDNAcn", "Coef"]),
                      TCX.se=c(b$results["age", "se"] / b$sdAge,   b$results["SexM", "se"],   b$results["DiagnosisAD", "se"],   b$results["DiagnosisPathologicAging", "se"],   b$results["DiagnosisPSP", "se"],   b$results["mtDNAcn", "se"]),
                      TCX.p= c(b$results["age", "p"],              b$results["SexM", "p"],    b$results["DiagnosisAD", "p"],    b$results["DiagnosisPathologicAging", "p"],    b$results["DiagnosisPSP", "p"],    b$results["mtDNAcn", "p"]),
                      
                      FP.b= c(c$results["age", "Coef"] / c$sdAge, c$results["SexM", "Coef"], c$results["pathoAD1", "Coef"], NA, NA, c$results["mtDNAcn", "Coef"]),
                      FP.se=c(c$results["age", "se"] / c$sdAge,   c$results["SexM", "se"],   c$results["pathoAD1", "se"],   NA, NA, c$results["mtDNAcn", "se"]),
                      FP.p= c(c$results["age", "p"],              c$results["SexM", "p"],    c$results["pathoAD1", "p"],    NA, NA, c$results["mtDNAcn", "p"]),
                      
                      stringsAsFactors=FALSE)
  
  table$DLPFC.b <- round(table$DLPFC.b, digits=3)
  table$TCX.b <- round(table$TCX.b, digits=3)
  table$FP.b <- round(table$FP.b, digits=3)
  table$DLPFC.se <- round(table$DLPFC.se, digits=3)
  table$TCX.se <- round(table$TCX.se, digits=3)
  table$FP.se <- round(table$FP.se, digits=3)
  table$DLPFC.p <- c(formatC(table$DLPFC.p[1], format="e", digits=1), formatC(table$DLPFC.p[2:6], digits=2))
  table$TCX.p <- c(formatC(table$TCX.p[1], format="e", digits=1), formatC(table$TCX.p[2:6], digits=2))
  table$FP.p <- c(formatC(table$FP.p[1], format="e", digits=1), formatC(table$FP.p[2:6], digits=2))
  
  write.csv(table, file="tables/tableS5.csv", row.names=FALSE)
}


generateFigure4()
generateFigureS4()
generateTableS5()
