library(meta)
library(mediation)


# Table 2 - APOE e4 effect on mtDNAcn adjusted for pathology in all three cohorts
getApoe4Analysis <- function (round=FALSE) {

  alleleN <- list()
  
  
  # ROSMAP DLPFC
  cn <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  cn <- cn[cn$BrainRegion == "DLPFC", ]
  # sum(!is.na(cn$apoe_genotype)) # 454
  cn$apoe4 <- 0
  cn$apoe4[cn$apoe_genotype %in% c("24", "34")] <- 1
  cn$apoe4[cn$apoe_genotype %in% c("44")] <- 2
  alleleN$dlpfc <- table(cn$apoe4)
  dlpfc <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + msex + age_death, data=cn)
  dlpfc.adj <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + gpath_sqrt + msex + age_death, data=cn)
  # summary(dlpfc)
  # summary(dlpfc.adj)
  
  
  # ROSMAP CB
  cn <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  cn <- cn[cn$BrainRegion == "CB", ]
  # sum(!is.na(cn$apoe_genotype)) # 236
  cn$apoe4 <- 0
  cn$apoe4[cn$apoe_genotype %in% c("24", "34")] <- 1
  cn$apoe4[cn$apoe_genotype %in% c("44")] <- 2
  alleleN$cb <- table(cn$apoe4)
  cb <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + msex + age_death, data=cn)
  cb.adj <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + gpath_sqrt + msex + age_death, data=cn)
  # summary(cb)
  # summary(cb.adj)
  
  
  # Mayo TCX
  cn <- read.csv("mayo.csv", colClasses=c(SampleID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  # nrow(cn) # 262
  cn$apoe4 <- 0
  cn$apoe4[cn$apoe_genotype %in% c("24", "34")] <- 1
  cn$apoe4[cn$apoe_genotype %in% c("44")] <- 2
  alleleN$tcx <- table(cn$apoe4)
  tcx <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + Sex + AgeStrataNumeric, data=cn)
  tcx.adj <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + Diagnosis + Sex + AgeStrataNumeric, data=cn)
  # summary(tcx)
  # summary(tcx.adj)
  
  
  # MSBB FP
  cn <- read.csv("msbb.csv", colClasses=c(SampleID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  cn <- cn[cn$RaceInferred == "W", ]
  # nrow(cn) # 270
  cn$apoe4 <- 0
  cn$apoe4[cn$apoe_genotype %in% c("24", "34")] <- 1
  cn$apoe4[cn$apoe_genotype %in% c("44")] <- 2
  alleleN$fp <- table(cn$apoe4)
  fp <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + Sex + AgeStrataNumeric, data=cn)
  fp.adj <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + pathoAD + Sex + AgeStrataNumeric, data=cn)
  # summary(fp)
  # summary(fp.adj)

    
  apoeTable <- data.frame(Region=c("DLPFC", "CB", "TCX", "FP", "Meta"),
                          Apoe4N=c(sapply(alleleN, paste, collapse="/"), paste(colSums(do.call(rbind, alleleN)), collapse="/")),
                          c=c(coefficients(summary(dlpfc))["apoe4", "Estimate"], coefficients(summary(cb))["apoe4", "Estimate"],
                              coefficients(summary(tcx))["apoe4", "Estimate"], coefficients(summary(fp))["apoe4", "Estimate"], NA),
                          se=c(coefficients(summary(dlpfc))["apoe4", "Std. Error"], coefficients(summary(cb))["apoe4", "Std. Error"],
                               coefficients(summary(tcx))["apoe4", "Std. Error"], coefficients(summary(fp))["apoe4", "Std. Error"], NA),
                          p=c(coefficients(summary(dlpfc))["apoe4", "Pr(>|t|)"], coefficients(summary(cb))["apoe4", "Pr(>|t|)"],
                              coefficients(summary(tcx))["apoe4", "Pr(>|t|)"], coefficients(summary(fp))["apoe4", "Pr(>|t|)"], NA),
                          c.adj=c(coefficients(summary(dlpfc.adj))["apoe4", "Estimate"], coefficients(summary(cb.adj))["apoe4", "Estimate"],
                                  coefficients(summary(tcx.adj))["apoe4", "Estimate"], coefficients(summary(fp.adj))["apoe4", "Estimate"], NA),
                          se.adj=c(coefficients(summary(dlpfc.adj))["apoe4", "Std. Error"], coefficients(summary(cb.adj))["apoe4", "Std. Error"],
                                   coefficients(summary(tcx.adj))["apoe4", "Std. Error"], coefficients(summary(fp.adj))["apoe4", "Std. Error"], NA),
                          p.adj=c(coefficients(summary(dlpfc.adj))["apoe4", "Pr(>|t|)"], coefficients(summary(cb.adj))["apoe4", "Pr(>|t|)"],
                                  coefficients(summary(tcx.adj))["apoe4", "Pr(>|t|)"], coefficients(summary(fp.adj))["apoe4", "Pr(>|t|)"], NA))
  
  ma <- metagen(TE=apoeTable$c[1:4],
                seTE=apoeTable$se[1:4],
                studlab=c("DLPFC", "CB", "TCX", "FP"),
                comb.fixed=FALSE,
                comb.random=TRUE,
                method.tau="REML",
                sm="MD")
  apoeTable$c[5] <- ma$TE.random
  apoeTable$se[5] <- ma$seTE.random
  apoeTable$p[5] <- ma$pval.random
  
  ma.adj <- metagen(TE=apoeTable$c.adj[1:4],
                    seTE=apoeTable$se.adj[1:4],
                    studlab=c("DLPFC", "CB", "TCX", "FP"),
                    comb.fixed=FALSE,
                    comb.random=TRUE,
                    method.tau="REML",
                    sm="MD")
  apoeTable$c.adj[5] <- ma.adj$TE.random
  apoeTable$se.adj[5] <- ma.adj$seTE.random
  apoeTable$p.adj[5] <- ma.adj$pval.random
  
  rownames(apoeTable) <- NULL
  if (round) {
    apoeTable$c <- round(apoeTable$c, digits=3)
    apoeTable$se <- round(apoeTable$se, digits=3)
    apoeTable$p[apoeTable$p >= 0.001] <- round(apoeTable$p, digits=3)[apoeTable$p >= 0.001]
    apoeTable$p[apoeTable$p < 0.001] <- formatC(apoeTable$p, format="e", digits=1)[apoeTable$p < 0.001]
    apoeTable$c.adj <- round(apoeTable$c.adj, digits=3)
    apoeTable$se.adj <- round(apoeTable$se.adj, digits=3)
    apoeTable$p.adj[apoeTable$p.adj >= 0.001] <- round(apoeTable$p.adj, digits=3)[apoeTable$p.adj >= 0.001]
    apoeTable$p.adj[apoeTable$p.adj < 0.001] <- formatC(apoeTable$p.adj, format="e", digits=1)[apoeTable$p.adj < 0.001]
  }
  
  return(apoeTable)
}


# Figure S3 - APOE e4 - tau - mtDNAcn mediation analysis in the DLPFC
# (no figure is plotted in R, just analysis)
mediationApoe4 <- function () {
  cn <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  cn <- cn[cn$BrainRegion == "DLPFC", ]
  # sum(!is.na(cn$apoe_genotype)) # 454
  
  cn$gpath_sqrt <- (cn$gpath_sqrt - mean(cn$gpath_sqrt)) / sd(cn$gpath_sqrt)
  cn$apoe4 <- 0
  cn$apoe4[cn$apoe_genotype %in% c("24", "34")] <- 1
  cn$apoe4[cn$apoe_genotype %in% c("44")] <- 2
  
  dlpfc.gpath <- lm(gpath_sqrt ~ apoe4 + PC1 + PC2 + PC3 + msex + age_death, data=cn)
  dlpfc.mtDNAcn <- lm(mtDNAcn ~ apoe4 + PC1 + PC2 + PC3 + msex + age_death, data=cn)
  dlpfc.full <- lm(mtDNAcn ~ apoe4 + gpath_sqrt + PC1 + PC2 + PC3 + msex + age_death, data=cn)

  set.seed(15032021)
  res <- mediate(model.m=dlpfc.gpath, model.y=dlpfc.full, treat="apoe4", mediator="gpath_sqrt")
  
  # apoe4 -> pathology
  summary(dlpfc.gpath)
  # apoe4        0.755129   0.092917   8.127 4.35e-15 ***
  
  # apoe4 -> mtDNAcn
  summary(dlpfc.mtDNAcn)
  # apoe4       -0.364029   0.099397  -3.662  0.00028 ***
  
  # pathology -> mtDNAcn
  summary(dlpfc.full)
  # apoe4       -0.204296   0.104501  -1.955   0.0512 .
  # gpath_sqrt  -0.211530   0.049654  -4.260 2.49e-05 ***
  
  
  # mediation
  summary(res)
  #                Estimate 95% CI Lower 95% CI Upper  p-value    
  # ACME             -0.160       -0.255        -0.08  <2e-16 ***
  # ADE              -0.203       -0.410         0.00   0.054 .  
  # Total Effect     -0.363       -0.553        -0.17   0.002 ** 
  # Prop. Mediated    0.443        0.186         1.01   0.002 ** 
  
  return(res)
}


apoe4Table <- getApoe4Analysis(round=TRUE)
write.csv(apoe4Table, file="tables/table2.csv", row.names = FALSE, quote=FALSE)

# getApoe4Analysis(round=FALSE)
#   Region     Apoe4N          c         se            p       c.adj     se.adj      p.adj
# 1  DLPFC  335/112/7 -0.3640285 0.09939676 2.798039e-04 -0.20429633 0.10450150 0.05121193
# 2     CB   182/56/4 -0.2923348 0.13434815 3.055589e-02 -0.22310786 0.14440576 0.12369595
# 3    TCX   183/70/9 -0.2798025 0.11408014 1.485004e-02 -0.12582116 0.12136367 0.30085622
# 4     FP  171/88/11 -0.1581928 0.10617806 1.374534e-01 -0.04829677 0.10362137 0.64153912
# 5   Meta 871/326/31 -0.2749370 0.05571390 8.023236e-07 -0.14123061 0.05768216 0.01434811