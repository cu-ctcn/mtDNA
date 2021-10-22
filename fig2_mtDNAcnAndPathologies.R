library(ggplot2)
library(rsq)
library(ggpubr)
library(cowplot)


# Calculates the statistics for Figure 2A and Table S3 (univariate regression for each pathology)
# test="F" resturn results from an overall F test of each variable
# test="t" resturn results from t-tests of each variable/factor level
getUniPathoStats <- function (mtcn, test="F", standardize=TRUE) {
  
  # Standardize continuous variables for better comparability of coefficients.
  if (standardize) {
    mtcn$cogn_global_lv <- (mtcn$cogn_global_lv - mean(mtcn$cogn_global_lv, na.rm=TRUE)) / sd(mtcn$cogn_global_lv, na.rm=TRUE)
    mtcn$cogng_random_slope <- (mtcn$cogng_random_slope - mean(mtcn$cogng_random_slope, na.rm=TRUE)) / sd(mtcn$cogng_random_slope, na.rm=TRUE)
    mtcn$age_death <- (mtcn$age_death - mean(mtcn$age_death, na.rm=TRUE)) / sd(mtcn$age_death, na.rm=TRUE)
    mtcn$educ <- (mtcn$educ - mean(mtcn$educ, na.rm=TRUE)) / sd(mtcn$educ, na.rm=TRUE)
    mtcn$amyloid_sqrt <- (mtcn$amyloid_sqrt - mean(mtcn$amyloid_sqrt, na.rm=TRUE)) / sd(mtcn$amyloid_sqrt, na.rm=TRUE)
    mtcn$tangles_sqrt <- (mtcn$tangles_sqrt - mean(mtcn$tangles_sqrt, na.rm=TRUE)) / sd(mtcn$tangles_sqrt, na.rm=TRUE)
    mtcn$gpath_sqrt <- (mtcn$gpath_sqrt - mean(mtcn$gpath_sqrt, na.rm=TRUE)) / sd(mtcn$gpath_sqrt, na.rm=TRUE)
    mtcn$Neu <- (mtcn$Neu - mean(mtcn$Neu, na.rm=TRUE)) / sd(mtcn$Neu, na.rm=TRUE)
    mtcn$Ast <- (mtcn$Ast - mean(mtcn$Ast, na.rm=TRUE)) / sd(mtcn$Ast, na.rm=TRUE)
    mtcn$Oli <- (mtcn$Oli - mean(mtcn$Oli, na.rm=TRUE)) / sd(mtcn$Oli, na.rm=TRUE)
    mtcn$Mic <- (mtcn$Mic - mean(mtcn$Mic, na.rm=TRUE)) / sd(mtcn$Mic, na.rm=TRUE)
    mtcn$End <- (mtcn$End - mean(mtcn$End, na.rm=TRUE)) / sd(mtcn$End, na.rm=TRUE)
  }
  
  resAnova <- data.frame()
  resMl <- data.frame()
  
  # patho AD, adjusted for age and gender
  mtcn$pathoAD <- factor(mtcn$pathoAD)
  ind <- !is.na(mtcn$pathoAD)
  fit <- lm(mtDNAcn ~ pathoAD + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="pathoAD",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="pathoAD",
                   region=mtcn$BrainRegion[1],
                   baseline="pathoAD0",
                   baselineN=table(mtcn[ind, "pathoAD"])["0"],
                   level="pathoAD1",
                   levelN=table(mtcn[ind, "pathoAD"])["1"],
                   coef=coefficients(summary(fit))["pathoAD1", "Estimate"],
                   p=coefficients(summary(fit))["pathoAD1", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # global patho, adjusted for age and gender
  ind <- !is.na(mtcn$gpath_sqrt)
  fit <- lm(mtDNAcn ~ gpath_sqrt + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="gpath_sqrt",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="gpath_sqrt",
                   region=mtcn$BrainRegion[1],
                   baseline=NA,
                   baselineN=NA,
                   level=NA,
                   levelN=NA,
                   coef=coefficients(summary(fit))["gpath_sqrt", "Estimate"],
                   p=coefficients(summary(fit))["gpath_sqrt", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # amyloid_sqrt, adjusted for age and gender
  ind <- !is.na(mtcn$amyloid_sqrt)
  fit <- lm(mtDNAcn ~ amyloid_sqrt + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  a <- anova(fit, fitr)
  rA <- data.frame(var="amyloid_sqrt",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="amyloid_sqrt",
                   region=mtcn$BrainRegion[1],
                   baseline=NA,
                   baselineN=NA,
                   level=NA,
                   levelN=NA,
                   coef=coefficients(summary(fit))["amyloid_sqrt", "Estimate"],
                   p=coefficients(summary(fit))["amyloid_sqrt", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # tangles_sqrt, adjusted for age and gender
  ind <- !is.na(mtcn$tangles_sqrt)
  fit <- lm(mtDNAcn ~ tangles_sqrt + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  a <- anova(fit, fitr)
  rA <- data.frame(var="tangles_sqrt",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="tangles_sqrt",
                   region=mtcn$BrainRegion[1],
                   baseline=NA,
                   baselineN=NA,
                   level=NA,
                   levelN=NA,
                   coef=coefficients(summary(fit))["tangles_sqrt", "Estimate"],
                   p=coefficients(summary(fit))["tangles_sqrt", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # TDP-43, adjusted for age and gender
  mtcn$tdp_stage4 <- factor(mtcn$tdp_stage4)
  ind <- !is.na(mtcn$tdp_stage4)
  fit <- lm(mtDNAcn ~ tdp_stage4 + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="tdp_stage4",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  levels <- c("tdp_stage41", "tdp_stage42", "tdp_stage43")
  rM <- data.frame(var=rep("tdp_stage4", 3),
                   region=rep(mtcn$BrainRegion[1], 3),
                   baseline=rep("tdp_stage40", 3),
                   baselineN=rep(table(mtcn[ind, "tdp_stage4"])["0"], 3),
                   level=levels,
                   levelN=as.vector(table(mtcn[ind, "tdp_stage4"])[c("1", "2", "3")]),
                   coef=coefficients(summary(fit))[levels, "Estimate"],
                   p=coefficients(summary(fit))[levels, "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # Lewy bodies, adjusted for age and gender
  mtcn$dlbdx <- factor(mtcn$dlbdx)
  ind <- !is.na(mtcn$dlbdx)
  fit <- lm(mtDNAcn ~ dlbdx + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="dlbdx",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  levels <- c("dlbdx1", "dlbdx2", "dlbdx3")
  rM <- data.frame(var=rep("dlbdx", 3),
                   region=rep(mtcn$BrainRegion[1], 3),
                   baseline=rep("dlbdx0", 3),
                   baselineN=rep(table(mtcn[ind, "dlbdx"])["0"], 3),
                   level=levels,
                   levelN=as.vector(table(mtcn[ind, "dlbdx"])[c("1", "2", "3")]),
                   coef=coefficients(summary(fit))[levels, "Estimate"],
                   p=coefficients(summary(fit))[levels, "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # CAA, adjusted for age and gender
  mtcn$caa_4gp <- factor(mtcn$caa_4gp)
  ind <- !is.na(mtcn$caa_4gp)
  fit <- lm(mtDNAcn ~ caa_4gp + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="caa_4gp",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  levels <- c("caa_4gp1", "caa_4gp2", "caa_4gp3")
  rM <- data.frame(var=rep("caa_4gp", 3),
                   region=rep(mtcn$BrainRegion[1], 3),
                   baseline=rep("caa_4gp0", 3),
                   baselineN=rep(table(mtcn[ind, "caa_4gp"])["0"], 3),
                   level=levels,
                   levelN=as.vector(table(mtcn[ind, "caa_4gp"])[c("1", "2", "3")]),
                   coef=coefficients(summary(fit))[levels, "Estimate"],
                   p=coefficients(summary(fit))[levels, "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # Cerebral atherosclerosis, adjusted for age and gender
  mtcn$cvda_4gp2 <- factor(mtcn$cvda_4gp2)
  ind <- !is.na(mtcn$cvda_4gp2)
  fit <- lm(mtDNAcn ~ cvda_4gp2 + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="cvda_4gp2",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  levels <- c("cvda_4gp21", "cvda_4gp22", "cvda_4gp23")
  rM <- data.frame(var=rep("cvda_4gp2", 3),
                   region=rep(mtcn$BrainRegion[1], 3),
                   baseline=rep("cvda_4gp20", 3),
                   baselineN=rep(table(mtcn[ind, "cvda_4gp2"])["0"], 3),
                   level=levels,
                   levelN=as.vector(table(mtcn[ind, "cvda_4gp2"])[c("1", "2", "3")]),
                   coef=coefficients(summary(fit))[levels, "Estimate"],
                   p=coefficients(summary(fit))[levels, "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # Arteriolosclerosis: Linear model, adjusted for age and gender
  mtcn$arteriol_scler <- factor(mtcn$arteriol_scler)
  ind <- !is.na(mtcn$arteriol_scler)
  fit <- lm(mtDNAcn ~ arteriol_scler + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="arteriol_scler",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  levels <- c("arteriol_scler1", "arteriol_scler2", "arteriol_scler3")
  rM <- data.frame(var=rep("arteriol_scler", 3),
                   region=rep(mtcn$BrainRegion[1], 3),
                   baseline=rep("arteriol_scler0", 3),
                   baselineN=rep(table(mtcn[ind, "arteriol_scler"])["0"], 3),
                   level=levels,
                   levelN=as.vector(table(mtcn[ind, "arteriol_scler"])[c("1", "2", "3")]),
                   coef=coefficients(summary(fit))[levels, "Estimate"],
                   p=coefficients(summary(fit))[levels, "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # Gross chronic infarcts, adjusted for age and gender
  mtcn$ci_num2_gct <- factor(mtcn$ci_num2_gct)
  ind <- !is.na(mtcn$ci_num2_gct)
  fit <- lm(mtDNAcn ~ ci_num2_gct + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="ci_num2_gct",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="ci_num2_gct",
                   region=mtcn$BrainRegion[1],
                   baseline="ci_num2_gct0",
                   baselineN=table(mtcn[ind, "ci_num2_gct"])["0"],
                   level="ci_num2_gct1",
                   levelN=table(mtcn[ind, "ci_num2_gct"])["1"],
                   coef=coefficients(summary(fit))["ci_num2_gct1", "Estimate"],
                   p=coefficients(summary(fit))["ci_num2_gct1", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # Chronic microinfarcts, adjusted for age and gender
  mtcn$ci_num2_mct <- factor(mtcn$ci_num2_mct)
  ind <- !is.na(mtcn$ci_num2_mct)
  fit <- lm(mtDNAcn ~ ci_num2_mct + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="ci_num2_mct",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="ci_num2_mct",
                   region=mtcn$BrainRegion[1],
                   baseline="ci_num2_mct0",
                   baselineN=table(mtcn[ind, "ci_num2_mct"])["0"],
                   level="ci_num2_mct1",
                   levelN=table(mtcn[ind, "ci_num2_mct"])["1"],
                   coef=coefficients(summary(fit))["ci_num2_mct1", "Estimate"],
                   p=coefficients(summary(fit))["ci_num2_mct1", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # Hippocampal sclerosis, adjusted for age and gender
  mtcn$hspath_typ <- factor(mtcn$hspath_typ)
  ind <- !is.na(mtcn$hspath_typ)
  fit <- lm(mtDNAcn ~ hspath_typ + age_death + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="hspath_typ",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="hspath_typ",
                   region=mtcn$BrainRegion[1],
                   baseline="hspath_typ0",
                   baselineN=table(mtcn[ind, "hspath_typ"])["0"],
                   level="hspath_typ1",
                   levelN=table(mtcn[ind, "hspath_typ"])["1"],
                   coef=coefficients(summary(fit))["hspath_typ1", "Estimate"],
                   p=coefficients(summary(fit))["hspath_typ1", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # cognition, adjusted for age, gender, and education
  ind <- !is.na(mtcn$cogn_global_lv)
  fit <- lm(mtDNAcn ~ cogn_global_lv + age_death + msex + educ, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex + educ, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="cogn_global_lv",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="cogn_global_lv",
                   region=mtcn$BrainRegion[1],
                   baseline=NA,
                   baselineN=NA,
                   level=NA,
                   levelN=NA,
                   coef=coefficients(summary(fit))["cogn_global_lv", "Estimate"],
                   p=coefficients(summary(fit))["cogn_global_lv", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # cogn decline, adjusted for age and gender
  ind <- !is.na(mtcn$cogng_random_slope)
  fit <- lm(mtDNAcn ~ cogng_random_slope + age_death + msex + educ, data=mtcn)
  fitr <- lm(mtDNAcn ~ age_death + msex + educ, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="cogng_random_slope",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="cogng_random_slope",
                   region=mtcn$BrainRegion[1],
                   baseline=NA,
                   baselineN=NA,
                   level=NA,
                   levelN=NA,
                   coef=coefficients(summary(fit))["cogng_random_slope", "Estimate"],
                   p=coefficients(summary(fit))["cogng_random_slope", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # age_death, adjusted for gpath and gender
  ind <- !is.na(mtcn$age_death)
  fit <- lm(mtDNAcn ~ age_death + gpath_sqrt + msex, data=mtcn)
  fitr <- lm(mtDNAcn ~ gpath_sqrt + msex, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="age_death",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="age_death",
                   region=mtcn$BrainRegion[1],
                   baseline=NA,
                   baselineN=NA,
                   level=NA,
                   levelN=NA,
                   coef=coefficients(summary(fit))["age_death", "Estimate"],
                   p=coefficients(summary(fit))["age_death", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  # msex, adjusted for gpath and gender
  mtcn$msex <- factor(mtcn$msex)
  ind <- !is.na(mtcn$msex)
  fit <- lm(mtDNAcn ~ msex + gpath_sqrt + age_death, data=mtcn)
  fitr <- lm(mtDNAcn ~ gpath_sqrt + age_death, data=mtcn[ind, ])
  a <- anova(fit, fitr)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitr, adj = FALSE, type="v")
  rA <- data.frame(var="msex",
                   region=mtcn$BrainRegion[1],
                   rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                   p=a$`Pr(>F)`[2],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  rM <- data.frame(var="msex",
                   region=mtcn$BrainRegion[1],
                   baseline="msex0",
                   baselineN=table(mtcn[ind, "msex"])["0"],
                   level="msex1",
                   levelN=table(mtcn[ind, "msex"])["1"],
                   coef=coefficients(summary(fit))["msex1", "Estimate"],
                   p=coefficients(summary(fit))["msex1", "Pr(>|t|)"],
                   n=sum(ind),
                   stringsAsFactors=FALSE)
  resAnova <- rbind(resAnova, rA)
  resMl <- rbind(resMl, rM)
  
  if (test == "F") {
    rownames(resAnova) <- NULL
    return(resAnova)
  } else {
    rownames(resMl) <- NULL
    return(resMl)
  }
}


# Calculates the statistics for Figure 2C (mtDNAcn ~ cell type proportion)
getCellTypeStats <- function (mtcn, standardize=TRUE) {
  
  # Standardize continuous variables for better comparability of coefficients.
  if (standardize) {
    mtcn$cogn_global_lv <- (mtcn$cogn_global_lv - mean(mtcn$cogn_global_lv, na.rm=TRUE)) / sd(mtcn$cogn_global_lv, na.rm=TRUE)
    mtcn$cogng_random_slope <- (mtcn$cogng_random_slope - mean(mtcn$cogng_random_slope, na.rm=TRUE)) / sd(mtcn$cogng_random_slope, na.rm=TRUE)
    mtcn$age_death <- (mtcn$age_death - mean(mtcn$age_death, na.rm=TRUE)) / sd(mtcn$age_death, na.rm=TRUE)
    mtcn$educ <- (mtcn$educ - mean(mtcn$educ, na.rm=TRUE)) / sd(mtcn$educ, na.rm=TRUE)
    mtcn$amyloid_sqrt <- (mtcn$amyloid_sqrt - mean(mtcn$amyloid_sqrt, na.rm=TRUE)) / sd(mtcn$amyloid_sqrt, na.rm=TRUE)
    mtcn$tangles_sqrt <- (mtcn$tangles_sqrt - mean(mtcn$tangles_sqrt, na.rm=TRUE)) / sd(mtcn$tangles_sqrt, na.rm=TRUE)
    mtcn$gpath_sqrt <- (mtcn$gpath_sqrt - mean(mtcn$gpath_sqrt, na.rm=TRUE)) / sd(mtcn$gpath_sqrt, na.rm=TRUE)
    mtcn$Neu <- (mtcn$Neu - mean(mtcn$Neu, na.rm=TRUE)) / sd(mtcn$Neu, na.rm=TRUE)
    mtcn$Ast <- (mtcn$Ast - mean(mtcn$Ast, na.rm=TRUE)) / sd(mtcn$Ast, na.rm=TRUE)
    mtcn$Oli <- (mtcn$Oli - mean(mtcn$Oli, na.rm=TRUE)) / sd(mtcn$Oli, na.rm=TRUE)
    mtcn$Mic <- (mtcn$Mic - mean(mtcn$Mic, na.rm=TRUE)) / sd(mtcn$Mic, na.rm=TRUE)
    mtcn$End <- (mtcn$End - mean(mtcn$End, na.rm=TRUE)) / sd(mtcn$End, na.rm=TRUE)
  }
  
  results <- data.frame()
  mtcn <- mtcn[!is.na(mtcn$Neu), ]
  
  stopifnot(all(is.na(mtcn$Neu) == FALSE))
  fit <- lm(mtDNAcn ~ Neu + age_death + msex, data=mtcn)
  fitR <- lm(mtDNAcn ~ age_death + msex, data=mtcn)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitR, adj = FALSE, type="v")
  r <- data.frame(var="Neurons",
                  tissue=mtcn$BrainRegion[1],
                  c=coefficients(summary(fit))["Neu", "Estimate"],
                  cf.low=confint(fit, "Neu", level=0.95)[1],
                  cf.high=confint(fit, "Neu", level=0.95)[2],
                  p=coefficients(summary(fit))["Neu", "Pr(>|t|)"],
                  rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                  n=sum(!is.na(mtcn$Neu)),
                  stringsAsFactors=FALSE)
  results <- rbind(results, r)
  
  stopifnot(all(is.na(mtcn$Ast) == FALSE))
  fit <- lm(mtDNAcn ~ Ast + age_death + msex, data=mtcn)
  fitR <- lm(mtDNAcn ~ age_death + msex, data=mtcn)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitR, adj = FALSE, type="v")
  r <- data.frame(var="Astrocytes",
                  tissue=mtcn$BrainRegion[1],
                  c=coefficients(summary(fit))["Ast", "Estimate"],
                  cf.low=confint(fit, "Ast", level=0.95)[1],
                  cf.high=confint(fit, "Ast", level=0.95)[2],
                  p=coefficients(summary(fit))["Ast", "Pr(>|t|)"],
                  rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                  n=sum(!is.na(mtcn$Neu)),
                  stringsAsFactors=FALSE)
  results <- rbind(results, r)
  
  stopifnot(all(is.na(mtcn$Oli) == FALSE))
  fit <- lm(mtDNAcn ~ Oli + age_death + msex, data=mtcn)
  fitR <- lm(mtDNAcn ~ age_death + msex, data=mtcn)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitR, adj = FALSE, type="v")
  r <- data.frame(var="Oligodendrocytes",
                  tissue=mtcn$BrainRegion[1],
                  c=coefficients(summary(fit))["Oli", "Estimate"],
                  cf.low=confint(fit, "Oli", level=0.95)[1],
                  cf.high=confint(fit, "Oli", level=0.95)[2],
                  p=coefficients(summary(fit))["Oli", "Pr(>|t|)"],
                  rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                  n=sum(!is.na(mtcn$Mic)),
                  stringsAsFactors=FALSE)
  results <- rbind(results, r)
  
  stopifnot(all(is.na(mtcn$Mic) == FALSE))
  fit <- lm(mtDNAcn ~ Mic + age_death + msex, data=mtcn)
  fitR <- lm(mtDNAcn ~ age_death + msex, data=mtcn)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitR, adj = FALSE, type="v")
  r <- data.frame(var="Microglia",
                  tissue=mtcn$BrainRegion[1],
                  c=coefficients(summary(fit))["Mic", "Estimate"],
                  cf.low=confint(fit, "Mic", level=0.95)[1],
                  cf.high=confint(fit, "Mic", level=0.95)[2],
                  p=coefficients(summary(fit))["Mic", "Pr(>|t|)"],
                  rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                  n=sum(!is.na(mtcn$Mic)),
                  stringsAsFactors=FALSE)
  results <- rbind(results, r)
  
  
  stopifnot(all(is.na(mtcn$End) == FALSE))
  fit <- lm(mtDNAcn ~ End + age_death + msex, data=mtcn)
  fitR <- lm(mtDNAcn ~ age_death + msex, data=mtcn)
  rsqF <- rsq(fit, adj = FALSE, type="v")
  rsqR <- rsq(fitR, adj = FALSE, type="v")
  r <- data.frame(var="Endothelial cells",
                  tissue=mtcn$BrainRegion[1],
                  c=coefficients(summary(fit))["End", "Estimate"],
                  cf.low=confint(fit, "End", level=0.95)[1],
                  cf.high=confint(fit, "End", level=0.95)[2],
                  p=coefficients(summary(fit))["End", "Pr(>|t|)"],
                  rsq=1 - ((1 - rsqF)/(1 - rsqR)),
                  n=sum(!is.na(mtcn$End)),
                  stringsAsFactors=FALSE)
  results <- rbind(results, r)
  
  #rm("mtcn", envir=.GlobalEnv)
  return(results)
}


# Fig 2A - Variance explained by each pathology for each brain region in ROSMAP cohort
plotVarExplByPathologies <- function (csv="rosmap.csv") {
  
  tissueCols <- c("DLPFC"="royalblue", "PCC"="deepskyblue", "CB"="forestgreen")
  
  mt <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- split(mt, mt$BrainRegion)
  res <- getUniPathoStats(mt$DLPFC, test="F", standardize=TRUE)
  res <- rbind(res, getUniPathoStats(mt$PCC, test="F", standardize=TRUE))
  res <- rbind(res, getUniPathoStats(mt$CB, test="F", standardize=TRUE))
  
  res$Region <- factor(res$region, levels=c("DLPFC", "PCC", "CB"))
  
  labels <- c(pathoAD=paste("Pathologic AD\nn=(", 
                            res$n[res$var == "pathoAD" & res$Region =="DLPFC"], ", ",
                            res$n[res$var == "pathoAD" & res$Region == "PCC"], ", ",
                            res$n[res$var == "pathoAD" & res$Region == "CB"], ")", sep=""),
              gpath_sqrt=paste("Global AD path. score\nn=(",
                               res$n[res$var == "gpath_sqrt" & res$Region =="DLPFC"], ", ",
                               res$n[res$var == "gpath_sqrt" & res$Region == "PCC"], ", ",
                               res$n[res$var == "gpath_sqrt" & res$Region == "CB"], ")", sep=""),
              amyloid_sqrt=paste("Amyloid\nn=(",
                                 res$n[res$var == "amyloid_sqrt" & res$Region =="DLPFC"], ", ",
                                 res$n[res$var == "amyloid_sqrt" & res$Region == "PCC"], ", ",
                                 res$n[res$var == "amyloid_sqrt" & res$Region == "CB"], ")", sep=""),
              tangles_sqrt=paste("Tau\nn=(",
                                 res$n[res$var == "tangles_sqrt" & res$Region =="DLPFC"], ", ",
                                 res$n[res$var == "tangles_sqrt" & res$Region == "PCC"], ", ",
                                 res$n[res$var == "tangles_sqrt" & res$Region == "CB"], ")", sep=""),
              dlbdx=paste("Lewy bodies\nn=(",
                          res$n[res$var == "dlbdx" & res$Region =="DLPFC"], ", ",
                          res$n[res$var == "dlbdx" & res$Region == "PCC"], ", ",
                          res$n[res$var == "dlbdx" & res$Region == "CB"], ")", sep=""),
              tdp_stage4=paste("TDP-43\nn=(",
                               res$n[res$var == "tdp_stage4" & res$Region =="DLPFC"], ", ",
                               res$n[res$var == "tdp_stage4" & res$Region == "PCC"], ", ",
                               res$n[res$var == "tdp_stage4" & res$Region == "CB"], ")", sep=""),
              ci_num2_gct=paste("Gross chronic infarcts\nn=(",
                                res$n[res$var == "ci_num2_gct" & res$Region =="DLPFC"], ", ",
                                res$n[res$var == "ci_num2_gct" & res$Region == "PCC"], ", ",
                                res$n[res$var == "ci_num2_gct" & res$Region == "CB"], ")", sep=""),
              ci_num2_mct=paste("Chronic microinfarcts\nn=(",
                                res$n[res$var == "ci_num2_mct" & res$Region =="DLPFC"], ", ",
                                res$n[res$var == "ci_num2_mct" & res$Region == "PCC"], ", ",
                                res$n[res$var == "ci_num2_mct" & res$Region == "CB"], ")", sep=""),
              cvda_4gp2=paste("Cerebral atherosclerosis\nn=(",
                              res$n[res$var == "cvda_4gp2" & res$Region =="DLPFC"], ", ",
                              res$n[res$var == "cvda_4gp2" & res$Region == "PCC"], ", ",
                              res$n[res$var == "cvda_4gp2" & res$Region == "CB"], ")", sep=""),
              arteriol_scler=paste("Arteriolosclerosis\nn=(",
                                   res$n[res$var == "arteriol_scler" & res$Region =="DLPFC"], ", ",
                                   res$n[res$var == "arteriol_scler" & res$Region == "PCC"], ", ",
                                   res$n[res$var == "arteriol_scler" & res$Region == "CB"], ")", sep=""),
              caa_4gp=paste("CAA\nn=(",
                            res$n[res$var == "caa_4gp" & res$Region =="DLPFC"], ", ",
                            res$n[res$var == "caa_4gp" & res$Region == "PCC"], ", ",
                            res$n[res$var == "caa_4gp" & res$Region == "CB"], ")", sep=""),
              hspath_typ=paste("Hippocampal sclerosis\nn=(",
                               res$n[res$var == "hspath_typ" & res$Region =="DLPFC"], ", ",
                               res$n[res$var == "hspath_typ" & res$Region == "PCC"], ", ",
                               res$n[res$var == "hspath_typ" & res$Region == "CB"], ")", sep=""),
              cogn_global_lv=paste("Cognition\nn=(",
                                   res$n[res$var == "cogn_global_lv" & res$Region =="DLPFC"], ", ",
                                   res$n[res$var == "cogn_global_lv" & res$Region == "PCC"], ", ",
                                   res$n[res$var == "cogn_global_lv" & res$Region == "CB"], ")", sep=""),
              cogng_random_slope=paste("Cognitive decline\nn=(",
                                       res$n[res$var == "cogng_random_slope" & res$Region =="DLPFC"], ", ",
                                       res$n[res$var == "cogng_random_slope" & res$Region == "PCC"], ", ",
                                       res$n[res$var == "cogng_random_slope" & res$Region == "CB"], ")", sep=""),
              msex=paste("Sex\nn=(",
                         res$n[res$var == "msex" & res$Region =="DLPFC"], ", ",
                         res$n[res$var == "msex" & res$Region == "PCC"], ", ",
                         res$n[res$var == "msex" & res$Region == "CB"], ")", sep=""),
              age_death=paste("Age\nn=(",
                              res$n[res$var == "age_death" & res$Region =="DLPFC"], ", ",
                              res$n[res$var == "age_death" & res$Region == "PCC"], ", ",
                              res$n[res$var == "age_death" & res$Region == "CB"], ")", sep=""))
  
  res <- res[res$var != "educ", ]
  res$Variable <- factor(labels[res$var], levels=labels)
  
  res$Significance <- ""
  res$Significance[res$p <= 0.1 ] <- "⋅"
  res$Significance[res$p <= 0.05 ] <- "*"
  res$Significance[res$p <= 0.01 ] <- "**"
  res$Significance[res$p <= 0.001 ] <- "***"
  
  p <- ggplot(data=res, aes(x=Variable, y=rsq, fill=Region)) +
    geom_bar(stat="identity", width=0.8, position=position_dodge(width=0.8)) +
    geom_text(aes(y=rsq, label=Significance, group=Region), position=position_dodge(width=0.8), vjust=-0.1, size=3) +
    scale_fill_manual(values=tissueCols) +
    xlab("") + ylab(expression(paste("mtDNAcn partial R"^2))) +
    theme_light() + theme(axis.text.x=element_text(angle=35, vjust=1, hjust=1), legend.position=c(0.001, 0.996), legend.justification = c(0,1))
  
  return(p)
}


# Fig 2B and Fig S2A - boxplot showing mtDNAcn by TDP-43 stages in the PCC (Fig 2B) or DLPFC (Fig S2A)
plotTdpDetails <- function (csv="rosmap.csv", region="PCC") {
  
  tissueCols <- c("DLPFC"="royalblue", "PCC"="deepskyblue", "CB"="forestgreen")
  stopifnot(region %in% names(tissueCols))
  
  mt <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  data <- split(mt, mt$BrainRegion)[[region]]
  data <- data[!is.na(data$tdp_stage4), ]
  
  labels <- c("0"=paste("None\n(n=", sum(data$tdp_stage4 == "0"), ")", sep=""),
              "1"=paste("Amygdala\n(n=", sum(data$tdp_stage4 == "1"), ")", sep=""),
              "2"=paste("Limbic\n(n=", sum(data$tdp_stage4 == "2"), ")", sep=""),
              "3"=paste("Neocortical\n(n=", sum(data$tdp_stage4 == "3"), ")", sep=""))
  compare <- list(c(labels[1], labels[2]), c(labels[1], labels[3]), c(labels[1], labels[4]))
  
  data$tdp <- factor(labels[as.character(data$tdp_stage4)], levels=labels)
  
  p <- ggplot(data, aes(x=tdp, y=mtDNAcn)) +
    geom_boxplot(fill=tissueCols[region], outlier.shape=NA) +
    geom_jitter(width=0.1, size=1) +
    stat_compare_means(comparisons=compare, method="t.test") +
    xlab("TDP-43 pathology") +
    ylab(paste("mtDNAcn in", region)) +
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  
  return(p)
}


# Fig 2C - Variance of mtDNAcn explained by cell type proportions
plotVarExplByCellTypes <- function (csv="rosmap.csv") {
  
  mt <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- split(mt, mt$BrainRegion)
  
  res <- getCellTypeStats(mt$DLPFC)
  res$Variable <- res$var # paste(res$var, "\nn=", res$n, sep="") # 327 samples have Neu/RNA and mtDNAcn
  res$Variable <- factor(res$Variable, levels=res$Variable)
  res$Significance <- ""
  res$Significance[res$p <= 0.1 ] <- "⋅"
  res$Significance[res$p <= 0.05 ] <- "*"
  res$Significance[res$p <= 0.01 ] <- "**"
  res$Significance[res$p <= 0.001 ] <- "***"
  
  p <- ggplot(data=res, aes(x=Variable, y=rsq)) +
    geom_bar(stat="identity", width=0.8, position=position_dodge(width=0.8), fill="royalblue") +
    geom_text(aes(y=rsq, label=Significance), position=position_dodge(width=0.8), vjust=-0.1, size=3) +
    xlab("") + ylab(expression(paste("mtDNAcn partial R"^2~"in DLPFC"))) +
    theme_light() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  
  return(p)
}


# Fig 2D and Fig S2B - Forest plot showing coefficients from multivariable regression with mtDNAcn as outcome
# (Fig 2D) incl. neuronal proportion or (Fig S2B) excluding neuronal proportion
plotMultiVarForest <- function (csv="rosmap.csv", Neu=FALSE, standardize=TRUE) {
  
  mt <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  rosmap <- split(mt, mt$BrainRegion)[["DLPFC"]]
  
  # Standardize continuous variables for better comparability of coefficients.
  if (standardize) {
    rosmap$cogn_global_lv <- (rosmap$cogn_global_lv - mean(rosmap$cogn_global_lv, na.rm=TRUE)) / sd(rosmap$cogn_global_lv, na.rm=TRUE)
    rosmap$cogng_random_slope <- (rosmap$cogng_random_slope - mean(rosmap$cogng_random_slope, na.rm=TRUE)) / sd(rosmap$cogng_random_slope, na.rm=TRUE)
    rosmap$age_death <- (rosmap$age_death - mean(rosmap$age_death, na.rm=TRUE)) / sd(rosmap$age_death, na.rm=TRUE)
    rosmap$educ <- (rosmap$educ - mean(rosmap$educ, na.rm=TRUE)) / sd(rosmap$educ, na.rm=TRUE)
    rosmap$amyloid_sqrt <- (rosmap$amyloid_sqrt - mean(rosmap$amyloid_sqrt, na.rm=TRUE)) / sd(rosmap$amyloid_sqrt, na.rm=TRUE)
    rosmap$tangles_sqrt <- (rosmap$tangles_sqrt - mean(rosmap$tangles_sqrt, na.rm=TRUE)) / sd(rosmap$tangles_sqrt, na.rm=TRUE)
    rosmap$gpath_sqrt <- (rosmap$gpath_sqrt - mean(rosmap$gpath_sqrt, na.rm=TRUE)) / sd(rosmap$gpath_sqrt, na.rm=TRUE)
    rosmap$Neu <- (rosmap$Neu - mean(rosmap$Neu, na.rm=TRUE)) / sd(rosmap$Neu, na.rm=TRUE)
    rosmap$Ast <- (rosmap$Ast - mean(rosmap$Ast, na.rm=TRUE)) / sd(rosmap$Ast, na.rm=TRUE)
    rosmap$Oli <- (rosmap$Oli - mean(rosmap$Oli, na.rm=TRUE)) / sd(rosmap$Oli, na.rm=TRUE)
    rosmap$Mic <- (rosmap$Mic - mean(rosmap$Mic, na.rm=TRUE)) / sd(rosmap$Mic, na.rm=TRUE)
    rosmap$End <- (rosmap$End - mean(rosmap$End, na.rm=TRUE)) / sd(rosmap$End, na.rm=TRUE)
  }
  
  # Simplify factors with 4 levels to 2 levels to reduce number of model parameters.
  # The simpler model with binary factors has a larger BIC.
  rosmap$tdp_simple <- factor(rosmap$tdp_stage4 == 2 | rosmap$tdp_stage4 == 3)
  rosmap$dlb_simple <- factor(rosmap$dlbdx == 2 | rosmap$dlbdx == 3)
  rosmap$cvda_simple <- factor(rosmap$cvda_4gp2 == 2 | rosmap$cvda_4gp2 == 3)
  rosmap$caa_simple <- factor(rosmap$caa_4gp == 2 | rosmap$caa_4gp == 3)
  rosmap$arteriol_scler_simple <- factor(rosmap$arteriol_scler == 2 | rosmap$arteriol_scler == 3)
  
  rosmap$msex <- factor(rosmap$msex)
  rosmap$ci_num2_gct <- factor(rosmap$ci_num2_gct)
  rosmap$ci_num2_mct <- factor(rosmap$ci_num2_mct)
  rosmap$hspath_typ <- factor(rosmap$hspath_typ)
  
  
  ind <- !is.na(rosmap$mtDNAcn) & !is.na(rosmap$age_death) & !is.na(rosmap$msex) & !is.na(rosmap$amyloid_sqrt) & !is.na(rosmap$tangles_sqrt) &
    !is.na(rosmap$tdp_simple) & !is.na(rosmap$dlb_simple) & !is.na(rosmap$ci_num2_gct) & !is.na(rosmap$ci_num2_mct) & !is.na(rosmap$cvda_simple) &
    !is.na(rosmap$arteriol_scler_simple) & !is.na(rosmap$caa_simple) & !is.na(rosmap$hspath_typ)
  
  if (Neu) {
    ind <- ind &  !is.na(rosmap$Neu)
    rosmap <- rosmap[ind, ]
    model <- mtDNAcn ~ age_death + msex + amyloid_sqrt + tangles_sqrt + tdp_simple + dlb_simple + ci_num2_gct + ci_num2_mct +
                       cvda_simple + arteriol_scler_simple + caa_simple + hspath_typ + Neu
  } else {
    rosmap <- rosmap[ind, ]
    model <- mtDNAcn ~ age_death + msex + amyloid_sqrt + tangles_sqrt + tdp_simple + dlb_simple + ci_num2_gct + ci_num2_mct +
                       cvda_simple + arteriol_scler_simple + caa_simple + hspath_typ
  }
  
  # Note: There shoudn't be anymore NAs in any variable
  n <- c(age_death=sum(!is.na(rosmap$age_death)),
         msex1=sum(rosmap$msex == "1"),
         amyloid_sqrt=sum(!is.na(rosmap$amyloid_sqrt)),
         tangles_sqrt=sum(!is.na(rosmap$tangles_sqrt)),
         tdp_simpleTRUE=sum(rosmap$tdp_simple == TRUE),
         dlb_simpleTRUE=sum(rosmap$dlb_simple == TRUE),
         ci_num2_gct1=sum(rosmap$ci_num2_gct == "1"),
         ci_num2_mct1=sum(rosmap$ci_num2_mct == "1"), 
         cvda_simpleTRUE=sum(rosmap$cvda_simple == TRUE),
         arteriol_scler_simpleTRUE=sum(rosmap$arteriol_scler_simple == TRUE),
         caa_simpleTRUE=sum(rosmap$caa_simple == TRUE),
         hspath_typ1=sum(rosmap$hspath_typ == "1"),
         Neu=sum(!is.na(rosmap$Neu)))
  
  labels <- c(age_death="Age",
              educ="Education",
              msex1=paste("Sex\n(male, n=", n["msex1"], ")", sep=""),
              amyloid_sqrt="Amyloid",
              tangles_sqrt="Tau",
              dlb_simpleTRUE=paste("Lewy bodies\n(limbic or neocortical, n=", n["dlb_simpleTRUE"], ")", sep=""),
              tdp_simpleTRUE=paste("TDP-43\n(limbic or neocortical, n=", n["tdp_simpleTRUE"], ")", sep=""),
              ci_num2_gct1=paste("Gross chronic infarcts\n(one or more, n=", n["ci_num2_gct1"], ")", sep=""),
              ci_num2_mct1=paste("Chronic microinfarcts\n(one or more, n=", n["ci_num2_mct1"], ")", sep=""),
              cvda_simpleTRUE=paste("Cerebral atherosclerosis\n(moderate or severe, n=", n["cvda_simpleTRUE"], ")", sep=""),
              arteriol_scler_simpleTRUE=paste("Arteriolosclerosis\n(moderate or severe, n=", n["arteriol_scler_simpleTRUE"], ")", sep=""),
              caa_simpleTRUE=paste("Cerebral amyloid angiopathy\n(moderate or severe, n=", n["caa_simpleTRUE"], ")", sep=""),
              hspath_typ1=paste("Hippocampal sclerosis\n(present, n=", n["hspath_typ1"], ")", sep=""),
              mtDNAcn="mtDNAcn",
              Neu="Proportion of neurons")
  labels <- labels[length(labels):1]
  
  fit <- lm(model, rosmap)
  data <- as.data.frame(coefficients(summary(fit)))
  stopifnot(all(rownames(data) == row.names(confint(fit))))
  data <- cbind(data, as.data.frame(confint(fit, level=0.95)))
  
  colnames(data)[5:6] <- c("lower", "upper")
  data <- data[-1, ]
  data$Label <- factor(labels[rownames(data)], levels=labels)
  data$significant <- data$`Pr(>|t|)` <= 0.05
  data$n <- n[rownames(data)]
  
  p <- ggplot(data, aes(x=Estimate, xmin=lower, xmax=upper, y=Label, col=significant)) +
    geom_pointrange(show.legend=FALSE) + 
    geom_vline(xintercept=0, lty=2) +
    scale_color_manual(values=c("FALSE"="black", "TRUE"="royalblue")) +
    xlab("Effect size (standard deviations of mtDNAcn)") +
    ylab("") +
    theme_light()
  
  return(list(p=p, data=data))
}


# Table S3 - details from the univariate regression models (just reformating of getUniPathoStats() output)
generateUniPathoTable <- function(csv="rosmap.csv") {

  mt <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                 stringsAsFactors=FALSE)
  mt <- split(mt, mt$BrainRegion)

  tables <- list()
  for (region in names(mt)) {
    f <- getUniPathoStats(mt[[region]], test="F", standardize=TRUE)
    t <- getUniPathoStats(mt[[region]], test="t", standardize=TRUE)
    
    # Patho AD
    ind.f <- f$var == "pathoAD"
    ind.t <- t$var == "pathoAD"
    df <- data.frame(Var=c("Pathologic AD", "No AD", "AD"),
                     n=c(f$n[ind.f], t$baselineN[ind.t], t$levelN[ind.t]),
                     r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-"),
                     p_F=c(formatC(f$p[ind.f], digits=2), "-", "-"),
                     b=c("-", "-", formatC(t$coef[ind.t], format="f", digits=2)),
                     p_t=c("-", "-", formatC(t$p[ind.t], digits=2)),
                     stringsAsFactors=FALSE)
    
    # Global pathology
    ind.f <- f$var == "gpath_sqrt"
    ind.t <- t$var == "gpath_sqrt"
    df <- rbind(df, data.frame(Var="Global path. score",
                               n=f$n[ind.f],
                               r2_part=formatC(f$rsq[ind.f], format="f", digits=3),
                               p_F=formatC(f$p[ind.f], digits=2),
                               b=formatC(t$coef[ind.t], format="f", digits=2),
                               p_t=formatC(t$p[ind.t], digits=2)),
                               stringsAsFactors=FALSE)
    
    # Amyloid
    ind.f <- f$var == "amyloid_sqrt"
    ind.t <- t$var == "amyloid_sqrt"
    df <- rbind(df, data.frame(Var="Amyloid",
                               n=f$n[ind.f],
                               r2_part=formatC(f$rsq[ind.f], format="f", digits=3),
                               p_F=formatC(f$p[ind.f], digits=2),
                               b=formatC(t$coef[ind.t], format="f", digits=2),
                               p_t=formatC(t$p[ind.t], digits=2)),
                               stringsAsFactors=FALSE)
    
    # Tau
    ind.f <- f$var == "tangles_sqrt"
    ind.t <- t$var == "tangles_sqrt"
    df <- rbind(df, data.frame(Var="Tau",
                               n=f$n[ind.f],
                               r2_part=formatC(f$rsq[ind.f], format="f", digits=3),
                               p_F=formatC(f$p[ind.f], digits=2),
                               b=formatC(t$coef[ind.t], format="f", digits=2),
                               p_t=formatC(t$p[ind.t], digits=2)),
                               stringsAsFactors=FALSE)
    
    # TDP-43
    ind.f <- f$var == "tdp_stage4"
    ind.t.1 <- !is.na(t$level) & t$level == "tdp_stage41"
    ind.t.2 <- !is.na(t$level) & t$level == "tdp_stage42"
    ind.t.3 <- !is.na(t$level) & t$level == "tdp_stage43"
    df <- rbind(df, data.frame(Var=c("TDP-43", "none", "amygdala", "limbic", "neocortical"),
                               n=c(f$n[ind.f], t$baselineN[ind.t.1], t$levelN[ind.t.1], t$levelN[ind.t.2], t$levelN[ind.t.3]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-", "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-", "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t.1], format="f", digits=2), formatC(t$coef[ind.t.2], format="f", digits=2), formatC(t$coef[ind.t.3], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t.1], digits=2), formatC(t$p[ind.t.2], digits=2), formatC(t$p[ind.t.3], digits=2))),
                               stringsAsFactors=FALSE)
    
    # Lewy bodies
    ind.f <- f$var == "dlbdx"
    ind.t.1 <- !is.na(t$level) & t$level == "dlbdx1"
    ind.t.2 <- !is.na(t$level) & t$level == "dlbdx2"
    ind.t.3 <- !is.na(t$level) & t$level == "dlbdx3"
    df <- rbind(df, data.frame(Var=c("Lewy bodies", "none", "nigral", "limbic", "neocortical"),
                               n=c(f$n[ind.f], t$baselineN[ind.t.1], t$levelN[ind.t.1], t$levelN[ind.t.2], t$levelN[ind.t.3]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-", "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-", "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t.1], format="f", digits=2), formatC(t$coef[ind.t.2], format="f", digits=2), formatC(t$coef[ind.t.3], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t.1], digits=2), formatC(t$p[ind.t.2], digits=2), formatC(t$p[ind.t.3], digits=2))),
                               stringsAsFactors=FALSE)
    
    # Cerebral amyloid angiopathy
    ind.f <- f$var == "caa_4gp"
    ind.t.1 <- !is.na(t$level) & t$level == "caa_4gp1"
    ind.t.2 <- !is.na(t$level) & t$level == "caa_4gp2"
    ind.t.3 <- !is.na(t$level) & t$level == "caa_4gp3"
    df <- rbind(df, data.frame(Var=c("Cerebral amyloid angiopathy", "none", "mild", "moderate", "severe"),
                               n=c(f$n[ind.f], t$baselineN[ind.t.1], t$levelN[ind.t.1], t$levelN[ind.t.2], t$levelN[ind.t.3]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-", "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-", "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t.1], format="f", digits=2), formatC(t$coef[ind.t.2], format="f", digits=2), formatC(t$coef[ind.t.3], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t.1], digits=2), formatC(t$p[ind.t.2], digits=2), formatC(t$p[ind.t.3], digits=2))),
                               stringsAsFactors=FALSE)
    
    # Cerebral atherosclerosis
    ind.f <- f$var == "cvda_4gp2"
    ind.t.1 <- !is.na(t$level) & t$level == "cvda_4gp21"
    ind.t.2 <- !is.na(t$level) & t$level == "cvda_4gp22"
    ind.t.3 <- !is.na(t$level) & t$level == "cvda_4gp23"
    df <- rbind(df, data.frame(Var=c("Cerebral atherosclerosis", "none", "mild", "moderate", "severe"),
                               n=c(f$n[ind.f], t$baselineN[ind.t.1], t$levelN[ind.t.1], t$levelN[ind.t.2], t$levelN[ind.t.3]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-", "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-", "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t.1], format="f", digits=2), formatC(t$coef[ind.t.2], format="f", digits=2), formatC(t$coef[ind.t.3], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t.1], digits=2), formatC(t$p[ind.t.2], digits=2), formatC(t$p[ind.t.3], digits=2))),
                               stringsAsFactors=FALSE)
   
    # Arteriolosclerosis
    ind.f <- f$var == "arteriol_scler"
    ind.t.1 <- !is.na(t$level) & t$level == "arteriol_scler1"
    ind.t.2 <- !is.na(t$level) & t$level == "arteriol_scler2"
    ind.t.3 <- !is.na(t$level) & t$level == "arteriol_scler3"
    df <- rbind(df, data.frame(Var=c("Arteriolosclerosis", "none", "mild", "moderate", "severe"),
                               n=c(f$n[ind.f], t$baselineN[ind.t.1], t$levelN[ind.t.1], t$levelN[ind.t.2], t$levelN[ind.t.3]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-", "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-", "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t.1], format="f", digits=2), formatC(t$coef[ind.t.2], format="f", digits=2), formatC(t$coef[ind.t.3], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t.1], digits=2), formatC(t$p[ind.t.2], digits=2), formatC(t$p[ind.t.3], digits=2))),
                               stringsAsFactors=FALSE)
     
    # Gross chronic infarcts
    ind.f <- f$var == "ci_num2_gct"
    ind.t <- t$var == "ci_num2_gct"
    df <- rbind(df, data.frame(Var=c("Gross chronic infarcts", "None", "One or more"),
                               n=c(f$n[ind.f], t$baselineN[ind.t], t$levelN[ind.t]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t], digits=2))),
                               stringsAsFactors=FALSE)
    
    # Chronic microinfarcts
    ind.f <- f$var == "ci_num2_mct"
    ind.t <- t$var == "ci_num2_mct"
    df <- rbind(df, data.frame(Var=c("Chronic microinfarcts", "None", "One or more"),
                               n=c(f$n[ind.f], t$baselineN[ind.t], t$levelN[ind.t]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t], digits=2))),
                               stringsAsFactors=FALSE)
    
    # Hippocampal sclerosis
    ind.f <- f$var == "hspath_typ"
    ind.t <- t$var == "hspath_typ"
    df <- rbind(df, data.frame(Var=c("Hippocampal sclerosis", "Not present", "Present"),
                               n=c(f$n[ind.f], t$baselineN[ind.t], t$levelN[ind.t]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t], digits=2))),
                               stringsAsFactors=FALSE)
    
    # Global cognitive score
    ind.f <- f$var == "cogn_global_lv"
    ind.t <- t$var == "cogn_global_lv"
    df <- rbind(df, data.frame(Var="Global cognition",
                               n=f$n[ind.f],
                               r2_part=formatC(f$rsq[ind.f], format="f", digits=3),
                               p_F=formatC(f$p[ind.f], digits=2),
                               b=formatC(t$coef[ind.t], format="f", digits=2),
                               p_t=formatC(t$p[ind.t], digits=2)),
                               stringsAsFactors=FALSE)
    
    # Slope of cognitive decline
    ind.f <- f$var == "cogng_random_slope"
    ind.t <- t$var == "cogng_random_slope"
    df <- rbind(df, data.frame(Var="Cognitive decline",
                               n=f$n[ind.f],
                               r2_part=formatC(f$rsq[ind.f], format="f", digits=3),
                               p_F=formatC(f$p[ind.f], digits=2),
                               b=formatC(t$coef[ind.t], format="f", digits=2),
                               p_t=formatC(t$p[ind.t], digits=2)),
                               stringsAsFactors=FALSE)
    
    # Age
    ind.f <- f$var == "age_death"
    ind.t <- t$var == "age_death"
    df <- rbind(df, data.frame(Var="Age",
                               n=f$n[ind.f],
                               r2_part=formatC(f$rsq[ind.f], format="f", digits=3),
                               p_F=formatC(f$p[ind.f], digits=2),
                               b=formatC(t$coef[ind.t], format="f", digits=2),
                               p_t=formatC(t$p[ind.t], digits=2)),
                               stringsAsFactors=FALSE)
    
    # Sex
    ind.f <- f$var == "msex"
    ind.t <- t$var == "msex"
    df <- rbind(df, data.frame(Var=c("Sex", "Female", "Male"),
                               n=c(f$n[ind.f], t$baselineN[ind.t], t$levelN[ind.t]),
                               r2_part=c(formatC(f$rsq[ind.f], format="f", digits=3), "-", "-"),
                               p_F=c(formatC(f$p[ind.f], digits=2), "-", "-"),
                               b=c("-", "-", formatC(t$coef[ind.t], format="f", digits=2)),
                               p_t=c("-", "-", formatC(t$p[ind.t], digits=2))),
                               stringsAsFactors=FALSE)
    
    colnames(df)[-1] <- paste(region, colnames(df)[-1], sep=".")
    tables[[region]] <- df
  }
  
  stopifnot(all(tables$DLPFC$Var == tables$PCC$Var))
  stopifnot(all(tables$DLPFC$Var == tables$CB$Var))
  supTable <- cbind(cbind(tables$DLPFC, tables$PCC[, -1]), tables$CB[, -1])
  
  return(supTable)
}


# Fig 2 panel
generateFigure2 <- function() {
  pa <- plotVarExplByPathologies(csv="rosmap.csv")
  pb <- plotTdpDetails(csv="rosmap.csv", region="PCC") + ylim(c(-2.6, 3.3))
  pc <- plotVarExplByCellTypes(csv="rosmap.csv")
  pdfull <- plotMultiVarForest(csv="rosmap.csv", Neu=TRUE, standardize=TRUE)
  
  pNeu <- paste("p == ", gsub("e-05", "%*%10^-5", formatC(pdfull$data["Neu", "Pr(>|t|)"], digits=2)), sep="")
  pTau <- paste("p == ", formatC(pdfull$data["tangles_sqrt", "Pr(>|t|)"], digits=2), sep="")
  pd <- pdfull$p + annotate("text", x=0.28, y=1.5, label=pNeu, parse=TRUE, col="royalblue") +
          annotate("text", x=-0.29, y=10.4, label=pTau, parse=TRUE, col="royalblue")

  p <- plot_grid(plot_grid(pa, labels=c("A")),
                 plot_grid(pb, pc, pd, labels=c("B", "C", "D"), ncol=3, rel_widths = c(0.24, 0.19, 0.57)),
                 nrow=2, rel_heights=c(0.45, 0.55))
  
  ggsave("figures/figure2.pdf", plot=p, width=10.5, height=7.5, device=cairo_pdf)
}


# Fig S2 panel
generateFigureS2 <- function() {
  pa <- plotTdpDetails(csv="rosmap.csv", region="DLPFC") + ylim(c(-2.9, 4.8))
  pbfull <- plotMultiVarForest(csv="rosmap.csv", Neu=FALSE, standardize=TRUE)
  pTau <- paste("p == ", formatC(pbfull$data["tangles_sqrt", "Pr(>|t|)"], digits=2), sep="")
  pb <- pbfull$p + annotate("text", x=-0.28, y=9.4, label=pTau, parse=TRUE, col="royalblue")
  
  p <- plot_grid(pa, pb, labels=c("A", "B"), ncol=2, rel_widths = c(0.33, 0.66))
  
  ggsave("figures/figureS2.pdf", plot=p, width=10, height=5, device=cairo_pdf)
}


generateFigure2()
generateFigureS2()
tableS3 <- generateUniPathoTable()
write.csv(tableS3, file="tables/tableS3.csv", row.names=FALSE)
