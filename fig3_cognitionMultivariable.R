library(ggplot2)
library(cowplot)


# Figure 3A-B - forest plot showing coefficients from multivariable model with cognition as outcome
# for 3A set Neu=FALSE, for 3B set Neu=TRUE
plotForestCognition <- function (csv="rosmap.csv", Neu=FALSE, standardize=TRUE) {
  
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
  
  
  ind <- !is.na(rosmap$mtDNAcn) & !is.na(rosmap$age_death) & !is.na(rosmap$msex) & !is.na(rosmap$educ) & !is.na(rosmap$amyloid_sqrt) & !is.na(rosmap$tangles_sqrt) &
    !is.na(rosmap$tdp_simple) & !is.na(rosmap$dlb_simple) & !is.na(rosmap$ci_num2_gct) & !is.na(rosmap$ci_num2_mct) & !is.na(rosmap$cvda_simple) &
    !is.na(rosmap$arteriol_scler_simple) & !is.na(rosmap$caa_simple) & !is.na(rosmap$hspath_typ)
  
  if (Neu) {
    ind <- ind &  !is.na(rosmap$Neu)
    rosmap <- rosmap[ind, ]
    model <- cogn_global_lv ~ age_death + msex + educ + amyloid_sqrt + tangles_sqrt + tdp_simple + dlb_simple + ci_num2_gct +
      ci_num2_mct + cvda_simple + arteriol_scler_simple + caa_simple + hspath_typ + mtDNAcn + Neu
    modelRed <- cogn_global_lv ~ age_death + msex + educ + amyloid_sqrt + tangles_sqrt + tdp_simple + dlb_simple + ci_num2_gct +
      ci_num2_mct + cvda_simple + arteriol_scler_simple + caa_simple + hspath_typ + Neu
  } else {
    rosmap <- rosmap[ind, ]
    model <- cogn_global_lv ~ age_death + msex + educ + amyloid_sqrt + tangles_sqrt + tdp_simple + dlb_simple + ci_num2_gct +
      ci_num2_mct + cvda_simple + arteriol_scler_simple + caa_simple + hspath_typ + mtDNAcn
    modelRed <- cogn_global_lv ~ age_death + msex + educ + amyloid_sqrt + tangles_sqrt + tdp_simple + dlb_simple + ci_num2_gct +
      ci_num2_mct + cvda_simple + arteriol_scler_simple + caa_simple + hspath_typ
  }
  
  # Note: There shoudn't be anymore NAs in any variable
  n <- c(age_death=sum(!is.na(rosmap$age_death)),
         msex1=sum(rosmap$msex == "1"),
         educ=sum(!is.na(rosmap$educ)),
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
         mtDNAcn=sum(!is.na(rosmap$mtDNAcn)),
         Neu=sum(!is.na(rosmap$Neu)))
  
  labels <- c(age_death="Age",
              msex1=paste("Sex\n(male, n=", n["msex1"], ")", sep=""),
              educ="Education",
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
  fitR <- lm(modelRed, rosmap)
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
    xlab("Effect size (standard deviations of cognition)") +
    ylab("") +
    theme_light()
  
  return(list(p=p, data=data, fit=fit, fitR=fitR))
}


# Figure 3 panel
generateFigure3 <- function() {
  pa <- plotForestCognition(csv="rosmap.csv", standardize=TRUE, Neu=FALSE)
  pb <- plotForestCognition(csv="rosmap.csv", standardize=TRUE, Neu=TRUE)
  
  # add p-values to figure A
  pEdu <- paste("p == ", formatC(pa$data["educ", "Pr(>|t|)"], digits=2), sep="")
  pTau <- paste("p == ", gsub("e-", "%*%10^-", formatC(pa$data["tangles_sqrt", "Pr(>|t|)"], digits=2)), sep="")
  pLewy <- paste("p == ", gsub("e-0", "%*%10^-", formatC(pa$data["dlb_simpleTRUE", "Pr(>|t|)"], digits=2)), sep="")
  pGI <- paste("p == ", gsub("e-0", "%*%10^-", formatC(pa$data["ci_num2_gct1", "Pr(>|t|)"], digits=2)), sep="")
  pHs <- paste("p == ", gsub("e-0", "%*%10^-", formatC(pa$data["hspath_typ1", "Pr(>|t|)"], digits=2)), sep="")
  pmtDNAcn <- paste("p == ", formatC(pa$data["mtDNAcn", "Pr(>|t|)"], digits=2), sep="")
  plotA <- pa$p + 
    annotate("text", x=0.135, y=12.4, label=pEdu, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.69, y=10.4, label=pTau, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.59, y=9.4, label=pLewy, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.51, y=7.4, label=pGI, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.74, y=2.4, label=pHs, parse=TRUE, col="royalblue") +
    annotate("text", x=0.10, y=1.4, label=pmtDNAcn, parse=TRUE, col="royalblue")
  
  # add p-values to figure B
  pEdu <- paste("p == ", formatC(pb$data["educ", "Pr(>|t|)"], digits=2), sep="")
  pTau <- paste("p == ", gsub("e-", "%*%10^-", formatC(pb$data["tangles_sqrt", "Pr(>|t|)"], digits=2)), sep="")
  pLewy <- paste("p == ", gsub("e-0", "%*%10^-", formatC(pb$data["dlb_simpleTRUE", "Pr(>|t|)"], digits=2)), sep="")
  pGI <- paste("p == ", gsub("e-0", "%*%10^-", formatC(pb$data["ci_num2_gct1", "Pr(>|t|)"], digits=2)), sep="")
  pHs <- paste("p == ", gsub("e-0", "%*%10^-", formatC(pb$data["hspath_typ1", "Pr(>|t|)"], digits=2)), sep="")
  pmtDNAcn <- paste("p == ", formatC(pb$data["mtDNAcn", "Pr(>|t|)"], digits=2), sep="")
  plotB <- pb$p + 
    annotate("text", x=0.135, y=13.4, label=pEdu, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.65, y=11.4, label=pTau, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.58, y=10.4, label=pLewy, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.51, y=8.4, label=pGI, parse=TRUE, col="royalblue") +
    annotate("text", x=-0.65, y=3.4, label=pHs, parse=TRUE, col="royalblue") +
    annotate("text", x=0.12, y=2.4, label=pmtDNAcn, parse=TRUE, col="royalblue")
  
  p <- plot_grid(plotA, plotB, nrow=1, labels=c("A", "B"))
  ggsave("figures/figure3.pdf", plot=p, width=11, height=5, device=cairo_pdf)
  
  # summary(pa$fit)$r.squared - summary(pa$fitR)$r.squared # 0.01801835
  # summary(pb$fit)$r.squared - summary(pb$fitR)$r.squared # 0.008106487
  
  # length(pa$fit$residuals) # n=393
  # length(pb$fit$residuals) # n=287
  
}

generateFigure3()
