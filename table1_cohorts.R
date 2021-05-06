library(censReg)

# Table 1 - Characteristics of the ROSMAP cohort
generateRosmapTable <- function (csv="rosmap.csv") {
  
  rosmap <- read.csv(csv, colClasses=c(ProjID="character", apoe_genotype="character"),
                     stringsAsFactors=FALSE)
  regions <- split(rosmap, rosmap$BrainRegion)
  regions <- regions[c("DLPFC", "PCC", "CB")]
  
  rmTable <- as.data.frame(matrix("", nrow=14, ncol=4), stringsAsFactors=FALSE)
  colnames(rmTable) <- c("Variable", names(regions))
  
  # Sex 
  rmTable[1, "Variable"] <- "Sex (male)"
  rmTable[1, "DLPFC"] <- paste(sum(regions[["DLPFC"]]$msex), " (", round(sum(regions[["DLPFC"]]$msex)/nrow(regions[["DLPFC"]]), digits=3) * 100, "%)", sep="")
  rmTable[1, "PCC"] <- paste(sum(regions[["PCC"]]$msex), " (", round(sum(regions[["PCC"]]$msex)/nrow(regions[["PCC"]]), digits=3) * 100, "%)", sep="")
  rmTable[1, "CB"] <- paste(sum(regions[["CB"]]$msex), " (", round(sum(regions[["CB"]]$msex)/nrow(regions[["CB"]]), digits=3) * 100, "%)", sep="")
  
  # Age
  rmTable[2, "Variable"] <- "Age (years)"
  rmTable[2, "DLPFC"] <- paste(round(mean(regions[["DLPFC"]]$age_death), digits=1), " (", round(sd(regions[["DLPFC"]]$age_death), digits=1), ")", sep="")
  rmTable[2, "PCC"] <- paste(round(mean(regions[["PCC"]]$age_death), digits=1), " (", round(sd(regions[["PCC"]]$age_death), digits=1), ")", sep="")
  rmTable[2, "CB"] <- paste(round(mean(regions[["CB"]]$age_death), digits=1), " (", round(sd(regions[["CB"]]$age_death), digits=1), ")", sep="")
  
  # Pathologic AD diagnosis
  rmTable[3, "Variable"] <- "Pathologic AD"
  rmTable[3, "DLPFC"] <- paste(sum(regions[["DLPFC"]]$pathoAD), " (", round(sum(regions[["DLPFC"]]$pathoAD)/nrow(regions[["DLPFC"]]), digits=3) * 100, "%)", sep="")
  rmTable[3, "PCC"] <- paste(sum(regions[["PCC"]]$pathoAD), " (", round(sum(regions[["PCC"]]$pathoAD)/nrow(regions[["PCC"]]), digits=3) * 100, "%)", sep="")
  rmTable[3, "CB"] <- paste(sum(regions[["CB"]]$pathoAD), " (", round(sum(regions[["CB"]]$pathoAD)/nrow(regions[["CB"]]), digits=3) * 100, "%)", sep="")
  
  # Amyloid
  rmTable[4, "Variable"] <- "Amyloid (% area affected)"
  rmTable[4, "DLPFC"] <- paste(round(mean(regions[["DLPFC"]]$amyloid_sqrt^2, na.rm=TRUE), digits=1), " (", round(sd(regions[["DLPFC"]]$amyloid_sqrt^2, na.rm=TRUE), digits=1), ")", sep="")
  rmTable[4, "PCC"] <- paste(round(mean(regions[["PCC"]]$amyloid_sqrt^2, na.rm=TRUE), digits=1), " (", round(sd(regions[["PCC"]]$amyloid_sqrt^2, na.rm=TRUE), digits=1), ")", sep="")
  rmTable[4, "CB"] <- paste(round(mean(regions[["CB"]]$amyloid_sqrt^2, na.rm=TRUE), digits=1), " (", round(sd(regions[["CB"]]$amyloid_sqrt^2, na.rm=TRUE), digits=1), ")", sep="")
  
  # Tau
  rmTable[5, "Variable"] <- "Tau (% area affected)"
  rmTable[5, "DLPFC"] <- paste(round(mean(regions[["DLPFC"]]$tangles_sqrt^2, na.rm=TRUE), digits=1), " (", round(sd(regions[["DLPFC"]]$tangles_sqrt^2, na.rm=TRUE), digits=1), ")", sep="")
  rmTable[5, "PCC"] <- paste(round(mean(regions[["PCC"]]$tangles_sqrt^2, na.rm=TRUE), digits=1), " (", round(sd(regions[["PCC"]]$tangles_sqrt^2, na.rm=TRUE), digits=1), ")", sep="")
  rmTable[5, "CB"] <- paste(round(mean(regions[["CB"]]$tangles_sqrt^2, na.rm=TRUE), digits=1), " (", round(sd(regions[["CB"]]$tangles_sqrt^2, na.rm=TRUE), digits=1), ")", sep="")
  
  # TDP-43
  rmTable[6, "Variable"] <- "TDP-43 (none, amygdala, limbic, neocortical)"
  rmTable[6, "DLPFC"] <- paste(paste(table(regions[["DLPFC"]]$tdp_stage4), collapse=", "), " (",
                         paste(paste(round(table(regions[["DLPFC"]]$tdp_stage4)/sum(!is.na(regions[["DLPFC"]]$tdp_stage4)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[6, "PCC"] <- paste(paste(table(regions[["PCC"]]$tdp_stage4), collapse=", "), " (",
                             paste(paste(round(table(regions[["PCC"]]$tdp_stage4)/sum(!is.na(regions[["PCC"]]$tdp_stage4)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[6, "CB"] <- paste(paste(table(regions[["CB"]]$tdp_stage4), collapse=", "), " (",
                               paste(paste(round(table(regions[["CB"]]$tdp_stage4)/sum(!is.na(regions[["CB"]]$tdp_stage4)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  
  # Lewy bodies
  rmTable[7, "Variable"] <- "Lewy bodies (none, nigral, limbic, neocortical)"
  rmTable[7, "DLPFC"] <- paste(paste(table(regions[["DLPFC"]]$dlbdx), collapse=", "), " (",
                               paste(paste(round(table(regions[["DLPFC"]]$dlbdx)/sum(!is.na(regions[["DLPFC"]]$dlbdx)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[7, "PCC"] <- paste(paste(table(regions[["PCC"]]$dlbdx), collapse=", "), " (",
                             paste(paste(round(table(regions[["PCC"]]$dlbdx)/sum(!is.na(regions[["PCC"]]$dlbdx)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[7, "CB"] <- paste(paste(table(regions[["CB"]]$dlbdx), collapse=", "), " (",
                             paste(paste(round(table(regions[["CB"]]$dlbdx)/sum(!is.na(regions[["CB"]]$dlbdx)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  
  # Cerebral amyloid angiopathy
  rmTable[8, "Variable"] <- "Cerebral amyloid angiopathy (none, mild, moderate, severe)"
  rmTable[8, "DLPFC"] <- paste(paste(table(regions[["DLPFC"]]$caa_4gp), collapse=", "), " (",
                               paste(paste(round(table(regions[["DLPFC"]]$caa_4gp)/sum(!is.na(regions[["DLPFC"]]$caa_4gp)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[8, "PCC"] <- paste(paste(table(regions[["PCC"]]$caa_4gp), collapse=", "), " (",
                             paste(paste(round(table(regions[["PCC"]]$caa_4gp)/sum(!is.na(regions[["PCC"]]$caa_4gp)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[8, "CB"] <- paste(paste(table(regions[["CB"]]$caa_4gp), collapse=", "), " (",
                             paste(paste(round(table(regions[["CB"]]$caa_4gp)/sum(!is.na(regions[["CB"]]$caa_4gp)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  
  # Cerebral atherosclerosis
  rmTable[9, "Variable"] <- "Cerebral atherosclerosis (none, mild, moderate, severe)"
  rmTable[9, "DLPFC"] <- paste(paste(table(regions[["DLPFC"]]$cvda_4gp2), collapse=", "), " (",
                               paste(paste(round(table(regions[["DLPFC"]]$cvda_4gp2)/sum(!is.na(regions[["DLPFC"]]$cvda_4gp2)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[9, "PCC"] <- paste(paste(table(regions[["PCC"]]$cvda_4gp2), collapse=", "), " (",
                             paste(paste(round(table(regions[["PCC"]]$cvda_4gp2)/sum(!is.na(regions[["PCC"]]$cvda_4gp2)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[9, "CB"] <- paste(paste(table(regions[["CB"]]$cvda_4gp2), collapse=", "), " (",
                             paste(paste(round(table(regions[["CB"]]$cvda_4gp2)/sum(!is.na(regions[["CB"]]$cvda_4gp2)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  
  # Arteriolosclerosis
  rmTable[10, "Variable"] <- "Arteriolosclerosis (none, mild, moderate, severe)"
  rmTable[10, "DLPFC"] <- paste(paste(table(regions[["DLPFC"]]$arteriol_scler), collapse=", "), " (",
                                paste(paste(round(table(regions[["DLPFC"]]$arteriol_scler)/sum(!is.na(regions[["DLPFC"]]$arteriol_scler)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[10, "PCC"] <- paste(paste(table(regions[["PCC"]]$arteriol_scler), collapse=", "), " (",
                              paste(paste(round(table(regions[["PCC"]]$arteriol_scler)/sum(!is.na(regions[["PCC"]]$arteriol_scler)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  rmTable[10, "CB"] <- paste(paste(table(regions[["CB"]]$arteriol_scler), collapse=", "), " (",
                              paste(paste(round(table(regions[["CB"]]$arteriol_scler)/sum(!is.na(regions[["CB"]]$arteriol_scler)), digits=3)*100, "%", sep="") , collapse=", "), ")", sep="")
  
  # Gross chronic infarcts
  rmTable[11, "Variable"] <- "Gross chronic infarcts (one or more)"
  rmTable[11, "DLPFC"] <- paste(sum(regions[["DLPFC"]]$ci_num2_gct, na.rm=TRUE), " (", round(sum(regions[["DLPFC"]]$ci_num2_gct, na.rm=TRUE)/sum(!is.na(regions[["DLPFC"]]$ci_num2_gct)), digits=3) * 100, "%)", sep="")
  rmTable[11, "PCC"] <- paste(sum(regions[["PCC"]]$ci_num2_gct, na.rm=TRUE), " (", round(sum(regions[["PCC"]]$ci_num2_gct, na.rm=TRUE)/sum(!is.na(regions[["PCC"]]$ci_num2_gct)), digits=3) * 100, "%)", sep="")
  rmTable[11, "CB"] <- paste(sum(regions[["CB"]]$ci_num2_gct, na.rm=TRUE), " (", round(sum(regions[["CB"]]$ci_num2_gct, na.rm=TRUE)/sum(!is.na(regions[["CB"]]$ci_num2_gct)), digits=3) * 100, "%)", sep="")
  
  # Chronic microinfarcts
  rmTable[12, "Variable"] <- "Chronic microinfarcts (one or more)"
  rmTable[12, "DLPFC"] <- paste(sum(regions[["DLPFC"]]$ci_num2_mct, na.rm=TRUE), " (", round(sum(regions[["DLPFC"]]$ci_num2_mct, na.rm=TRUE)/sum(!is.na(regions[["DLPFC"]]$ci_num2_mct)), digits=3) * 100, "%)", sep="")
  rmTable[12, "PCC"] <- paste(sum(regions[["PCC"]]$ci_num2_mct, na.rm=TRUE), " (", round(sum(regions[["PCC"]]$ci_num2_mct, na.rm=TRUE)/sum(!is.na(regions[["PCC"]]$ci_num2_mct)), digits=3) * 100, "%)", sep="")
  rmTable[12, "CB"] <- paste(sum(regions[["CB"]]$ci_num2_mct, na.rm=TRUE), " (", round(sum(regions[["CB"]]$ci_num2_mct, na.rm=TRUE)/sum(!is.na(regions[["CB"]]$ci_num2_mct)), digits=3) * 100, "%)", sep="")
  
  # Hippocampal sclerosis
  rmTable[13, "Variable"] <- "Hippocampal sclerosis (present)"
  rmTable[13, "DLPFC"] <- paste(sum(regions[["DLPFC"]]$hspath_typ, na.rm=TRUE), " (", round(sum(regions[["DLPFC"]]$hspath_typ, na.rm=TRUE)/sum(!is.na(regions[["DLPFC"]]$hspath_typ)), digits=3) * 100, "%)", sep="")
  rmTable[13, "PCC"] <- paste(sum(regions[["PCC"]]$hspath_typ, na.rm=TRUE), " (", round(sum(regions[["PCC"]]$hspath_typ, na.rm=TRUE)/sum(!is.na(regions[["PCC"]]$hspath_typ)), digits=3) * 100, "%)", sep="")
  rmTable[13, "CB"] <- paste(sum(regions[["CB"]]$hspath_typ, na.rm=TRUE), " (", round(sum(regions[["CB"]]$hspath_typ, na.rm=TRUE)/sum(!is.na(regions[["CB"]]$hspath_typ)), digits=3) * 100, "%)", sep="")
  
  # Post mortem interval
  rmTable[14, "Variable"] <- "Post mortem interval (hours)"
  rmTable[14, "DLPFC"] <- paste(round(mean(regions[["DLPFC"]]$pmi, na.rm=TRUE), digits=1), " (", round(sd(regions[["DLPFC"]]$pmi, na.rm=TRUE), digits=1), ")", sep="")
  rmTable[14, "PCC"] <- paste(round(mean(regions[["PCC"]]$pmi, na.rm=TRUE), digits=1), " (", round(sd(regions[["PCC"]]$pmi, na.rm=TRUE), digits=1), ")", sep="")
  rmTable[14, "CB"] <- paste(round(mean(regions[["CB"]]$pmi, na.rm=TRUE), digits=1), " (", round(sd(regions[["CB"]]$pmi, na.rm=TRUE), digits=1), ")", sep="")
  
  colnames(rmTable) <- c("Variable", paste(names(regions), " (n=", sapply(regions, nrow), ")", sep=""))
  
  return(rmTable)
}


# Table S1 - Characteristics of the Mayo cohort
generateMayoTable <- function (csv="mayo.csv") {
  
  mayo <- read.csv(csv, colClasses=c(SampleID="character", apoe_genotype="character"),
                   stringsAsFactors=FALSE)
  
  mayoTable <- as.data.frame(matrix("", nrow=3, ncol=5), stringsAsFactors=FALSE)
  colnames(mayoTable) <- c("Variable", "Control", "AD", "PSP", "Path. aging")
  
  # Sex 
  mayoTable[1, "Variable"] <- "Sex (male)"
  mayoTable[1, "Control"] <- paste(sum(mayo$Sex[mayo$Diagnosis == "Control"] == "M"), " (", round(sum(mayo$Sex[mayo$Diagnosis == "Control"] == "M")/sum(mayo$Diagnosis == "Control"), digits=3) * 100, "%)", sep="")
  mayoTable[1, "AD"] <- paste(sum(mayo$Sex[mayo$Diagnosis == "AD"] == "M"), " (", round(sum(mayo$Sex[mayo$Diagnosis == "AD"] == "M")/sum(mayo$Diagnosis == "AD"), digits=3) * 100, "%)", sep="")
  mayoTable[1, "PSP"] <- paste(sum(mayo$Sex[mayo$Diagnosis == "PSP"] == "M"), " (", round(sum(mayo$Sex[mayo$Diagnosis == "PSP"] == "M")/sum(mayo$Diagnosis == "PSP"), digits=3) * 100, "%)", sep="")
  mayoTable[1, "Path. aging"] <- paste(sum(mayo$Sex[mayo$Diagnosis == "PathologicAging"] == "M"), " (", round(sum(mayo$Sex[mayo$Diagnosis == "PathologicAging"] == "M")/sum(mayo$Diagnosis == "PathologicAging"), digits=3) * 100, "%)", sep="")
  
  # Age - use simple Tobit regression to account for right-censored age at death
  mayoTable[2, "Variable"] <- "Age (years)"
  age <- mayo$AgeAtDeath[mayo$Diagnosis == "Control"]
  age[age == "90+"] <- 90
  age <- as.numeric(age)
  coefs <- coef(censReg(age ~ 1, left=-Inf, right=90))
  mayoTable[2, "Control"] <- paste(round(coefs["(Intercept)"], digits=1), " (", round(exp(coefs["logSigma"]), digits=1), ")", sep="")
  
  age <- mayo$AgeAtDeath[mayo$Diagnosis == "AD"]
  age[age == "90+"] <- 90
  age <- as.numeric(age)
  coefs <- coef(censReg(age ~ 1, left=-Inf, right=90))
  mayoTable[2, "AD"] <- paste(round(coefs["(Intercept)"], digits=1), " (", round(exp(coefs["logSigma"]), digits=1), ")", sep="")
  
  age <- mayo$AgeAtDeath[mayo$Diagnosis == "PSP"]
  stopifnot(all(age != "90+"))
  age <- as.numeric(age)
  mayoTable[2, "PSP"] <- paste(round(mean(age), digits=1), " (", round(sd(age), digits=1), ")", sep="")
  
  age <- mayo$AgeAtDeath[mayo$Diagnosis == "PathologicAging"]
  age[age == "90+"] <- 90
  age <- as.numeric(age)
  coefs <- coef(censReg(age ~ 1, left=-Inf, right=90))
  mayoTable[2, "Path. aging"] <- paste(round(coefs["(Intercept)"], digits=1), " (", round(exp(coefs["logSigma"]), digits=1), ")", sep="")
  
  # Post mortem interval
  mayoTable[3, "Variable"] <- "Post mortem interval (hours)"
  mayoTable[3, "Control"] <- paste(round(mean(mayo$PMI[mayo$Diagnosis == "Control"], na.rm=TRUE), digits=1), " (", round(sd(mayo$PMI[mayo$Diagnosis == "Control"], na.rm=TRUE), digits=1), ")", sep="")
  mayoTable[3, "AD"] <- paste(round(mean(mayo$PMI[mayo$Diagnosis == "AD"], na.rm=TRUE), digits=1), " (", round(sd(mayo$PMI[mayo$Diagnosis == "AD"], na.rm=TRUE), digits=1), ")", sep="")
  mayoTable[3, "PSP"] <- paste(round(mean(mayo$PMI[mayo$Diagnosis == "PSP"], na.rm=TRUE), digits=1), " (", round(sd(mayo$PMI[mayo$Diagnosis == "PSP"], na.rm=TRUE), digits=1), ")", sep="")
  stopifnot(all(is.na(mayo$PMI[mayo$Diagnosis == "PathologicAging"])))
  mayoTable[3, "Path. aging"] <- "-"
  
  colnames(mayoTable) <- c("Variable", paste(colnames(mayoTable)[2:5], " (n=", table(mayo$Diagnosis)[c("Control", "AD", "PSP", "PathologicAging")], ")", sep=""))
  
  return(mayoTable)
}


# Table S2 - Characteristics of the MSBB cohort
generateMsbbTable <- function (csv="msbb.csv") {
  
  msbb <- read.csv(csv, colClasses=c(SampleID="character", apoe_genotype="character"),
                   stringsAsFactors=FALSE)
  
  msbbTable <- as.data.frame(matrix("", nrow=3, ncol=3), stringsAsFactors=FALSE)
  colnames(msbbTable) <- c("Variable", "Control", "AD")
  
  # Sex 
  msbbTable[1, "Variable"] <- "Sex (male)"
  msbbTable[1, "Control"] <- paste(sum(msbb$Sex[msbb$pathoAD == 0] == "M"), " (", round(sum(msbb$Sex[msbb$pathoAD == 0] == "M")/sum(msbb$pathoAD == 0), digits=3) * 100, "%)", sep="")
  msbbTable[1, "AD"] <- paste(sum(msbb$Sex[msbb$pathoAD == 1] == "M"), " (", round(sum(msbb$Sex[msbb$pathoAD == 1] == "M")/sum(msbb$pathoAD == 1), digits=3) * 100, "%)", sep="")
  
  # Age - use simple Tobit regression to account for right-censored age at death
  msbbTable[2, "Variable"] <- "Age (years)"
  age <- msbb$AgeAtDeath[msbb$pathoAD == 0]
  age[age == "90+"] <- 90
  age <- as.numeric(age)
  coefs <- coef(censReg(age ~ 1, left=-Inf, right=90))
  msbbTable[2, "Control"] <- paste(round(coefs["(Intercept)"], digits=1), " (", round(exp(coefs["logSigma"]), digits=1), ")", sep="")
  age <- msbb$AgeAtDeath[msbb$pathoAD == 1]
  age[age == "90+"] <- 90
  age <- as.numeric(age)
  coefs <- coef(censReg(age ~ 1, left=-Inf, right=90))
  msbbTable[2, "AD"] <- paste(round(coefs["(Intercept)"], digits=1), " (", round(exp(coefs["logSigma"]), digits=1), ")", sep="")
  
  # Post mortem interval (given in minutes for MSBB)
  msbbTable[3, "Variable"] <- "Post mortem interval (hours)"
  msbbTable[3, "Control"] <- paste(round(mean(msbb$PMI[msbb$pathoAD == 0]/60, na.rm=TRUE), digits=1), " (", round(sd(msbb$PMI[msbb$pathoAD == 0]/60, na.rm=TRUE), digits=1), ")", sep="")
  msbbTable[3, "AD"] <- paste(round(mean(msbb$PMI[msbb$pathoAD == 1]/60, na.rm=TRUE), digits=1), " (", round(sd(msbb$PMI[msbb$pathoAD == 1]/60, na.rm=TRUE), digits=1), ")", sep="")
  
  colnames(msbbTable) <- c("Variable", paste(colnames(msbbTable)[2:3], " (n=", table(msbb$pathoAD)[c("0", "1")], ")", sep=""))
  
  return(msbbTable)
}


table1 <- generateRosmapTable()
write.csv(table1, file="tables/table1.csv", row.names=FALSE)
tableS1 <- generateMayoTable()
write.csv(tableS1, file="tables/tableS1.csv", row.names=FALSE)
tableS2 <- generateMsbbTable()
write.csv(tableS2, file="tables/tableS2.csv", row.names=FALSE)
