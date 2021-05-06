library(VariantAnnotation)
library(SummarizedExperiment)
library(meta)


# Reads in the genotypes for the mtDNAcn risk loci (primary analysis) from Longchamps et al.
# We don't read in rs28665408 (18:67516846_C/A) since we observed multiple indels at this position.
# We also don't read in rs3218221 (19:10676941_G/A) because it very rare (3 in DLPFC, 0 in CB, 0 in FP, 1 in TCX)
# Note that "ALLELE1" is the reference allele and used as effect allele in the original file. We change the effect
# allele to the alternative allele "ALLELE0".
# The original Table S6 from Longchamps et al. containing the risk loci can be downloaded here:
# https://www.biorxiv.org/content/biorxiv/early/2021/01/28/2021.01.25.428086/DC2/embed/media-2.xlsx?download=true
# The vcf files with the genotypes of the samples used by our study can be downloaded from the AD knowledge portal:
# ROSMAP: https://www.synapse.org/#!Synapse:syn11724057
# Mayo: https://www.synapse.org/#!Synapse:syn11724002
# MSBB: https://www.synapse.org/#!Synapse:syn11723899
getGwasGenotypes <- function(samples, study="rosmap") {
   
  gwas <- read.csv("data/LongchampsTableS6.csv", stringsAsFactors=FALSE, header=TRUE, skip=1)
  gwas <- gwas[gwas$GWAS == "primary", ]
  gwas <- gwas[gwas$SNP != "rs28665408", ] # exclude SNP with complex indels
  gwas <- gwas[gwas$SNP != "rs3218221", ]  # exclude rare SNP not present in most of our datasets
  gwas$BP <- as.numeric(gsub(",", "", gwas$BP))
  gwas$EFFECT.SIZE.ESTIMATE <- gwas$EFFECT.SIZE.ESTIMATE * -1
  gwas <- GRanges(IRanges(gwas$BP, width=1), seqnames=gwas$CHR, gwas=gwas)
  names(gwas) <- paste(as.character(seqnames(gwas)), ":", start(gwas), "_",
                       mcols(gwas)$gwas.ALLELE1, "/", mcols(gwas)$gwas.ALLELE0, sep="")

  if (study == "rosmap") {
    vcfDir <- "/mnt/mfs/ctcn/datasets/rosmap/wgs/ampad/variants/snvCombined/"
    filePrefix <- "DEJ_11898_B01_GRM_WGS_2017-05-15_"
  } else if (study == "mayo") {
    vcfDir <- "/mnt/mfs/ctcn/datasets/mayoRnaSeqStudy/wgs/ampad/variants/vcf/"
    filePrefix <- "SCH_11923_B01_GRM_WGS_2017-04-27_"
  } else if (study == "msbb") {
    vcfDir <- "/mnt/mfs/ctcn/datasets/msbb/wgs/ampad/variants/vcf/"
    filePrefix <- "SCH_11923_B02_GRM_WGS_2017-05-15_"
  } else {
    stop("study must be either rosmap, mayo, or msbb")
  }
  genotypes <- numeric()
  
  for (chr in seqlevels(gwas)) {
    tab <- TabixFile(paste(vcfDir, filePrefix, chr, ".recalibrated_variants.vcf.gz",sep=""))
    param <- ScanVcfParam(fixed=c("ALT", "FILTER"), samples=samples, geno="GT", which=gwas[seqnames(gwas) == chr])
    vcf <- readVcf(tab, param=param)
    
    snpNames <- paste(chr, ":", start(gwas[seqnames(gwas) == chr]), "_",
                      mcols(gwas[seqnames(gwas) == chr])$gwas.ALLELE1, "/",
                      mcols(gwas[seqnames(gwas) == chr])$gwas.ALLELE0, sep="")
    stopifnot(all(names(gwas[seqnames(gwas) == chr]) %in% names(vcf)))
    vcf <- vcf[names(gwas[seqnames(gwas) == chr])]
    
    gts <- sapply(as.list(1:nrow(vcf)), function(i) {
      altCode <- which(rowRanges(vcf)$ALT[[i]] == mcols(gwas[seqnames(gwas) == chr])$gwas.ALLELE0[i])
      stopifnot(length(altCode) == 1)
      gt <- rep(NA, ncol(vcf))
      gt[geno(vcf)$GT[i, ] == "0/0"] <- 0
      gt[geno(vcf)$GT[i, ] == paste("0", altCode, sep="/")] <- 1
      gt[geno(vcf)$GT[i, ] == paste(altCode, altCode, sep="/")] <- 2
      return(gt)
    })
    gts <- t(gts)
    colnames(gts) <- colnames(vcf)
    rownames(gts) <- names(gwas[seqnames(gwas) == chr])
    
    genotypes <- rbind(genotypes, gts)
  }
  
  se <- SummarizedExperiment(assays=list(gt=genotypes), rowRanges=gwas)
  return(se)
}


# Run a simple linear regression and return statistics
testSnps <- function(model, gt) {
  gwas <- apply(assay(gt), 1, function (snp) {
    data <- colData(gt)
    data$snp <- snp
    res <- summary(lm(model, data=data))
    
    coefs <- as.vector(t(res$coefficients[match(names(res$aliased), rownames(res$coefficients)), ]))
    raf <- 1 - (sum(data$snp, na.rm=TRUE) / (sum(!is.na(data$snp)) * 2))
    coefs <- c(coefs, raf)
    names(coefs) <- c(paste(rep(c("c", "se", "t", "p"), rep(length(res$aliased))),
                           rep(gsub("\\(|\\)", "", names(res$aliased)), rep(4, length(res$aliased))), sep="."),
                      "RAF")
    
    return(coefs)
  })
  gwas <- t(gwas)
  
  return(gwas[, c("c.snp", "se.snp", "t.snp", "p.snp", "RAF")])
}


# Test for an association between mtDNAcn and loci identified by Longchamps et al. in each brain region. Conduct meta-analysis.
replicateGwas <- function() {
  # rosmap
  rosmap <- read.csv("rosmap.csv", colClasses=c(ProjID="character", apoe_genotype="character"),
                     stringsAsFactors=FALSE)
  rownames(rosmap) <- rosmap$ProjID
   
  # rosmap DLPFC
  dlpfc <- rosmap[rosmap$BrainRegion == "DLPFC", ]
  nrow(dlpfc) # 454
  gt <- getGwasGenotypes(samples=dlpfc$SampleID, study="rosmap")
  stopifnot(all(colnames(gt) == dlpfc$SampleID))
  colnames(gt) <- dlpfc$ProjID
  colData(gt) <- DataFrame(dlpfc)
  model <- mtDNAcn ~ snp + PC1 + PC2 + PC3 + msex + age_death + gpath_sqrt
  gwas <- testSnps(model=model, gt=gt)
  stopifnot(all(rownames(gt) == rownames(gwas)))
  
  gwasTable <- as.data.frame(mcols(gt))
  colnames(gwas) <- paste("dlpfc", colnames(gwas), sep=".")
  gwasTable <- cbind(gwasTable, gwas)
  
 
  # rosmap CB
  cb <- rosmap[rosmap$BrainRegion == "CB", ]
  nrow(cb) # 242
  gt <- getGwasGenotypes(samples=cb$SampleID, study="rosmap")
  stopifnot(all(colnames(gt) == cb$SampleID))
  colnames(gt) <- cb$ProjID
  colData(gt) <- DataFrame(cb)
  model <- mtDNAcn ~ snp + PC1 + PC2 + PC3 + msex + age_death + gpath_sqrt
  gwas <- testSnps(model=model, gt=gt)
  stopifnot(all(rownames(gt) == rownames(gwas)))
  
  colnames(gwas) <- paste("cb", colnames(gwas), sep=".")
  stopifnot(all(rownames(gwas) == rownames(gwasTable)))
  gwasTable <- cbind(gwasTable, gwas)
  
  
  # mayo TCX
  tcx <- read.csv("mayo.csv", colClasses=c(SampleID="character", apoe_genotype="character"),
                  stringsAsFactors=FALSE)
  nrow(tcx) # 262
  gt <- getGwasGenotypes(samples=tcx$SampleID, study="mayo")
  stopifnot(all(colnames(gt) == tcx$SampleID))
  colData(gt) <- DataFrame(tcx)
  model <- mtDNAcn ~ snp + PC1 + PC2 + PC3 + Sex + AgeStrataNumeric + Diagnosis
  gwas <- testSnps(model=model, gt=gt)
  stopifnot(all(rownames(gt) == rownames(gwas)))
  
  colnames(gwas) <- paste("tcx", colnames(gwas), sep=".")
  stopifnot(all(rownames(gwas) == rownames(gwasTable)))
  gwasTable <- cbind(gwasTable, gwas)
  
  
  # msbb FP
  fp <- read.csv("msbb.csv", colClasses=c(SampleID="character", apoe_genotype="character"),
                   stringsAsFactors=FALSE)
  fp <- fp[fp$RaceInferred == "W", ]
  nrow(fp) # 270
  gt <- getGwasGenotypes(samples=fp$SampleID, study="msbb")
  stopifnot(all(colnames(gt) == fp$SampleID))
  colData(gt) <- DataFrame(fp)
  model <- mtDNAcn ~ snp + PC1 + PC2 + PC3 + Sex + AgeStrataNumeric + pathoAD
  gwas <- testSnps(model=model, gt=gt)
  stopifnot(all(rownames(gt) == rownames(gwas)))
  
  colnames(gwas) <- paste("fp", colnames(gwas), sep=".")
  stopifnot(all(rownames(gwas) == rownames(gwasTable)))
  gwasTable <- cbind(gwasTable, gwas)
   
  
  # meta gwas
  gwasTable$meta.r.c <- NA
  gwasTable$meta.r.se <- NA
  gwasTable$meta.r.p <- NA
  for (i in 1:nrow(gwasTable)) {
    ma <- metagen(TE=c(gwasTable$dlpfc.c.snp[i], gwasTable$cb.c.snp[i], gwasTable$tcx.c.snp[i], gwasTable$fp.c.snp[i]),
                  seTE=c(gwasTable$dlpfc.se.snp[i], gwasTable$cb.se.snp[i], gwasTable$tcx.se.snp[i], gwasTable$fp.se.snp[i]),
                  studlab=c("DLPFC", "CB", "TCX", "FP"),
                  comb.fixed=FALSE,
                  comb.random=TRUE,
                  method.tau="REML",
                  sm="MD")
    
    gwasTable$meta.r.c[i] <- ma$TE.random
    gwasTable$meta.r.se[i] <- ma$seTE.random
    gwasTable$meta.r.p[i] <- ma$pval.random
  }
  
  return(gwasTable)
}


# Supplementary Excel File 2 - targeted GWAS results
writeGwasTable <- function() {
  table <- replicateGwas()
  table <- table[c("gwas.SNP", "gwas.CHR", "gwas.BP", "gwas.ALLELE1", "gwas.ALLELE0",
                   "gwas.EFFECT.SIZE.ESTIMATE", "gwas.STANDARD.ERROR", "gwas.P.VALUE", "gwas.A1FREQ",
                   "dlpfc.c.snp", "dlpfc.se.snp", "dlpfc.p.snp", "dlpfc.RAF",
                   "cb.c.snp", "cb.se.snp", "cb.p.snp", "cb.RAF",
                   "tcx.c.snp", "tcx.se.snp", "tcx.p.snp", "tcx.RAF",
                   "fp.c.snp", "fp.se.snp", "fp.p.snp", "fp.RAF",
                   "meta.r.c", "meta.r.se", "meta.r.p")]
  colnames(table) <- c("rsID", "Chr", "Pos", "RefAllele", "AltAllele",
                       "gwas.c", "gwas.se", "gwas.p", "gwas.RefFreq",
                       "dlpfc.c", "dlpfc.se", "dlpfc.p", "dlpfc.RefFreq",
                       "cb.c", "cb.se", "cb.p", "cb.RefFreq",
                       "tcx.c", "tcx.se", "tcx.p", "tcx.RefFreq",
                       "fp.c", "fp.se", "fp.p", "fp.RefFreq",
                       "meta.c", "meta.se", "meta.p")
  refAlleleFreq <- ((table$dlpfc.RefFreq * 454) + (table$cb.RefFreq * 242) +
                   (table$tcx.RefFreq * 262) + (table$fp.RefFreq * 270)) / (454 + 242 + 262 + 270)
  ind <- refAlleleFreq > 0.95 | refAlleleFreq < 0.05
  sum(ind)  # 13 (2 SNPs were already removed when reading the VCF files)
  sum(!ind) # 81
  table <- table[!ind, ]
  table$meta.p.adj <- p.adjust(table$meta.p, method="bonferroni")
  write.csv(table, file="tables/SupExcelFile2.csv", row.names=FALSE)
}

writeGwasTable()