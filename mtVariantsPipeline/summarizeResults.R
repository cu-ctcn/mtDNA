library(VariantAnnotation)

variants <- read.table("/mnt/mfs/ctcn/datasets/rosmap/wgs/ampad/experimentalData/sampleSheet_1200.csv",
                       sep=",", header=TRUE, colClasses=c(ProjID="character"), stringsAsFactors=FALSE)


# Add haplogroups and contamination
variants$Contamination=NA
variants$MajorHg=NA
variants$MajorLevel=NA
variants$MinorHg=NA
variants$MinorLevel=NA
for (i in 1:nrow(variants)) {
  conta <- read.csv(file=file.path("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam",
                                   variants$SampleID[i], "haplochecker",
                                   paste(variants$SampleID[i], "_MT.contamination.txt", sep="")),
                    header=TRUE, stringsAsFactors=FALSE, sep="\t")
  variants$Contamination[i] <- conta$Contamination
  variants$MajorHg[i] <- conta$MajorHG
  variants$MajorLevel[i] <- conta$MajorLevel
  variants$MinorHg[i] <- conta$MinorHG
  variants$MinorLevel[i] <- conta$MinorLevel
}
haploGrp <- substr(variants$MajorHg, 1, 2)
indSupGrp <- grepl("^[A-Za-z]+$", haploGrp) & nchar(haploGrp) == 2 # Keep HV and JT as separate groups
haploGrp[!indSupGrp] <- substr(haploGrp[!indSupGrp], 1, 1)
variants$HaploGrp <- haploGrp


# Add number of variants
vcf <- readVcf("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz")
heteroplasmyTh <- 0.9
relAlFreq <- apply(geno(vcf)$AF, 2, function(x) {
  freq <- rep(0, length(x))
  ind <- sapply(lapply(x, is.na), all)
  freq[!ind] <- sapply(x[!ind], function (a) {
    sum(a[!is.na(a)])
  })
  return(freq)
})

relAlFreq <- relAlFreq[, variants$SampleID]
variants$NumHomoplasmy <- apply(relAlFreq >= heteroplasmyTh, 2, sum)
variants$NumHeteroplasmy <- apply(relAlFreq < heteroplasmyTh & relAlFreq != 0, 2, sum)


save(variants, file="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/values/variants.Rdata")
write.csv(variants, file="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/values/variants.csv",
          row.names=FALSE, quote=FALSE)
