# Run after collectMetrics.sh has been run for all samples.
# Collects results from haplochecker and gets the median
# autosomal coverage. The generated csv file is used by
# callVariants.sh to filter variants.

samples <- read.csv("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/bam/bamfiles.csv",
                    header=FALSE, stringsAsFactors=FALSE)[, 1]
projDir <- "/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2"

mtMetrics <- data.frame(Sample=samples,
                        Contamination=NA,
                        MajorHg=NA, MajorLevel=NA,
                        MinorHg=NA, MinorLevel=NA,
                        AutosomalCov=NA,
                        MtCov=NA,
                        stringsAsFactors=FALSE)

# Get contamination from haplochecker
for (i in 1:nrow(mtMetrics)) {
  conta <- read.csv(file=file.path(projDir, "bam", mtMetrics$Sample[i], "haplochecker",
                                   paste(mtMetrics$Sample[i], "_MT.contamination.txt", sep="")),
                    header=TRUE, stringsAsFactors=FALSE, sep="\t")
  mtMetrics$Contamination[i] <- conta$Contamination
  mtMetrics$MajorHg[i] <- conta$MajorHG
  mtMetrics$MajorLevel[i] <- conta$MajorLevel
  mtMetrics$MinorHg[i] <- conta$MinorHG
  mtMetrics$MinorLevel[i] <- conta$MinorLevel
}

# Add coverages (previously calculated for MT copy number using masked genomes)
load("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtCopyNumber/coverage/Rdata/copyNumbers.Rdata")
ind <- match(mtMetrics$Sample, copyNumbers$SampleID)
mtMetrics$AutosomalCov <- round(copyNumbers$autosomalCov[ind])
mtMetrics$MtCov <- round((copyNumbers$MTcn[ind] / 2) * copyNumbers$autosomalCov[ind])

write.csv(mtMetrics, file="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/mtMetrics.csv",
          row.names=FALSE, quote=FALSE)
