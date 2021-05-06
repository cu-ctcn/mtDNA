# Write a simple text file with all vcf files to be merged in the next step

load("/mnt/mfs/ctcn/datasets/rosmap/wgs/ampad/experimentalData/sampleSheet_1200.Rdata")
vcfFiles <- file.path("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf",
                      sampleSheet$SampleID, paste(sampleSheet$SampleID, "_final.vcf.gz", sep=""))
stopifnot(all(file.exists(vcfFiles)))
sampleSheet$vcfFiles=vcfFiles

writeLines(sampleSheet$vcfFiles, con=paste("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/vcfFileLists/vcfFiles.txt"))

sampleSheet$Tissue <- as.character(sampleSheet$Tissue)
sampleSheet$Tissue[is.na(sampleSheet$Tissue)] <- "NA"

samples <- split.data.frame(sampleSheet, sampleSheet$Tissue)
for (tissue in names(samples)) {
  writeLines(samples[[tissue]]$vcfFiles, con=paste("/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/vcfFileLists/vcfFiles_", tissue, ".txt", sep=""))
}
