#!/bin/bash

#$-l h_vmem=16G
#$-wd /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2
#$-o /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/logs/combineVariants.o
#$-e /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/logs/combineVariants.e
#$-N combineVariants


source /mnt/mfs/cluster/bin/HgrcPathSetup.sh

bcftools="/mnt/mfs/ctcn/tools/bcftools-1.9/bcftools"

mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged

# All vcf files together does not work - maybe because of too many open file handles?
# Merge vcf files for each tissue (brain region) separately
declare -a tissues=( AnteriorCaudate  OccipitalAssociationCortex Cerebellum PBMC DLPFC PosteriorCingulateCortex FrontalCortex FrontalPole WholeBlood NA)


################################################################
# Filter out all variants that do not have PASS before merging #
################################################################
for t in ${tissues[@]}
do
${bcftools} merge \
  --missing-to-ref \
  --apply-filters PASS \
  --filter-logic + \
  --merge none \
  --output /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_${t}.vcf \
  --output-type v \
  --file-list /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/vcfFileLists/vcfFiles_${t}.txt
bgzip /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_${t}.vcf
${bcftools} index /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_${t}.vcf.gz
done

# Merge all tissue vcfs
${bcftools} merge \
  --missing-to-ref \
  --apply-filters PASS \
  --filter-logic + \
  --merge none \
  --output /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf \
  --output-type v \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_AnteriorCaudate.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_OccipitalAssociationCortex.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_Cerebellum.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_PBMC.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_DLPFC.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_PosteriorCingulateCortex.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_FrontalCortex.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_FrontalPole.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_WholeBlood.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_NA.vcf.gz
bgzip /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf
${bcftools} index /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz
${bcftools} index -t /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants.vcf.gz

##########################################################################################
# No filter - keep all variants but set only PASS if all variant has PASS in all samples #
##########################################################################################
for t in ${tissues[@]}
do
${bcftools} merge \
  --missing-to-ref \
  --filter-logic + \
  --merge none \
  --output /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_${t}_all.vcf \
  --output-type v \
  --file-list /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/vcfFileLists/vcfFiles_${t}.txt
bgzip /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_${t}_all.vcf
${bcftools} index /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_${t}_all.vcf.gz
done

# Merge all tissue vcfs
${bcftools} merge \
  --missing-to-ref \
  --filter-logic + \
  --merge none \
  --output /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_all.vcf \
  --output-type v \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_AnteriorCaudate_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_OccipitalAssociationCortex_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_Cerebellum_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_PBMC_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_DLPFC_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_PosteriorCingulateCortex_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_FrontalCortex_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_FrontalPole_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_WholeBlood_all.vcf.gz \
  /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_NA_all.vcf.gz
bgzip /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_all.vcf
${bcftools} index /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_all.vcf.gz
${bcftools} index -t /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/merged/mtVariants_all.vcf.gz
