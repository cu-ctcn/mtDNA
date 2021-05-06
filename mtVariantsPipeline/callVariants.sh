#!/bin/bash

#$-l h_vmem=16G
#$-t 1-1200
#$-tc 100
#$-wd /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2
#$-N callVariants

source /mnt/mfs/cluster/bin/HgrcPathSetup.sh

offset=0
SGE_TASK_ID=$(( ${SGE_TASK_ID} + ${offset} ))

sampleName=`awk -F "," '{if (NR=='${SGE_TASK_ID}') print $1}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/bam/bamfiles.csv`

gatk="/mnt/mfs/ctcn/tools/gatk-4.1.2.0/gatk"
picard="/mnt/mfs/ctcn/tools/picard_2.20.3/picard.jar"
bcftools="/mnt/mfs/ctcn/tools/bcftools-1.9/bcftools"
tmpDir="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/temp"

mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf
mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/vcf/${sampleName}


#################
# Mutect2 on MT #
#################
## ${gatk} --java-options "-Xms6G" Mutect2 \
${gatk} --java-options "-Xms8G -Djava.io.tmpdir=${tmpDir}" Mutect2 \
  -R /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.fasta \
  -I bam/${sampleName}/${sampleName}_MT.bam \
  -O vcf/${sampleName}/${sampleName}_MT_raw.vcf \
  -L chrM:576-16024 \
  --annotation StrandBiasBySample \
  --mitochondria-mode \
  --max-reads-per-alignment-start 75 \
  --max-mnp-distance 0


#########################
# Mutect2 on shifted MT #
#########################
${gatk} --java-options "-Xms8G -Djava.io.tmpdir=${tmpDir}" Mutect2 \
  -R /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
  -I bam/${sampleName}/${sampleName}_MTshifted8k.bam \
  -O vcf/${sampleName}/${sampleName}_MTshifted8k_raw.vcf \
  -L chrM:8025-9144 \
  --annotation StrandBiasBySample \
  --mitochondria-mode \
  --max-reads-per-alignment-start 75 \
  --max-mnp-distance 0


##########################
# LiftoverAndCombineVcfs #
##########################
java -Xmx10G -jar ${picard} LiftoverVcf \
  I=vcf/${sampleName}/${sampleName}_MTshifted8k_raw.vcf \
  O=vcf/${sampleName}/${sampleName}_MTshiftedBack_raw.vcf \
  R=/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.fasta \
  CHAIN=/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/ShiftBack.chain \
  REJECT=vcf/${sampleName}/${sampleName}_rejected.vcf

java -Xmx10G -jar ${picard} MergeVcfs \
  I=vcf/${sampleName}/${sampleName}_MTshiftedBack_raw.vcf \
  I=vcf/${sampleName}/${sampleName}_MT_raw.vcf \
  O=vcf/${sampleName}/${sampleName}_final_raw.vcf

rm vcf/${sampleName}/${sampleName}_MT_raw.vcf
rm vcf/${sampleName}/${sampleName}_MTshifted8k_raw.vcf
rm vcf/${sampleName}/${sampleName}_MTshiftedBack_raw.vcf
rm vcf/${sampleName}/${sampleName}_MT_raw.vcf.idx
rm vcf/${sampleName}/${sampleName}_MTshifted8k_raw.vcf.idx
rm vcf/${sampleName}/${sampleName}_MTshiftedBack_raw.vcf.idx
# no variants should be rejected
if [ `grep -v \# vcf/${sampleName}/${sampleName}_rejected.vcf | wc -l` == 0 ]
then
  rm vcf/${sampleName}/${sampleName}_rejected.vcf
fi


####################
# MergeMutectStats #
####################
${gatk} --java-options "-Xmx10G -Djava.io.tmpdir=${tmpDir}" MergeMutectStats \
  --stats vcf/${sampleName}/${sampleName}_MT_raw.vcf.stats \
  --stats vcf/${sampleName}/${sampleName}_MTshifted8k_raw.vcf.stats \
  -O vcf/${sampleName}/${sampleName}_final_raw.vcf.stats

rm vcf/${sampleName}/${sampleName}_MT_raw.vcf.stats
rm vcf/${sampleName}/${sampleName}_MTshifted8k_raw.vcf.stats


###############
# FilterCalls #
###############
# autosomal_coverage: "Median coverage of the autosomes for filtering potential polymorphic NuMT variants"
# vaf_filter_thershold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
# f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
# contamination-estimate: "Level of the minor haplogroup estimated by haplochecker."

metricsSample=`awk -F "," '{if (NR=='$(( ${SGE_TASK_ID} + 1 ))') print $1}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/mtMetrics.csv`
if [ ${metricsSample} == ${sampleName} ]
then
  contamination=`awk -F "," '{if (NR=='$(( ${SGE_TASK_ID} + 1 ))') print $6}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/mtMetrics.csv`
  autosomalCov=`awk -F "," '{if (NR=='$(( ${SGE_TASK_ID} + 1 ))') print $7}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/scripts/mtMetrics.csv`
else
  exit
fi

${gatk} --java-options "-Xmx10G -Djava.io.tmpdir=${tmpDir}" FilterMutectCalls \
  -V vcf/${sampleName}/${sampleName}_final_raw.vcf \
  -R /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.fasta \
  -O vcf/${sampleName}/${sampleName}_final_filtered.vcf \
  --stats vcf/${sampleName}/${sampleName}_final_raw.vcf.stats \
  --max-alt-allele-count 4 \
  --mitochondria-mode \
  --autosomal-coverage ${autosomalCov} \
  --min-allele-fraction 0.01 \
  --f-score-beta 1 \
  --contamination-estimate ${contamination}

${gatk} --java-options "-Xmx10G -Djava.io.tmpdir=${tmpDir}" VariantFiltration \
  -V vcf/${sampleName}/${sampleName}_final_filtered.vcf \
  -O vcf/${sampleName}/${sampleName}_final.vcf \
  --mask /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/blacklist_sites.hg38.chrM.bed \
  --mask-name "blacklisted_site"

rm vcf/${sampleName}/${sampleName}_final_raw.vcf
rm vcf/${sampleName}/${sampleName}_final_raw.vcf.idx
rm vcf/${sampleName}/${sampleName}_final_filtered.vcf
rm vcf/${sampleName}/${sampleName}_final_filtered.vcf.idx


###################
# bgzip and index #
###################
bgzip vcf/${sampleName}/${sampleName}_final.vcf
${bcftools} index vcf/${sampleName}/${sampleName}_final.vcf.gz


SGE_TASK_ID=$(( ${SGE_TASK_ID} - ${offset} ))
mv callVariants.e*.${SGE_TASK_ID} logs/callVariants.e.${sampleName}
mv callVariants.o*.${SGE_TASK_ID} logs/callVariants.o.${sampleName}
