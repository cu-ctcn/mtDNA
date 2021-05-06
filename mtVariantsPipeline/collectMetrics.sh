#!/bin/bash

#$-l h_vmem=16G
#$-t 1-1200
#$-tc 60
#$-wd /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2
#$-N collectMetrics

source /mnt/mfs/cluster/bin/HgrcPathSetup.sh

offset=0
SGE_TASK_ID=$(( ${SGE_TASK_ID} + ${offset} ))

sampleName=`awk -F "," '{if (NR=='${SGE_TASK_ID}') print $1}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/bam/bamfiles.csv`

picard="/mnt/mfs/ctcn/tools/picard_2.20.3/picard.jar"
mitolib="/mnt/mfs/ctcn/tools/mitolib-0.1.2/mitolib-0.1.2.jar"
tmpDir="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mitoAnalyzer/temp"

mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}/picard
mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}/haplochecker


#####################
# CollectWgsMetrics #
#####################
java -Xmx8G -Djava.io.tmpdir=${tmpdir} -jar $picard CollectWgsMetrics \
  INPUT=bam/${sampleName}/${sampleName}_MT.bam \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.fasta \
  OUTPUT=bam/${sampleName}/picard/${sampleName}_MT.wgs_metrics.txt \
  USE_FAST_ALGORITHM=true \
  READ_LENGTH=151 \
  COVERAGE_CAP=250000 \
  INCLUDE_BQ_HISTOGRAM=true \
  THEORETICAL_SENSITIVITY_OUTPUT=bam/${sampleName}/picard/${sampleName}_MT.theoretical_sensitivity.txt


####################
# GetContamination #
####################
java -Xmx8G -jar ${mitolib} haplochecker \
  --in bam/${sampleName}/${sampleName}_MT.bam \
  --out bam/${sampleName}/haplochecker \
  --ref /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/resources/mtReference/MT.fasta \
  --QUAL 20 \
  --MAPQ 30 \
  --VAF 0.01 > bam/${sampleName}/haplochecker/haplochecker.log


SGE_TASK_ID=$(( ${SGE_TASK_ID} - ${offset} ))
mv collectMetrics.e*.${SGE_TASK_ID} logs/collectMetrics.e.${sampleName}
mv collectMetrics.o*.${SGE_TASK_ID} logs/collectMetrics.o.${sampleName}
