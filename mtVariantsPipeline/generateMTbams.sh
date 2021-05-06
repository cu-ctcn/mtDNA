#!/bin/bash

#$-l h_vmem=16G
#$-t 1-1200
#$-tc 60
#$-wd /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2
#$-N generateMtBams

source /mnt/mfs/cluster/bin/HgrcPathSetup.sh

offset=0
SGE_TASK_ID=$(( ${SGE_TASK_ID} + ${offset} ))

sampleName=`awk -F "," '{if (NR=='${SGE_TASK_ID}') print $1}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/bam/bamfiles.csv`
bamFile=`awk -F "," '{if (NR=='${SGE_TASK_ID}') print $2}' /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/bam/bamfiles.csv`

gatk="/mnt/mfs/ctcn/tools/gatk-4.1.2.0/gatk"
picard="/mnt/mfs/ctcn/tools/picard_2.20.3/picard.jar"
tmpDir="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/temp"

mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam
mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}
mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}/bwa
mkdir -p /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}/picard


###################
# SubsetBamToChrM #
###################
${gatk} --java-options "-Xmx10G -Djava.io.tmpdir=${tmpDir}" PrintReads \
  -L MT \
  --read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --read-filter MateUnmappedAndUnmappedReadFilter \
  -I ${bamFile} \
  -O bam/${sampleName}/${sampleName}.bam


#################################
# RevertSam (remove alignments) #
#################################
java -Xmx10G -jar -Djava.io.tmpdir=${tmpDir} -jar $picard RevertSam \
  INPUT=bam/${sampleName}/${sampleName}.bam \
  OUTPUT_BY_READGROUP=false \
  OUTPUT=bam/${sampleName}/${sampleName}_unmapped.bam \
  VALIDATION_STRINGENCY=LENIENT \
  ATTRIBUTE_TO_CLEAR=FT \
  ATTRIBUTE_TO_CLEAR=CO \
  SORT_ORDER=queryname \
  RESTORE_ORIGINAL_QUALITIES=false
rm bam/${sampleName}/${sampleName}.ba?


#####################
# Align reads to MT #
#####################
refFasta="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.fasta"
java -Xmx10G -jar ${picard} SamToFastq \
  INPUT=bam/${sampleName}/${sampleName}_unmapped.bam \
  FASTQ=/dev/stdout \
  INTERLEAVE=true \
  NON_PF=true | \
bwa mem -K 100000000 -p -v 3 -t 2 -Y ${refFasta} /dev/stdin - 2> >(tee /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}/bwa/MT.bwa.stderr.log >&2) | \
java -Xmx10G -jar ${picard} MergeBamAlignment \
  VALIDATION_STRINGENCY=SILENT \
  EXPECTED_ORIENTATIONS=FR \
  ATTRIBUTES_TO_RETAIN=X0 \
  ATTRIBUTES_TO_REMOVE=NM \
  ATTRIBUTES_TO_REMOVE=MD \
  ALIGNED_BAM=/dev/stdin \
  UNMAPPED_BAM=bam/${sampleName}/${sampleName}_unmapped.bam \
  OUTPUT=bam/${sampleName}/${sampleName}_MT_unmarked.bam \
  REFERENCE_SEQUENCE=${refFasta} \
  PAIRED_RUN=true \
  SORT_ORDER="unsorted" \
  IS_BISULFITE_SEQUENCE=false \
  ALIGNED_READS_ONLY=false \
  CLIP_ADAPTERS=false \
  MAX_RECORDS_IN_RAM=2000000 \
  ADD_MATE_CIGAR=true \
  MAX_INSERTIONS_OR_DELETIONS=-1 \
  PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
  PROGRAM_RECORD_ID="bwamem" \
  PROGRAM_GROUP_VERSION="0.7.16a-r1181" \
  PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 2 -Y ${refFasta}" \
  PROGRAM_GROUP_NAME="bwamem" \
  UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
  ALIGNER_PROPER_PAIR_FLAGS=true \
  UNMAP_CONTAMINANT_READS=true \
  ADD_PG_TAG_TO_READS=false

java -Xmx10G -jar ${picard} MarkDuplicates \
  INPUT=bam/${sampleName}/${sampleName}_MT_unmarked.bam \
  OUTPUT=bam/${sampleName}/${sampleName}_MT_marked.bam \
  METRICS_FILE=bam/${sampleName}/picard/${sampleName}_MT.marked_dup_metrics.txt \
  VALIDATION_STRINGENCY=SILENT \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  ASSUME_SORT_ORDER="queryname" \
  CLEAR_DT="false" \
  ADD_PG_TAG_TO_READS=false
rm bam/${sampleName}/${sampleName}_MT_unmarked.bam

java -Xmx10G -jar ${picard} SortSam \
  INPUT=bam/${sampleName}/${sampleName}_MT_marked.bam \
  OUTPUT=bam/${sampleName}/${sampleName}_MT.bam \
  SORT_ORDER="coordinate" \
  CREATE_INDEX=true \
  MAX_RECORDS_IN_RAM=300000
rm bam/${sampleName}/${sampleName}_MT_marked.bam


#############################
# Align reads to shifted MT #
#############################
refFasta="/mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/resources/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
java -Xmx10G -jar ${picard} SamToFastq \
  INPUT=bam/${sampleName}/${sampleName}_unmapped.bam \
  FASTQ=/dev/stdout \
  INTERLEAVE=true \
  NON_PF=true | \
bwa mem -K 100000000 -p -v 3 -t 2 -Y ${refFasta} /dev/stdin - 2> >(tee /mnt/mfs/ctcn/team/hklein/analyses/rosmap_WGS_MT/mtMutect2/bam/${sampleName}/bwa/MTshifted8k.bwa.stderr.log >&2) | \
java -Xmx10G -jar ${picard} MergeBamAlignment \
  VALIDATION_STRINGENCY=SILENT \
  EXPECTED_ORIENTATIONS=FR \
  ATTRIBUTES_TO_RETAIN=X0 \
  ATTRIBUTES_TO_REMOVE=NM \
  ATTRIBUTES_TO_REMOVE=MD \
  ALIGNED_BAM=/dev/stdin \
  UNMAPPED_BAM=bam/${sampleName}/${sampleName}_unmapped.bam \
  OUTPUT=bam/${sampleName}/${sampleName}_MTshifted8k_unmarked.bam \
  REFERENCE_SEQUENCE=${refFasta} \
  PAIRED_RUN=true \
  SORT_ORDER="unsorted" \
  IS_BISULFITE_SEQUENCE=false \
  ALIGNED_READS_ONLY=false \
  CLIP_ADAPTERS=false \
  MAX_RECORDS_IN_RAM=2000000 \
  ADD_MATE_CIGAR=true \
  MAX_INSERTIONS_OR_DELETIONS=-1 \
  PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
  PROGRAM_RECORD_ID="bwamem" \
  PROGRAM_GROUP_VERSION="0.7.16a-r1181" \
  PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 2 -Y ${refFasta}" \
  PROGRAM_GROUP_NAME="bwamem" \
  UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
  ALIGNER_PROPER_PAIR_FLAGS=true \
  UNMAP_CONTAMINANT_READS=true \
  ADD_PG_TAG_TO_READS=false

java -Xmx10G -jar ${picard} MarkDuplicates \
  INPUT=bam/${sampleName}/${sampleName}_MTshifted8k_unmarked.bam \
  OUTPUT=bam/${sampleName}/${sampleName}_MTshifted8k_marked.bam \
  METRICS_FILE=bam/${sampleName}/picard/${sampleName}_MTshifted8k.marked_dup_metrics.txt \
  VALIDATION_STRINGENCY=SILENT \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  ASSUME_SORT_ORDER="queryname" \
  CLEAR_DT="false" \
  ADD_PG_TAG_TO_READS=false
rm bam/${sampleName}/${sampleName}_MTshifted8k_unmarked.bam

java -Xmx10G -jar ${picard} SortSam \
  INPUT=bam/${sampleName}/${sampleName}_MTshifted8k_marked.bam \
  OUTPUT=bam/${sampleName}/${sampleName}_MTshifted8k.bam \
  SORT_ORDER="coordinate" \
  CREATE_INDEX=true \
  MAX_RECORDS_IN_RAM=300000
rm bam/${sampleName}/${sampleName}_MTshifted8k_marked.bam

rm bam/${sampleName}/${sampleName}_unmapped.bam



SGE_TASK_ID=$(( ${SGE_TASK_ID} - ${offset} ))
mv generateMtBams.e*.${SGE_TASK_ID} logs/generateMtBams.e.${sampleName}
mv generateMtBams.o*.${SGE_TASK_ID} logs/generateMtBams.o.${sampleName}
