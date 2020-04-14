#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="STAR_alignment"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=35G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Usage: sbatch star_updated.sh <ERR> <genref_dir> <alignment_dir>


# Check that we have our command line arguments
if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments: <ERR> <generef_dir> <alignment_dir>"
    exit
fi


# Read argument from command line
ERR=$1
REF_DIR=$2
ALIGN_DIR=$3


# check to see if reference genome directory exist
if [[ ! -d $REF_DIR ]]; then
    echo "Reference genome directory does not exist. Exiting... "
    exit
else
    echo "Reference genome directory exists..."
fi


# check to see if alignment output directory exist
if [[ ! -d $ALIGN_DIR ]]; then
    echo "Alignment output directory does not exist. Exiting... "
    exit
else
    echo "Alignment output directory exists..."
fi


# Load packages that we will need
spack load star@2.7.0e
spack load samtools@1.9% gcc@6.3.0


FILE1=/athena/angsd/scratch/jwh4001/project/fastq/${ERR}/${ERR}_1.fastq.gz
FILE2=/athena/angsd/scratch/jwh4001/project/fastq/${ERR}/${ERR}_2.fastq.gz

if [[ ! -f $FILE1 ]]; then 
    echo "$FILE1 does not exist. Exiting...";
    exit
else
    echo "$FILE1 exists..."
fi

if [[ ! -f $FILE2 ]]; then
    echo "$FILE2 does not exist. Exiting...";
    exit
else
    echo "$FILE2 exists..."
fi


# Run 2-pass STAR, if result not already present
# added --outReadsUnmapped to check a couple of cancer ones, can turn off
if [ ! -r ${ALIGN_DIR}/${ERR}.Aligned.sortedByCoord.out.bam ]; then
    STAR --runMode alignReads \
         --runThreadN 4 \
         --genomeDir ${REF_DIR} \
         --readFilesIn ${FILE1} ${FILE2} \
         --readFilesCommand zcat \
         --twopassMode Basic \
         --alignIntronMin 10 \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --outFileNamePrefix ${ALIGN_DIR}/${ERR}. \
	 --outReadsUnmapped Fastx \
         --outSAMtype BAM SortedByCoordinate
fi

# Samtools index on BAM file if index file doesn't already exist
if [ ! -r ${ALIGN_DIR}/${ERR}.Aligned.sortedByCoord.out.bam.bai ]; then
    samtools index ${ALIGN_DIR}/${ERR}.Aligned.sortedByCoord.out.bam
else
    echo "${ALIGN_DIR}/${ERR}.Aligned.sortedByCoord.out.bam.bai already exists..."
fi

exit
