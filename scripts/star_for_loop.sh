#!/bin/bash

# Usage: star.sh <ERR_file> <genref_dir> <alignment_dir_1> <alignment_dir_2>


# Check that we have our command line arguments
if [ "$#" -ne 4 ]; then
    echo "You must enter exactly 4 command line arguments: <ERR_file> <generef_dir> <alignment_dir_1> <alignment_dir_2>"
    exit
fi


# Read argument from command line
ERR_FILE=$1
REF_DIR=$2
ALIGN_DIR_1=$3
ALIGN_DIR_2=$4


# check that necessary files and directories exist
# check to see if ERR file exist
if [[ ! -r $ERR_FILE ]]; then
    echo "ERR file does not exist. Exiting..."
    exit
else
    echo "ERR file exists..."
fi

# check to see if reference genome directory exist
if [[ ! -d $REF_DIR ]]; then
    echo "Reference genome directory does not exist. Exiting... "
    exit
else
    echo "Reference genome directory exists..."
fi

# check to see if alignment output directory 1 exist
if [[ ! -d $ALIGN_DIR_1 ]]; then
    echo "Alignment output directory 1 does not exist. Exiting... "
    exit
else
    echo "Alignment output directory 1 exists..."
fi

# check to see if alignment output directory 2 exist
if [[ ! -d $ALIGN_DIR_2 ]]; then
    echo "Alignment output directory 2 does not exist. Exiting... "
    exit
else
    echo "Alignment output directory 2 exists..."
fi


# Load packages that we will need
spack load star@2.7.0e
spack load samtools@1.9% gcc@6.3.0


# for loop iterates over ERR list and runs 2-pass STAR
for recs in `cat $ERR_FILE`;
    do
    FILE1=fastq/${recs}/${recs}_1.fastq.gz
    FILE2=fastq/${recs}/${recs}_2.fastq.gz    

    if [[ ! -f $FILE1 ]]; then 
        echo "$FILE1 does not exist. Continuing to next entry...";
        continue
    else
        echo "$FILE1 exists..."
    fi

    if [[ ! -f $FILE2 ]]; then
        echo "$FILE2 does not exist. Continuing to next entry...";
        continue
    else
        echo "$FILE2 exists..."
    fi

    # Run 2-pass STAR, if result not already present
    if [ ! -r ${ALIGN_DIR_2}/${recs}.Aligned.sortedByCoord.out.bam ]; then
        #run first pass with files going to $ALIGN_DIR_1 and no SAM/BAM output
        STAR --runMode alignReads \
             --runThreadN 4 \
             --genomeDir ${REF_DIR} \
             --readFilesIn ${FILE1} ${FILE2} \
             --readFilesCommand zcat \
             --alignIntronMin 10 \
 	     --outFilterMultimapNmax 20 \
 	     --alignSJoverhangMin 8 \
             --outFileNamePrefix ${ALIGN_DIR_1}/${recs}. \
             --outSAMtype None
        #run second pass with files going to $ALIGN_DIR_2 using SJ.out.tab file from $ALIGN_DIR_1
        STAR --runMode alignReads \
             --runThreadN 4 \
             --genomeDir ${REF_DIR} \
             --sjdbFileChrStartEnd ${ALIGN_DIR_1}/${recs}.SJ.out.tab \
             --readFilesIn ${FILE1} ${FILE2} \
             --readFilesCommand zcat \
             --alignIntronMin 10 \
	     --outFilterMultimapNmax 20 \
	     --alignSJoverhangMin 8 \
             --outFileNamePrefix ${ALIGN_DIR_2}/${recs}. \
             --outSAMtype BAM SortedByCoordinate
    fi

    # Samtools index on BAM file
    if [ -r ${ALIGN_DIR_2}/${recs}.Aligned.sortedByCoord.out.bam ]; then
        samtools index ${ALIGN_DIR_2}/${recs}.Aligned.sortedByCoord.out.bam
    else
        echo "${ALIGN_DIR_2}/${recs}.Aligned.sortedByCoord.out.bam does not exist..."
    fi;

done

exit
