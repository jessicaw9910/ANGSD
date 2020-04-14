#! / bin / bash

# Usage : star.sh <fastq_file_1> <fastq_file_2> <genref_dir> <alignment_dir>
# Check that we have our command line argument (s)

if [ "$#" -ne 4 ]; then
    echo "You must enter exactly 2 command line arguments"
    exit
fi

# Read argument from command line
# check that necessary files and directories exist

FILE1=$1
FILE2=$2
REF_DIR=$3
ALIGN_DIR=$4

if [[ ! -f $FILE1 ]]; then 
    echo "$FILE1 does not exist . Exiting ... ";
    exit
else
    echo "$FILE1 exists ..."
fi


if [[ ! -f $FILE2 ]]; then
    echo "$FILE2 does not exist . Exiting ... ";
    exit
else
    echo "$FILE2 exists ..."
fi


#check to see if reference genome directory exist
if [[ ! -d $REF_DIR ]]; then
    echo "Reference genome directory does not exist . Exiting ... "
    exit
else
    echo "Reference genome directory exists ..."
fi


#check to see if alignment output directory exist
if [[ ! -d $ALIGN_DIR ]]; then
    echo "Alignment output directory does not exist . Exiting ... "
    exit
else
    echo "Alignment output directory exists ..."
fi


# Load packages that we will need
spack load star@2.7.0e
spack load samtools@1.9% gcc@6.3.0

f1="$(basename -- $FILE1)"
F1="${f1%.fastq.gz}"

f2="$(basename -- $FILE2)"
F2="${f2%.fastq.gz}"

FCOMBO="${F1%_1}"


# Run STAR , if result not already present
if [ ! -r ${WORKING_DIRECTORY}/${F1}.Aligned.sortedByCoord.out.bam ]; then
    STAR --runMode alignReads \
         --genomeDir ${REF_DIR} \
         --readFilesIn ${FILE1} ${FILE2} \
         --readFilesCommand zcat \
         --alignIntronMin 10 \
	 --outFilterMultimapNmax 20 \
	 --alignSJoverhangMin 8 \
         --outFileNamePrefix ${ALIGN_DIR}/${FCOMBO}. \
         --outSAMtype BAM SortedByCoordinate
fi


# Samtools index on BAM file
if [ -r ${ALIGN_DIR}/${F1}.Aligned.sortedByCoord.out.bam ]; then
    samtools index ${ALIGN_DIR}/${F1}.Aligned.sortedByCoord.out.bam
fi

if [ -r ${ALIGN_DIR}/${F2}.Aligned.sortedByCoord.out.bam ]; then
    samtools index ${ALIGN_DIR}/${F2}.Aligned.sortedByCoord.out.bam
fi

exit
