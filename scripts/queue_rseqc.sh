#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="batch_rseqc"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=1G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Usage: sbatch queue_rseqc.sh <ERR_list>


# Check that we have our command line arguments
if [ "$#" -ne 1 ]; then
    echo "You must enter exactly 1 command line arguments: <ERR_list>"
    exit
fi


# Read argument from command line
ERR_FILE=$1

# check to see if ERR list exist
if [[ ! -r $ERR_FILE ]]; then
    echo "ERR list does not exist. Exiting... "
    exit
else
    echo "ERR list exists..."
fi


BED_FILE="/home/jwh4001/angsd/jwh4001/project/hg38_RefSeq.bed"
BAM_DIR="/home/jwh4001/angsd/jwh4001/project/alignment"
OUT_DIR="/home/jwh4001/angsd/jwh4001/project/QC/rseqc"


# run rseqc if file doesn't already exist
for ERR in `cat $ERR_FILE`;
    do
    if [ ! -r ${OUT_DIR}/${ERR}.read_distribution.txt ]; then
        sbatch rseqc.sh ${BED_FILE} ${BAM_DIR}/${ERR}.Aligned.sortedByCoord.out.bam \
        ${OUT_DIR}/${ERR}.read_distribution.txt
    else
        echo "${OUT_DIR}/${ERR}.read_distribution.txt already exists..."
    fi;
done

exit
