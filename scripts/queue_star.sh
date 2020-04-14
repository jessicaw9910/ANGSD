#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="batch_STAR"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=1G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Usage: sbatch queue_star.sh <ERR_list>


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


REF_DIR=/athena/angsd/scratch/jwh4001/project/hg38_STARindex/

# set align directory to test for additional cancer samples, remove normally
ALIGN_DIR=/athena/angsd/scratch/jwh4001/project/alignment/test/


for ERR in `cat $ERR_FILE`;
    do
    if [ ! -r ${ALIGN_DIR}/${ERR}.Aligned.sortedByCoord.out.bam ]; then
        sbatch star_updated.sh ${ERR} ${REF_DIR} ${ALIGN_DIR}
    else
        echo "${ALIGN_DIR}/${ERR}.Aligned.sortedByCoord.out.bam already exists..."
    fi;
done

exit
