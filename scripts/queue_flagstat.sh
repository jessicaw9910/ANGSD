#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="batch_flagstat"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=1G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Usage: sbatch queue_flagstat.sh <ERR_list>


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


OUT_DIR=/athena/angsd/scratch/jwh4001/project/QC/flagstat


# run fastqc for first paired-end read
for ERR in `cat $ERR_FILE`;
    do
    if [ ! -r ${OUT_DIR}/${ERR}.flagstat.txt ]; then
        sbatch flagstat.sh ${ERR} ${OUT_DIR}
    else
        echo "${OUT_DIR}/${ERR}.flagstat.txt already exists..."
    fi;
done

exit
