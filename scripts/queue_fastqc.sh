#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="batch_fastqc"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=1G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Usage: sbatch queue_fastqc.sh <ERR_list>


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


REF_DIR=/athena/angsd/scratch/jwh4001/project/fastq
OUT_DIR=/athena/angsd/scratch/jwh4001/project/QC/fastqc


# run fastqc for first paired-end read
for ERR in `cat $ERR_FILE`;
    do
    if [ ! -r ${OUT_DIR}/${ERR}_1_fastqc.html ]; then
        sbatch fastqc.sh ${ERR}_1 ${REF_DIR}/${ERR} ${OUT_DIR}
    else
        echo "${OUT_DIR}/${ERR}_1_fastqc.html already exists..."
    fi;
done

# run fastqc for second paired-end read
for ERR in `cat $ERR_FILE`;
    do
    if [ ! -r ${OUT_DIR}/${ERR}_2_fastqc.html ]; then
        sbatch fastqc.sh ${ERR}_2 ${REF_DIR}/${ERR} ${OUT_DIR}
    else
        echo "${OUT_DIR}/${ERR}_2_fastqc.html already exists..."
    fi;
done

exit
