#!/bin/bash

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="fastcq"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=25G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

# Usage : fastqc.sh <ERR> <REF_DIR> <OUT_DIR>

# Check that we have our command line arguments
if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments: <ERR> <REF_DIR> <OUT_DIR>"
    exit
fi

# Read argument from command line
ERR=$1
REF_DIR=$2
OUT_DIR=$3

spack load fastqc

fastqc -o ${OUT_DIR} ${REF_DIR}/${ERR}.fastq.gz

exit
