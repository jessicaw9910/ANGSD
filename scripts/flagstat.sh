#!/bin/bash

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="flagstat"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=5G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

# Usage : flagstat.sh <ERR> <OUT_DIR>

# Check that we have our command line arguments
if [ "$#" -ne 2 ]; then
    echo "You must enter exactly 2 command line arguments: <ERR> <OUT_DIR>"
    exit
fi

# Read argument from command line
ERR=$1
OUT_DIR=$2

spack load samtools@1.9%gcc@6.3.0

samtools flagstat ../alignment/${ERR}.Aligned.sortedByCoord.out.bam > ${OUT_DIR}/${ERR}.flagstat.txt

echo "${OUT_DIR}/${ERR}.flagstat.txt created..."

exit
