#!/bin/bash

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="rseqc"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=25G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

# Usage : fastqc.sh <BED_FILE> <BAM_FILE> <OUT_FILE>

# Check that we have our command line arguments
if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments: <BED_FILE> <BAM_FILE> <OUT_FILE>"
    exit
fi

# Read argument from command line
BED_FILE=$1
BAM_FILE=$2
OUT_FILE=$3


RSEQC_IMAGE="/athena/angsd/scratch/simg/rseqc-3.0.1.simg"

cd ~

spack load singularity@2.6.0

singularity exec $RSEQC_IMAGE read_distribution.py -r $BED_FILE -i $BAM_FILE >> $OUT_FILE


exit
