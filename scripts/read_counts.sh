#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="read_counts"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=1G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue


# Usage: read_counts

var=$(date +'%y.%m.%d')

spack load subread

featureCounts -p -a ../hg38.99.gtf.gz ../read_counts/project_${var}.txt ../alignment/*.bam

exit
