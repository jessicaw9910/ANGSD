# ANGSD Final Project

## counts

Contains read count files assembled over time

## files

Contains key files needed to run scripts

## QC

Contains quality control outputs

## scripts

Contains scripts needed to process data including:
+ **get_files_queue** - an sbatch script to download all FASTQ files from the SRA
+ **queue_star** - takes ERR_list and creates individual sbatch commands to run STAR alignment
+ **star_updated** - run STAR alignment using basic 2-pass