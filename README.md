# ANGSD Final Project

## Effects of Smoking on Gene Expression in South Korean Patients with Lung Adenocarcinoma
This project analyzes data from 36 lung adenocarcinoma patients undergoing lung resections in South Korea.  A total of 63 samples were collected from these 36 patients, including 27 paired-normal biopsies and 36 tumor biopsies.

The original publication referenced is [Seo et al., *Genome Research*, 2012.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3483540/pdf/2109.pdf).  The data are available [here](https://www.ebi.ac.uk/ena/data/view/PRJEB2784) in the ENA.

## counts

Contains interim and final read count files assembled via featureCounts

+ **project_20.03.14.txt**
+ **project_20.03.14.txt.summary**
+ **project_20.04.04.txt**
+ **project_20.04.04.txt.summary**
+ **project_20.04.07.txt**
+ **project_20.04.09.txt** - final count used in downstream analysis

## files

Contains files needed to run scripts and key outputs

+ **Backup_v01.csv** - input needed for ANGSD_webscraping.ipynb taken from [GEO](www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419)
+ **Backup_v02.csv** - output of ANGSD_webscraping.ipynb
+ **Backup_v02.xlsm** - modified to provide additional information
+ **col_names.txt** - needed to run Rmd analysis
+ **demographics.txt** - age at diagnosis, smoking status, gender, and stage; used in Rmd
+ **experiment_list.txt** - needed for file downloading scripts
+ **script_table.csv** - table with description of scripts and links to GitHub repo for Rmd
+ **sig_genes.csv** - output of Rmd DESeq2 analysis
+ **smoking_status.txt** - smoking status annotation file; not currently in use

## QC

Contains quality control outputs

+ **bamqc/** - folder containing BamQC reports for all 63 sample alignments
+ **multiqc_report_2020.03.15.html** - MultiQC report as of 3/15
+ **multiqc_report_2020.04.12 (PE1).html** - MultiQC report as of 4/12 for paired-end 1 reads
+ **multiqc_report_2020.04.12 (PE2).html** - MultiQC report as of 4/12 for paired-end 2 reads

## scripts

Contains scripts needed to process data, including:

+ **ANGSD_webscraping.ipynb** - allows collection of per sample smoking and cancer vs. paired-normal data
+ **fastqc.sh** - runs FastQC on a single FASTQ file
+ **flagstat.sh** - runs Flagstat on a single BAM file
+ **geneBody.sh** - runs RSeQC geneBody_coverage (took too long to run to incorporate into output)
+ **get_files.sh** - an sbatch script to download all FASTQ files from the SRA (doesn't check for continuation)
+ **get_files_queue** - an sbatch script to download all FASTQ files from the SRA (does check for continuation)
+ **queue_fastqc.sh** - queues the FastQC jobs individually from ERR_list
+ **queue_flagstat.sh** - queues the Flagstat jobs individually from ERR_list
+ **queue_star.sh** - takes ERR_list and creates individual sbatch commands to run STAR alignment
+ **queue_geneBody.sh** - queues geneBody_coverage jobs individually from ERR_list
+ **queue_rseqc.sh** - queues read_distribution jobs individually from ERR_list
+ **queue_star.sh** - queues STAR jobs individually from ERR_list
+ **read_counts.sh** - runs FeatureCounts (could not get this to work)
+ **rseqc.sh** - runs RSeQC read_distribution 
+ **star.sh** - basic 1-pass STAR script
+ **star_for_loop.sh** - runs STAR 1-pass as for loop
+ **star_updated.sh** - run STAR alignment using basic 2-pass

## unmapped

Contains unmapped reads from cancer sample ERR164585 and BLAST alignment output

+ **BLAST_unmapped_1_subset.txt** - BLAST alignment output for last 10 unmapped reads from paired-end 1
+ **BLAST_unmapped_2_subset.txt** - BLAST alignment output for last 10 unmapped reads from paired-end 1
+ **ERR164585.Unmapped.out.mate1_subset.fasta** - ouput of last 10 unmapped reads in FASTA format from paired-end 1
+ **ERR164585.Unmapped.out.mate2_subset.fasta** - ouput of last 10 unmapped reads in FASTA format from paired-end 2