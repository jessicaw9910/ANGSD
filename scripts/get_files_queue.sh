#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="fastq_download"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=2G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=jwh4001@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

# Arguments no longer from command line

ERR_LIST=ERR_list.txt
FTP_1=ftp_1.txt
FTP_2=ftp_2.txt

for recs in `cat $ERR_LIST`;
    do
    if [[ ! -f ${recs}/${recs}_1.fastq.gz ]]; then
        LINK="$(egrep ${recs} $FTP_1)"
        wget -P `pwd`/${recs}/ ftp://${LINK}
    else
        echo "${recs}/${recs}_1.fastq.gz exists"
        LINK="$(egrep ${recs} $FTP_1)"
        wget -P `pwd`/${recs}/ -c ftp://${LINK}
    fi
    if [[ ! -f ${recs}/${recs}_2.fastq.gz ]]; then
        LINK="$(egrep ${recs} $FTP_2)"
        wget -P `pwd`/${recs}/ ftp://${LINK}
    else
        echo "${recs}/${recs}_2.fastq_gz exists"
        LINK="$(egrep ${recs} $FTP_2)"
        wget -P `pwd`/${recs}/ -c ftp://${LINK}

    fi;
done

exit
