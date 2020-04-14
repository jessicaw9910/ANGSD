#! / bin / bash

# Usage : get_files.sh <ERR_list> <ftp_1> <ftp_2>
# Check that we have our command line argument (s)

if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments"
    exit
fi

# Read argument from command line

ERR_LIST=$1
FTP_1=$2
FTP_2=$3

for recs in `cat $ERR_LIST`;
    do
    if [[ ! -f ${recs}/${recs}_1.fastq.gz ]]; then
        LINK="$(egrep ${recs} $FTP_1)"
        wget -P `pwd`/${recs}/ ftp://${LINK}
    else
        echo "${recs}/${recs}_1.fastq.gz exists"
    fi
    if [[ ! -f ${recs}/${recs}_2.fastq.gz ]]; then
        LINK="$(egrep ${recs} $FTP_2)"
        wget -P `pwd`/${recs}/ ftp://${LINK}
    else
        echo "${recs}/${recs}_2.fastq_gz exists"
    fi;
done

exit
