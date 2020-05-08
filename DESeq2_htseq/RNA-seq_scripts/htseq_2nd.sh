#!/bin/sh
#sh htseq.sh 2>&1 | tee htseq.log

FMiD='/fstorage_data/mkcheon/RNA-seq/Ref/HISAT2/grch38_tran'

#for i in `cat ../../../Prefetch/SRX.list`
for i in `cat SRX1.list`
    do
    echo $i

    ibam=$i'.2nd.bam'

    #htseq-count -q -r pos -f bam ../$ibam $FMiD/Homo_sapiens.GRCh38.93.gtf > $i'_raw_counts.txt'

    htseq-count -s no -f bam ../$ibam $FMiD/Homo_sapiens.GRCh38.93.gtf > $i'_raw_counts.txt'

#   break
    done



