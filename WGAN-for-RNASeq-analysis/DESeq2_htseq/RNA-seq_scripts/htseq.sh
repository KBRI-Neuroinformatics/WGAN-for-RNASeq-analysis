#!/bin/sh
#sh htseq.sh 2>&1 | tee htseq.log

FMiD='/storage_data/mkcheon/RNA-seq/Ref/HISAT2/grcm38_tran'

#for i in `cat ../../Prefetch/SRX.list`
for i in `cat SRX-1.list`
    do
    echo $i

    ibam=$i'.bam'

    #htseq-count -q -r pos -f bam ../$ibam $FMiD/Mus_musculus.GRCm38.93.gtf > $i'_raw_counts.txt'

    htseq-count -s no -f bam ../$ibam $FMiD/Mus_musculus.GRCm38.93.gtf > $i'_raw_counts.txt'

#   break
    done



