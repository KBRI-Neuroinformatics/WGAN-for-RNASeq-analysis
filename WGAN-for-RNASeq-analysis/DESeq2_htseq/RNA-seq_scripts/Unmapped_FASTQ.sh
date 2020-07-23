#!/bin/sh

FMiD='/storage_data/mkcheon/RNA-seq/Ref/HISAT2/grcm38_tran'
TrimD='/storage_data/mkcheon/RNA-seq/GSE104775_5xFAD_TREM2_Cortex_72/Trim'

for i in `cat ../../Prefetch/SRX.list`
#for i in `cat SRX1.list`
    do
    jc=$i
    echo " * "$jc

    i1=$i'.unmapped_1P.fq.gz'
    i2=$i'.unmapped_2P.fq.gz'
    ibam=$i'.bam'
    iunm=$i'.unmapped.bam'

    echo "Unmapped reads"
    samtools view -b -f 4 ../$ibam > $iunm

    echo "Bam to FASTQ"
    samtools fastq -1 $i1 -2 $i2 $iunm

    #break
    done
