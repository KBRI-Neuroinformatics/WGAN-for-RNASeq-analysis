#!/bin/sh
#sh hisat2.sh 2>&1 | tee hisat2.log

FMiD='/storage_data/mkcheon/RNA-seq/Ref/HISAT2/grcm38_tran'
TrimD='/storage_data/mkcheon/RNA-seq/GSE104775_5xFAD_TREM2_Cortex_72/Trim'

for i in `cat ../Prefetch/SRX.list`
    do
    jc=$i
    echo " * "$jc

    i1=$i'_1P.fq.gz'
    i2=$i'_2P.fq.gz'
    isam=$i'.sam'
    ibam=$i'.bam'
    ibai=$i'.bai'

    #echo $i1, $i2, $isam, $ibam, $ibai 

    #hisat2
    echo "    -> Hisat2"
    hisat2 -p 8 --dta-cufflinks -x $FMiD/genome_tran -1 $TrimD/$i1 -2 $TrimD/$i2 -S $isam

    #samtools
    echo "    -> Samtools_sort"    
    samtools sort -@ 8 -o $ibam $isam

    #samtools
    echo "    -> Samtools_index"    
    samtools index -b $ibam $ibai

    #remove .sam
    rm $isam

    #break
    done
