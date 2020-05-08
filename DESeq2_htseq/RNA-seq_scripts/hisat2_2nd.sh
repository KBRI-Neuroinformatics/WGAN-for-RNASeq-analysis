#!/bin/sh
#sh hisat2.sh 2>&1 | tee hisat2.log

FMiD='/storage_data/mkcheon/RNA-seq/Ref/HISAT2/grch38_tran'
UnmD='/storage_data/mkcheon/RNA-seq/GSE104775_5xFAD_TREM2_Cortex_72/Hisat2_GRCm38.93_with_Known_dta-cuff/2nd_Alignments_Unmapped'

for i in `cat ../../Prefetch/SRX.list`
    do
    jc=$i
    echo " * "$jc

    i1=$i'.unmapped_1P.fq.gz'
    i2=$i'.unmapped_2P.fq.gz'

    isam=$i'.2nd.sam'
    ibam=$i'.2nd.bam'
    ibai=$i'.2nd.bai'

    #echo $i1, $i2, $isam, $ibam, $ibai 

    #hisat2
    echo "    -> Hisat2"
    hisat2 -p 8 --dta-cufflinks -x $FMiD/genome_tran -1 $UnmD/$i1 -2 $UnmD/$i2 -S $isam

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
