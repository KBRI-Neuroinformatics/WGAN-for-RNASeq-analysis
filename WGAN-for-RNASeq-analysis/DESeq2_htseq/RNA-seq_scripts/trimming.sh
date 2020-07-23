#!/bin/bash

for i in `cat ../FASTQ/fastq1.list`
    do
    #echo $i
    ifp=${i%_*z}
    ifpfz=${ifp}.'fq.gz'

    echo '---', $ifp
    #echo $i, $ifp, $ifpfz


    java -jar /home/neuroinform/bin/Trimmomatic-0.38/trimmomatic-0.38.jar PE -basein ../FASTQ/$i -baseout $ifpfz ILLUMINACLIP:/home/neuroinform/bin/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    #break
    done
