import pandas as pd
from os import listdir
from os.path import isfile, join

path = "/storage_data/mkcheon/RNA-seq/GSE104775_5xFAD_TREM2_Cortex_72/Hisat2_GRCm38.93_with_Known_dta-cuff/htseq/"

textfiles = [f for f in listdir(path) if isfile(join(path, f)) and 'raw_counts.txt' in f ]
textfiles = sorted(textfiles)
print(len(textfiles))

Rawcount = pd.read_csv(path+"SRX3266465_raw_counts.txt", sep='\t', header=None)
countdf = pd.DataFrame()
countdf['gene'] = Rawcount[0]

ar = []
for f in textfiles :
    Rawcount = pd.read_csv(path+f, sep='\t', header=None)
    #ar.extend(list(Rawcount[0]))
    countdf[f.split('_')[0]] = Rawcount[1]
print(countdf.head())

countdf.to_csv(path+"Raw_counts_from_HTSeq.csv", sep=',', encoding='utf-8', index=False)


