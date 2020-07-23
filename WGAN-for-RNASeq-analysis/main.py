import argparse
import numpy as np
import preprocess
import model
from model import WGAN_GP
import scipy.stats
import pandas as pd
import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0"

parser = argparse.ArgumentParser()
parser.add_argument('--refpath', type=str, default='83.93_reference.tab', help='path of the reference')
parser.add_argument('--genepath', type=str, default='DESeq2_htseq/AD7M_over_WT7M_q005_baseMean100_PMean.csv', help='path of DEG results(for gene filtering)')
parser.add_argument('--datapath', type=str, default='DESeq2_htseq/GSE104775_RLD_w_HG.csv', help='path of DEG results')
args = parser.parse_args()

#1. Preprocess dataset
##Gene filtering and augmentation
ex_gnames, ex_glists = preprocess.gene_filtering(args.refpath, args.genepath)
ex_value, ex_gene, ex_label = preprocess.extract_filtered_gene(args.datapath, ex_gnames)
metadatas = preprocess.get_metadata(ex_label)
egene, rlds, aug_ages, aug_types, sp_size = preprocess.data_augmentation(ex_value, ex_gene, metadatas)
#egene, rlds, aug_ages, aug_types, sp_size = preprocess.gaussian_augmentation(ex_value, ex_gene, metadatas, 141)
#print (egene.shape, rlds.shape, aug_ages.shape, aug_types.shape)

##Rescaling(normalization)
rcond_list, augcond_list = preprocess.indexing(aug_types, np.int(rlds.shape[0]/sp_size), sp_size)
re_rld, rld_mean, max_rld_std, std_re_rld = preprocess.rescaling(rlds, augcond_list)

##Split data
xtr, xte = preprocess.data_split('data_split.npz', re_rld, test_ratio=0.1)
#xtr, xte = preprocess.data_split('gaussian_data_split.npz', re_rld, test_ratio=0.1)
n_tr, n_te = len(xtr), len(xte)
print ('xtr.shape:', xtr.shape, 'xte.shape:', xte.shape)

#2. Train WGAN+GP
wgangp = WGAN_GP(xtr, n_te)
dloss, gloss, genx = wgangp.train()
