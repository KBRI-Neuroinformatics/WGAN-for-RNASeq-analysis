
[![DOI](https://zenodo.org/badge/255605797.svg)](https://zenodo.org/badge/latestdoi/255605797)

## A practical application of generative adversarial networks for RNA-seq analysis to predict the molecular progress of Alzheimer's disease


---
### [Prerequisites]
* __Create virtual environment__  

          virtualenv WGAN_rnaseq -p python3
          
          source WGAN_rnaseq/bin/activate  

* __Install packages__  
     * Default:CPU
     * If you want to install tensorflow with GPU support, please visit tensorflow.org to set up your environment. Then replace tensorflow with tensorflow-gpu in requirements.txt

               pip install -r requirements.txt

### [Usage]
* __Training networks__
     
          python main.py --refpath 83.93_reference.tab --genepath DESeq2_htseq/AD7M_over_WT7M_q005_baseMean100_PMean.csv' --datapath DESeq2_htseq/GSE104775_RLD_w_HG.csv  

* __Analysis using WGAN_for_rnaseq_analysis.ipynb__
---
#### Reference code
>https://github.com/luslab/scRNAseq-WGAN-GP/blob/master/scripts/WGAN-GP_minimal.py
