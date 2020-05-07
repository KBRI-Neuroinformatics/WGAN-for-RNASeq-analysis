[![DOI](https://zenodo.org/badge/255605797.svg)](https://zenodo.org/badge/latestdoi/255605797)
## A practical application of generative adversarial networks for RNA-seq analysis to predict the molecular progress of Alzheimer's disease  
-----
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

    * __--refpath__ and __--genepath__ are used to obtain filtered gene information result from DEG analysis
    * __--datapath__ is the path of the merged data for HTSeq results of individual samples, and this file is used as input data
    * After running main.py, you can get two .npz files(data_split.npz and weights/training_results.npz) and weights folder, and these output files will be used in following the analysis part   
    
* __Analysis using WGAN_for_rnaseq_analysis.ipynb__
    ![image](https://user-images.githubusercontent.com/57948381/80566407-e64a1e80-8a2d-11ea-8ca7-a5ffaeb193cf.png)
    * The left figure shows the loss graph for training and the right figure shows average correlation of generated samples during training. High correlation means that the generator can create realistic samples. You can see when the performance of generator is good by looking at the right figure. We gave attention from 75K to 125K as the range of optimization then verified the performance of network in various ways.

    ![image](https://user-images.githubusercontent.com/57948381/80568186-8190c300-8a31-11ea-8da9-3424d7484ecc.png)
    * tSNE plot according to different epochs and each red/blue/orange dot means training/testing/generated samples. The more generator generates realistic samples, the more the distribution of generated samples becomes similar to the distribution of real samples.
    * Just change ep_ary to the array you want to check in the plot.tsne function.
    
    ![image](https://user-images.githubusercontent.com/57948381/80574499-ea316d00-8a3c-11ea-85c9-d24a16cea96b.png)
    * Even if the average correlation of samples generated is high, you can see variations of generated values from 75K to 125K. A detailed analysis of this will be handled later using heatmaps, and we set 100K weights as basis then analyse the remaining parts.
    * Just change kwd, condidx and spidx you want to check. We have six conditions(AD2M/4M/7M and WT2M/4M/7M) and each condition has six samples. So 'condidx=2 and spidx=3' means third sample of AD7M. The spidx starts from 1 but condidx starts from 0. (e.g., to check sixth of WT4M, set condidx and spidx 4 and 6 respectively.)

    ![image](https://user-images.githubusercontent.com/57948381/80571019-d551db00-8a36-11ea-8a10-eea735fc33bd.png)
    * In addition to tSNE analysis, validate the performance of generator again by comparing distribution of sample values and correlation coefficient. Each two distribution between real and generate samples looks similar, that implied this model is properly trained.
    
    ---
    * Remaining parts, we set 100K as basis for more analysis.
    * To set different epoch for analysis, change fWTAD_rz idx. Our study utilizes the results of the analysis.latent_interpolation function from 75K to 125K, so fWTAD_rz[50] means latent variables that generate realistic samples at 100K. (i.e., fWTAD_rz[0] means latent variables that generate a realistic samples at 75K)
    * If the range for analysis.latent_interpolation is different with raw codes, change it flexibly according to your environment.
    ![image](https://user-images.githubusercontent.com/57948381/80583056-ac3b4580-8a4a-11ea-8711-a64298d70128.png)
    
    
    ![image](https://user-images.githubusercontent.com/57948381/80579533-2e287000-8a45-11ea-8b33-309bc815f009.png)
    * Compare the generated samples using a latent variable that generates realistic samples to real samples. The index of condition starts from zero, so the value of 2 means third of six conditions(i.e., AD_7M).
   
    ![image](https://user-images.githubusercontent.com/57948381/80573546-1e0b9300-8a3b-11ea-9363-bfc348bf77f9.png)
    * Compare the generated value to real value of a specific gene, Apoe in this figure. Each left bar represents raw expression value and right bar represents generated expression value using a latent variable that generates a realistic sample. And each scatter point represents the Pearson's correlation coefficient between real and generated samples using the optimal latent variable with strong correlation coefficient for real values. Each scatter point represents the Pearson's correlation coefficient between real and generated samples, and based on these coefficients we obtained optimal latent vectors which can generate realistic samples.
    * Just change kwd to a specific gene you want to check.
    
    ![image](https://user-images.githubusercontent.com/57948381/80600173-f7fbe800-8a66-11ea-8b55-c5417f9bc496.png)
    * Simulate virtual disease progression using average latent variables from 75K to 125K for all augmented 7M samples because 7M samples have more distinct character between two conditions, WT and AD. The line represents simulated gene expression values as Alzheimer's progresses and the shaded area represents standard deviation of simulated values of epochs.
    * These transition curves can be represented as raw scales, and each point is the value of real samples we have.(The green line means 7M samples)
    * Just change flist to a list of genes you want to check. 
    
    ![image](https://user-images.githubusercontent.com/57948381/80604509-c5ed8480-8a6c-11ea-802f-66065bffb399.png)
    * We further performed gene ontology and pathway analysis using WebGestalt with genes used to train this network, that is, filtered genes from DEG analysis on 7M samples. You can use these transition curves to analyze any pathways you are interested in. In addition to the transition curve, we can also represent same values using heatmaps.
    * If you have Webgestalt results, you only need to change upplist and dnplist to a pathway you are interested in, but otherwise you have to create your own code by referencing our code.
    

---
#### Reference code
>https://github.com/luslab/scRNAseq-WGAN-GP/blob/master/scripts/WGAN-GP_minimal.py
