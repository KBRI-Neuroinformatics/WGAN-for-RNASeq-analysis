import numpy as np
import pandas as pd

def gene_filtering(refpath, genepath):    
    #get enid and gid from refpath
    dat = np.genfromtxt(refpath, delimiter='\t', dtype='str')[1:]
    ref_enid, ref_gid = dat[:,0], dat[:,1]

    #get genelist from genepath
    gnames = np.array(pd.read_csv(genepath))[:,0]
    #extract mouse genes
    gnames = np.intersect1d(ref_enid, gnames)
    
    #get genelists
    glists = []
    for g in range(len(gnames)):
        idx = np.where(gnames[g]==ref_enid)[0][0]
        glists.append(ref_gid[idx])
    glists = np.array(glists)
    return gnames, glists

def extract_filtered_gene(datapath, gnames):
    #Load dataset
    rld = np.genfromtxt(datapath, delimiter=',', dtype='str')
    label, gene, value = rld[0][1:], rld[1:,0], rld[1:,1:].T

    label = np.array([x.replace('"','') for x in label])
    gene = np.array([x.replace('"','') for x in gene])
    value = value.astype(np.float32)

    #Extract filtered genes
    eidx = []
    for g in range(len(gnames)):
        idx = np.where(gnames[g]==gene)[0][0]
        eidx.append(idx)

    evalue, egene = value[:, eidx], gnames
    print ('# filtered genes:', egene.shape,  ', data.shape :', evalue.shape)
    return evalue, egene, label

def get_metadata(label): #Fit your data 
    ages, genders, types = [], [], []
    for l in range(len(label)):
        sidx = 2 #AD, TA, WT
        if len(label[l])==5: #T
            sidx = 1
        types.append(label[l][:sidx])
        ages.append(label[l][sidx:sidx+2])
        genders.append(label[l][sidx+2:sidx+3])
    ages, genders, types = np.array(ages), np.array(genders), np.array(types)
    regions = np.array(['CX']*len(ages))
    return [ages, genders, types, regions]

def data_augmentation(evalue, egene, metainfo, save=True):
    ages, genders, types, regions = metainfo
    ad_idx, wt_idx = np.where(types=='AD')[0], np.where(types=='WT')[0]

    ext_idx = np.array([ad_idx, wt_idx])
    ext_values = np.array([evalue[ad_idx], evalue[wt_idx]])
    ext_ages = np.array([ages[ad_idx], ages[wt_idx]])
    ext_regions = np.array([regions[ad_idx], regions[wt_idx]])
    ext_genders = np.array([genders[ad_idx], genders[wt_idx]])
    ext_types = np.array([types[ad_idx], types[wt_idx]])

    aug_values = []
    aug_ages, aug_types = [], []

    ranges = (np.arange(9)+1)*.1
    sp_size = np.int(len(ext_ages[0])/len(np.unique(ages)))
    
    for t in range(len(ext_idx)): #types
        #print ('types:', t)
        tmpx = ext_values[t]
        tmpa, tmpt = ext_ages[t], ext_types[t]

        for a in range(len(np.unique(ages))): #ages
            #print ('  ages:', a)
            tmp_value = tmpx[(a*sp_size):(a*sp_size)+sp_size]
            tmp_ages = tmpa[(a*sp_size):(a*sp_size)+sp_size]
            tmp_types = tmpt[(a*sp_size):(a*sp_size)+sp_size]

            for s in range(sp_size): #samples
                for s_ in range(s, sp_size):
                    tmp1, tmp2 = tmp_value[s], tmp_value[s_]
                    if(s==s_): #rawdata
                        aug_values.append(tmp1)
                        aug_ages.append(tmp_ages[s])
                        aug_types.append('org'+tmp_types[s])
                    else: #augmentation
                        for r in ranges:
                            augx = ((tmp1*r)+(tmp2*(1.-r)))
                            aug_values.append(augx)
                            aug_ages.append(tmp_ages[s])
                            aug_types.append(tmp_types[s])

    aug_values, aug_ages, aug_types = np.array(aug_values), np.array(aug_ages), np.array(aug_types)
    aug_genders, aug_regions = np.concatenate(ext_genders), np.concatenate(ext_regions)
    print ('augmented_data.shape :', aug_values.shape)
    
    #save
    if save==True:
        np.savez('augmented_input.npz', genelist=egene, values=aug_values, ages=aug_ages, types=aug_types, genders=aug_genders, regions=aug_regions)
        print ('Augmented input is saved!')
    
    return [egene, aug_values, aug_ages, aug_types, sp_size]

def gaussian_augmentation(evalue, egene, metainfo, naug, save=True):
    ages, genders, types, regions = metainfo
    ad_idx, wt_idx = np.where(types=='AD')[0], np.where(types=='WT')[0]

    ext_idx = np.array([ad_idx, wt_idx])
    ext_values = np.array([evalue[ad_idx], evalue[wt_idx]])
    ext_ages = np.array([ages[ad_idx], ages[wt_idx]])
    ext_regions = np.array([regions[ad_idx], regions[wt_idx]])
    ext_genders = np.array([genders[ad_idx], genders[wt_idx]])
    ext_types = np.array([types[ad_idx], types[wt_idx]])

    aug_values = []
    aug_ages, aug_types = [], []

    ranges = (np.arange(9)+1)*.1
    sp_size = np.int(len(ext_ages[0])/len(np.unique(ages)))
    
    for t in range(len(ext_idx)): #types
        tmpx = ext_values[t]
        tmpa, tmpt = ext_ages[t], ext_types[t]
        
        for a in range(len(np.unique(ages))): #ages
            tmp_value = tmpx[(a*sp_size):(a*sp_size)+sp_size]
            tmp_ages = tmpa[(a*sp_size):(a*sp_size)+sp_size]
            tmp_types = tmpt[(a*sp_size):(a*sp_size)+sp_size]

            mean, cov = np.average(tmp_value, axis=0), np.cov(tmp_value, rowvar=False)
            
            augx = np.random.multivariate_normal(mean, cov, naug-sp_size)
            auga = np.array([tmp_ages[0] for _ in range(naug-sp_size)])
            augt = np.array([tmp_types[0] for _ in range(naug-sp_size)])
            
            aug_values.append(tmp_value)
            aug_ages.append(tmp_ages)
            aug_types.append(np.char.add(['org']*sp_size, tmp_types))
            aug_values.append(augx)
            aug_ages.append(auga)
            aug_types.append(augt)
            
    aug_values, aug_ages, aug_types = np.concatenate(aug_values, axis=0), np.concatenate(aug_ages, axis=0), np.concatenate(aug_types, axis=0)
    aug_genders, aug_regions = np.concatenate(ext_genders), np.concatenate(ext_regions)
    print ('gaussian augmented_data.shape :', aug_values.shape)
    
    #save
    if save==True:
        np.savez('gaussian_augmented_input.npz', genelist=egene, values=aug_values, ages=aug_ages, types=aug_types, genders=aug_genders, regions=aug_regions)
        print ('Gaussian augmented input is saved!')
    
    return [egene, aug_values, aug_ages, aug_types, sp_size]

def indexing(types, naug, sp_size):
    idx = [idx for idx, t in enumerate(types) if 'org' in t]
    ridx, augidx = [], []
    for t in range(np.int(len(idx)/sp_size)):
        ridx.append(idx[sp_size*t:sp_size*(t+1)])
        augidx.append(np.arange(naug*t,naug*(t+1)))
    AD2M, AD4M, AD7M = ridx[:3]
    WT2M, WT4M, WT7M = ridx[3:] 
    augAD2M, augAD4M, augAD7M = augidx[:3]
    augWT2M, augWT4M, augWT7M = augidx[3:]
    return [AD2M, AD4M, AD7M, WT2M, WT4M, WT7M], [augAD2M, augAD4M, augAD7M, augWT2M, augWT4M, augWT7M]

def rescaling(values, clist):
    #get mean and max_std
    rld_mean = np.average(values, axis=0)
    rld_std = []
    for cond in range(len(clist)):
        rld_std.append(np.std(values[clist[cond]], axis=0)) #per genes
    max_rld_std = np.max(rld_std, axis=0)

    #rescaling
    re_rld = (values-rld_mean)/max_rld_std
    std_re_rld = np.std(re_rld.flatten())
    re_rld = re_rld/(3.918*std_re_rld)+0.5
    return re_rld, rld_mean, max_rld_std, std_re_rld

def data_split(sname, data, test_ratio=0.1):
    ndata = data.flatten()
    
    xidx = np.arange(len(data))
    teidx = np.random.choice(len(data), int(len(data)*test_ratio), False)
    tridx = np.setdiff1d(xidx, teidx)
    
    dat_tr, dat_te = data[tridx], data[teidx]
    np.savez(sname, tridx=tridx, teidx=teidx)
    print ('Save data :', sname)
    return dat_tr, dat_te