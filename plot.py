import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from sklearn.manifold import TSNE
import seaborn as sns
import scipy.stats
import analysis

sp_size = 6
title_ary = ['AD2M', 'AD4M', 'AD7M', 'WT2M', 'WT4M', 'WT7M']
    
def cond_scatter(dat, cond1, cond2=0):
    if type(cond2)==int:
        cond2=cond1
    plt.figure(figsize=[16,20])
    for y in range(sp_size):
        for x in range(sp_size):
            plt.subplot(sp_size,sp_size,y*sp_size+x+1)
            plt.scatter(dat[cond1[y],:], dat[cond2[x],:])

def loss_graph(dloss, gloss, figsize=0):
    if type(figsize)!=int:  
        plt.figure(figsize=figsize)
    else:
        plt.figure()
    plt.plot(dloss, c='blue', label='discriminator')
    plt.plot(gloss, c='red', label='generator')
    plt.legend(loc=1, fontsize=10, bbox_to_anchor=(1.01,1.01))
    plt.xlabel('Epoch', fontsize=12)
    plt.ylabel('Loss', fontsize=12)
    plt.xticks(np.arange(0, 225000, 25000), ['0', '25K','50K','75K', '100K', '125K', '150K', '175K', '200K'], fontsize=12)
    plt.plot()
        
def avg_corr_graph(corr):
    plt.figure(figsize=(7,5))
    plt.plot(corr, c='blue')
    plt.ylim(0.7, 1.)
    plt.xticks([0,100,200,300,400], ['0','50K','100K','150K','200K'], fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Average correlation for 84 test data', fontsize=12)
    

def tsne(xtr, xte, genx, ep_ary):
    n_tr, n_te = len(xtr), len(xte)
    plt.figure(figsize=(16,3))
    for ep in range(len(ep_ary)):
        tsne = TSNE()
        tmp = np.concatenate((xtr, xte, genx[ep_ary[ep]]), axis=0)
        tmp_tsne = tsne.fit_transform(tmp)
        
        plt.subplot(1,4,ep+1)
        tit = str(500*(ep_ary[ep]+1))+' epoch'
        plt.title(tit)
        plt.scatter(tmp_tsne[:n_tr,0], tmp_tsne[:n_tr,1], alpha=.5, c='red')
        plt.scatter(tmp_tsne[n_tr:(n_tr+n_te),0], tmp_tsne[n_tr:(n_tr+n_te),1], alpha=.5, c='blue')
        plt.scatter(tmp_tsne[(n_tr+n_te):(n_tr+(n_te*2)),0], tmp_tsne[(n_tr+n_te):(n_tr+(n_te*2)),1], alpha=.5, c='orange')
        plt.xlabel('tSNE dimension1', fontsize=13)
        if ep==0:
            plt.ylabel('tSNE dimension2', fontsize=13)
        plt.xlim(-60, 60)
        plt.ylim(-60, 60)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
    plt.subplots_adjust(wspace=0.23)
    plt.legend(['Train', 'Test', 'Generated'], ncol=3, bbox_to_anchor=(-.5, -.3), fontsize=14)
    
def histograms(dat, bins, suptitle, label, figsize=(8,4)): #bins
    col = ['blue', 'red']
    if bins>100:
        xmin, xmax = -.1, 1.1
        yticks = [0,2000,4000,6000,8000,10000]
        lyticks = ['0','2K','4K','6K','8K','10K']
        
    else:
        xmin, xmax = -0.8, 1.
        yticks = [0, 4000, 8000, 12000, 16000]
        lyticks = ['0', '4K', '8K', '12K', '16K']
           
    if type(figsize)!=int:
        plt.figure(figsize=figsize)
    else:
        plt.figure()
        
    plt.suptitle(suptitle, fontsize=15, y=1.03)
    for l in range(len(label)):
        plt.subplot(1,2,(l+1))
        plt.hist(dat[l].flatten(), bins=bins, color=col[l]) 
        plt.xlim(xmin, xmax)
        plt.xticks(fontsize=12)
        if l==0:
            plt.yticks(yticks, lyticks, fontsize=12)
            plt.ylabel('Counts', fontsize=12)
        else:
            plt.yticks(yticks,[])
        plt.xlabel(suptitle[13:], fontsize=12)
        plt.title(label[l], fontsize=12)
    plt.subplots_adjust(wspace=.05)
    
def compare_realx_and_genx(realx, genx, clist, idx):
    plt.figure(figsize=(16,3))
    suptitle = 'Six'+title_ary[idx]+' samples'
    plt.suptitle(suptitle, fontsize=12)
    for i in range(len(clist[idx])):
        plt.subplot(1,6,i+1)
        plt.subplots_adjust(wspace=.3)
        plt.scatter(realx[clist[idx][i]], genx[(6*idx)+i])
        plt.plot([0,1], color='red')
        plt.xlim(-1, 2.4)
        plt.ylim(-1, 2.4)
        plt.yticks(np.arange(-1,3,1), fontsize=10)
        plt.xticks(np.arange(-1,3,1), fontsize=10)
        plt.grid(alpha=.3)
        plt.xlabel('Real_x', fontsize=12)
        if i==0:
            plt.ylabel('Generated_x', fontsize=12)
            
def compare_real_and_gen_gene(realx, genx, gencorr, clists, glists, kwd) :    
    kwdidx = np.where(glists==kwd)[0][0]
    xrange = np.concatenate([np.arange(len(clists))+1]*len(clists))
    clist = np.array(['coral', 'red', 'maroon', 'turquoise', 'royalblue', 'navy'])

    fig = plt.figure(figsize=(16,6))
    ax1 = fig.add_subplot(111)
    #scatter
    ax1.set_ylim(0.8, 1.0)
    ax1.set_yticks(np.arange(0.8, 1.05, 0.05))
    ax1.set_yticklabels(np.arange(0.8, 1.05, 0.05).astype(np.float16), fontsize=12)
    ax1.set_xlim(-0.5, len(gencorr))
    ax1.set_xticks(np.arange(len(gencorr)))
    ax1.set_xticklabels(xrange, fontsize=14)
    ax1.grid(linestyle='--')
    for ns in range(gencorr.shape[0]):
        cls = ns//6
        tx = np.full(gencorr[ns].shape, ns)
        ax1.scatter(tx, gencorr[ns], c=clist[cls], edgecolors='black', linewidth='0.2')
    #bar
    ax2 = ax1.twinx()
    ay2 = ax1.twiny()
    ay2.set_xlim(-0.5, len(gencorr))
    ax2.set_ylim(0., 2.5)
    ax2.set_yticklabels(np.arange(0., 3., 0.5), fontsize=14)
    for ns in range(genx.shape[0]):
        cls = ns//6
        ax2.bar(ns, realx[clists[cls][ns%6],kwdidx], width=0.2, color=clist[cls])
        ax2.bar(ns+.25, genx[ns,kwdidx], width=0.2, color='tan')

    #secondary ticks
    ay2.set_xticks(np.array([-0.5, 2.7, 5.5, 8.7, 11.5, 14.7, 17.5, 20.7, 23.5, 26.7, 29.5, 32.7, 36]))
    ay2.set_xticklabels(['','AD 2M','', 'AD 4M', '', 'AD 7M', '', 'WT 2M', '', 'WT 4M', '', 'WT 7M', ''])
    ay2.xaxis.set_ticks_position('bottom')
    ay2.spines['bottom'].set_position(('outward', 40))

    #label
    ax1.set_ylabel('Correlation between \nreal and generated samples', fontsize=14)
    ax2.set_ylabel('Rescaled RLD ('+kwd+')', fontsize=14)
    ax1.set_xlabel('Samples', fontsize=14)
    ax1.yaxis.set_label_coords(-0.055, 0.5)
    ax2.yaxis.set_label_coords(1.05, 0.5)
    ax1.xaxis.set_label_coords(.5, -0.25)
    plt.tight_layout()
    plt.xticks(fontsize=15)
    plt.show()

def genx_through_training(kwd, realx, genx, clist, cidx, spidx, glists):
    
    kwdidx = np.where(glists==kwd)[0][0]
    suptitle = title_ary[cidx]+'-'+str(spidx)+'th GE: '+kwd
    clist = clist[cidx]
    tpidx=-1 if cidx<3 else 0
    genx = genx[:,cidx,spidx-1,tpidx]
    
    plt.figure(figsize=(16,4))
    rx, gx = np.arange(len(clist)), np.arange(genx.shape[0])+10

    plt.bar(rx, realx[clist, kwdidx], width=.5)
    plt.bar(gx, genx[:, kwdidx], width=.5)
    plt.bar(rx[spidx-1], realx[clist[spidx-1], kwdidx], width=.5, color='turquoise')

    plt.legend(['Real','Generated'], ncol=2, bbox_to_anchor=(0., 0.18, 1.005, 1.0), fontsize=14)
    plt.ylim(0., 1.7)
    plt.title(suptitle, fontsize=14)
    plt.xlabel('Epoch', fontsize=14)
    plt.ylabel('Rescaled RLD', fontsize=14)
    plt.yticks(np.arange(0, 2., 0.4), fontsize=12)
    plt.grid(alpha=0.3)
    plt.xticks([2,10,110], (title_ary[cidx],'75K','125K'), fontsize=12)
    plt.plot()
    
def transition_curves(genx_A, genx_S, glists, flist):
    ncol = len(flist)
    ageidx = 2 #7M
    plt.figure(figsize=(6*ncol,5))
    for fl in range(len(flist)):
        plt.subplot(1,ncol,fl+1)
        plt.subplots_adjust(wspace=0.2)

        for gs in range(len(flist[fl])):
            idx_ = np.where(glists==flist[fl][gs])[0][0]
            plt.plot(genx_A[ageidx][:,idx_], label=flist[fl][gs], linewidth=3)
            y1, y2 = genx_A[ageidx][:,idx_]-genx_S[ageidx][:,idx_], genx_A[ageidx][:,idx_]+genx_S[ageidx][:,idx_]
            plt.fill_between(range(len(genx_A[ageidx][:,idx_])), y1, y2, alpha=0.2)
        plt.xticks([0,100], ['WT','AD'], fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylabel('Rescaled RLD', fontsize=14)
        plt.legend(fontsize=12, loc='upper right')
    plt.plot()
    
def transition_curves_all(realx, genx_A, genx_S, rld_mean, std_re_rld, max_rld_std, clist, glists, flist):
    age_ary = np.array([['2M', '4M', '7M'],['blue', 'orange', 'green']])
    X_WT2M, X_WT4M, X_WT7M = np.full(6,-2), np.full(6,0), np.full(6,2)
    X_AD2M, X_AD4M, X_AD7M = X_WT2M+100, X_WT4M+100, X_WT7M+100
    tickary = np.array([[X_WT2M, X_AD2M], [X_WT4M, X_AD4M], [X_WT7M, X_AD7M]])

    for fl in range(len(flist)):
        plt.figure(figsize=(24,5))
        for j in range(len(flist[fl])):
            plt.subplot(1,4,j+1)
            plt.subplots_adjust(wspace=0.27)

            idx_ = np.where(glists==flist[fl][j])[0][0]
            imin, imax = np.min(realx[:,idx_]), np.max(realx[:,idx_])

            AA, BB = 2*1.959*std_re_rld*max_rld_std[idx_], rld_mean[idx_]
            for ag in range(age_ary.shape[1]):
                plt.scatter(tickary[ag][0], realx[clist[3+ag],idx_], c=age_ary[1][ag])
                plt.scatter(tickary[ag][1], realx[clist[ag],idx_], c=age_ary[1][ag])
            for i in range(len(genx_A)):
                x = BB+AA*(genx_A[i][:,idx_]-.5) #transform to rawscale
                plt.plot(x)
                y1, y2 = x-genx_S[i][:,idx_], x+genx_S[i][:,idx_]
                plt.fill_between(range(len(genx_A[i][:,idx_])), y1, y2, alpha=.2)
            plt.title(flist[fl][j], fontsize=14)
            plt.ylabel('RLD', fontsize=14)
            plt.xticks([0,100], ['WT','AD'], fontsize=14)
            plt.yticks(fontsize=14)
        if fl==0:
            plt.legend(['2M', '4M', '7M'], bbox_to_anchor=(0.24, 0.07, 0.8, 1.15), ncol=3, fontsize=12)

def heatmap(corr):
    p = corr.shape[1]
    corr_ = np.zeros((p,p))
    for g1 in range(p):
        for g2 in range(p):
            corr_[g1,g2]=scipy.stats.pearsonr(corr[:,g1],corr[:,g2])[0]
    gheat_ = sns.clustermap(corr_, metric='euclidean', method='complete', xticklabels="", yticklabels="")
    gheat_.ax_row_dendrogram.set_visible(False)
    gheat_.ax_col_dendrogram.set_visible(False)
    
def transition_curves_with_plist(plist, texts, tres, genx_A, glists):
    for nl in range(len(plist)):
        plt.figure(figsize=(6,5))
        if plist[nl][:2]=='WP':
            idx_=0
        elif plist[nl][:2]=='mm':
            idx_=1
        elif plist[nl][:2]=='GO':
            idx_=2
        
        rgs, rsd, ri = analysis.get_egidx(texts[idx_], tres[idx_][0], tres[idx_][1], tres[idx_][2], np.where(tres[idx_][1]==plist[nl])[0][0], glists)
        #plot
        for g in range(len(ri)):
            plt.plot(genx_A[2][:,ri[g]], label=glists[ri[g]])
            plt.title(rgs+' : '+rsd)
            plt.ylabel('Rescaled RLD')
            plt.xticks([0,100], ['WT','AD'])
            
def heatmap_with_plist(plist, texts, tres, genx_A, glists, avg_mode=0):
    avg_exp = []
    for nl in range(len(plist)):
        #plt.figure(figsize=(6,5))
        if plist[nl][:2]=='WP':
            idx_=0
        elif plist[nl][:2]=='mm':
            idx_=1
        elif plist[nl][:2]=='GO':
            idx_=2

        rgs, rsd, ri = analysis.get_egidx(texts[idx_], tres[idx_][0], tres[idx_][1], tres[idx_][2], np.where(tres[idx_][1]==plist[nl])[0][0], glists)
        gexp, glabel = genx_A[2][:,ri].T, glists[ri]
        #plot
        fig = sns.clustermap(gexp, method='complete', cbar_kws=dict(label='Z-score', ticks=np.arange(-2,3,1)), vmin=-2, vmax=2, figsize=(10, len(ri)/1.5), z_score=0, col_cluster=False, row_cluster=True, xticklabels="", yticklabels=glabel, cmap='RdBu_r')
        fig.ax_row_dendrogram.set_visible(False)
        fig.cax.set_position([1.03,.2,.03,.45])
        fig.ax_heatmap.set_title(rgs+' : '+rsd)
        fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(), rotation=0)
        
        if avg_mode==1:
            tmp = np.array(fig.data2d)
            avg_exp.append(np.average(tmp,axis=0))
            plt.close()
    if avg_mode==1:
        avgexp = np.array(avg_exp)
        fig = sns.clustermap(avgexp, method='complete', cbar_kws=dict(label='Z-score', ticks=np.arange(-2., 3., 1)), vmin=-2, vmax=2, col_cluster=False, row_cluster=True, figsize=(10, len(avgexp)/1.5), xticklabels="", yticklabels=plist, cmap='RdBu_r')
        fig.ax_row_dendrogram.set_visible(False)
        fig.cax.set_position([1.07, .2, .03, .45])

def FDR_ER(ud, ttitle, mat_fdr, mat_er):
    if ud==0:
        col='Reds'
    elif ud==1:
        col='Blues'
    plt.figure()
    fig = sns.clustermap(mat_fdr, mask=(mat_fdr==0), linecolor='black', linewidths=0.5, cbar_kws=dict(label='-log10(FDR)', ticks=np.arange(0,10,2)), figsize=(2, len(ttitle)/1.5), vmin=0., vmax=8, row_cluster=False, col_cluster=False, xticklabels=['Up','P1','P2','P3'], yticklabels=ttitle, cmap=col)
    fig.ax_row_dendrogram.set_visible(False)
    fig.cax.set_position([5, .2, .1, .1])
    
    plt.figure()
    fig = sns.clustermap(mat_er, mask=(mat_er==0), linecolor='black', linewidths=0.5, cbar_kws=dict(label='log2(ER)', ticks=np.arange(0,8,2)), figsize=(2, len(ttitle)/1.5), vmin=0., vmax=4, row_cluster=False, col_cluster=False, xticklabels=['Up','P1','P2','P3'], yticklabels=ttitle, cmap=col)
    fig.ax_row_dendrogram.set_visible(False)
    fig.cax.set_position([5, .2, .1, .1])

def fold_changes(dpath, ext_glist, basis):
    ad2m, ad4m, ad7m  = np.arange(6), np.arange(6,12), np.arange(12,18)
    wt2m, wt4m, wt7m = np.arange(18,24), np.arange(24,30), np.arange(30,36)
    aidx = np.array([ad2m, ad4m, ad7m, wt2m, wt4m, wt7m])

    dat_ddsn = pd.read_csv(dpath)
    h_ddsn = dat_ddsn.columns
    dat_ddsn = np.array(dat_ddsn)
    eid_ddsn, gid_ddsn, cont_ddsn = dat_ddsn[:,0], dat_ddsn[:,1], dat_ddsn[:,2:].astype(np.float32)
    #print (dat_ddsn.shape, h_ddsn.shape, eid_ddsn.shape, gid_ddsn.shape, cont_ddsn.shape)
    gname = gid_ddsn

    #extract cont of ext_glist
    ext_gidx = []
    for i in range(len(ext_glist)):
        idx_ = np.where(ext_glist[i]==gname)[0][0]
        ext_gidx.append(idx_)
    #print (ext_gidx)

    cont = []
    for i in range(len(aidx)):
        cont.append(cont_ddsn[:, aidx[i]])
    cont = np.array(cont)
    #print (cont.shape)
    ext_cont = cont[:, ext_gidx, :]
    cont = ext_cont

    #get basis norm
    cont_avg = np.average(cont, axis=2)[basis]
    
    #norm
    ncont = cont.copy()
    for i in range(ncont.shape[0]):
        ncont[i]/=cont_avg.reshape((len(cont_avg),1))

    #get avg, std
    ncont_avg = np.average(ncont, axis=2)
    ncont_std = np.std(ncont, axis=2)
    
    #plot_gene levels
    xt = np.arange(1,ncont_avg.shape[1]*3+1,3)
    plt.figure(figsize=(20,9))
    plt.grid(True)
    
    plt.bar(xt-1, ncont_avg[4], width=0.25, color='black', yerr=ncont_std[4], capsize=3) #WT4M
    plt.bar(xt-0.5, ncont_avg[1], width=0.25, color='red', yerr=ncont_std[1], capsize=3) #AD4M
    plt.bar(xt, ncont_avg[5], width=0.25, color='gray', yerr=ncont_std[5], capsize=3) #WT7M
    plt.bar(xt+0.5, ncont_avg[2], width=0.25, color='blue', yerr=ncont_std[2], capsize=3) #AD7M

    plt.xlim(-1, ncont_avg.shape[1]*3)
    plt.ylim(0., 2.)
    plt.xticks(xt-0.25, gname, fontsize=16)
    plt.yticks([0., 0.5, 1., 1.5, 2.], fontsize=16)
    plt.legend(['WT4M', 'AD4M', 'WT7M', 'AD7M'], fontsize=16)
    plt.ylabel('Gene levels\n(fold change relative to WT4M)', fontsize=20)
    plt.title('GSE104775_about_WP103_15G', fontsize=20)
    #plt.savefig('test.png', bbox_inches='tight')
    plt.show()
    
    #plot_fold changes
    xt = np.arange(1,ncont.shape[0]*3+1,3)
    for gidx in range(ncont.shape[1]):
        plt.figure()
        plt.scatter((np.arange(0,6)*0.05)+xt[0], ncont[3][gidx], marker='o', color='royalblue', edgecolor='black', s=70)
        plt.errorbar(xt[0], np.average(ncont[3][gidx]), np.std(ncont[3][gidx]), color='black', capsize=5)

        plt.scatter((np.arange(0,6)*0.05)+xt[1], ncont[0][gidx], marker='o', color='orangered', edgecolor='black', s=70)
        plt.errorbar(xt[1], np.average(ncont[0][gidx]), np.std(ncont[0][gidx]), color='black', capsize=5)

        plt.scatter((np.arange(0,6)*0.05)+xt[2], ncont[4][gidx], marker='o', color='royalblue', edgecolor='black', s=70)
        plt.errorbar(xt[2], np.average(ncont[4][gidx]), np.std(ncont[4][gidx]), color='black', capsize=5)

        plt.scatter((np.arange(0,6)*0.05)+xt[3], ncont[1][gidx], marker='o', color='orangered', edgecolor='black', s=70)
        plt.errorbar(xt[3], np.average(ncont[1][gidx]), np.std(ncont[1][gidx]), color='black', capsize=5)

        plt.scatter((np.arange(0,6)*0.05)+xt[4], ncont[5][gidx], marker='o', color='royalblue', edgecolor='black', s=70)
        plt.errorbar(xt[4], np.average(ncont[5][gidx]), np.std(ncont[5][gidx]), color='black', capsize=5)

        plt.scatter((np.arange(0,6)*0.05)+xt[5], ncont[2][gidx], marker='o', color='orangered', edgecolor='black', s=70)
        plt.errorbar(xt[5], np.average(ncont[2][gidx]), np.std(ncont[2][gidx]), color='black', capsize=5)

        plt.xlim(0., 17)
        if gidx==1:
            plt.ylim(0., 2.5)
            plt.yticks([0.,0.5,1.,1.5,2.,2.5])
        else:
            plt.ylim(0.,2.)
            plt.yticks([0.,0.5,1.,1.5,2.])
        plt.xticks(xt, (['WT2M', 'AD2M', 'WT4M', 'AD4M', 'WT7M', 'AD7M']), fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylabel('Fold change', fontsize=12)
        plt.title(ext_glist[gidx], fontsize=14)
        #plt.savefig(sub_path+'/Fig_of_WP103_7G/bigver/GSE104775_WP103_16G_'+ext_glist[gidx]+'.png', bbox_inches='tight')
        plt.show() 

def pattern_of_curves(dat, clists, path_, tv, cond):
    dat_gname, dat_value = dat[0], dat[1:].astype(np.float32)
    
    #auto module
    fclist, fcidx, orderidx = [], [], []
    for i in range(len(clists)):
        flist = os.listdir(path_)
        clist, cidx = [], []
        for item in flist:
            if (item.find(clists[i]) is not -1):
                if item[-3:] != 'csv':
                    orderidx.append(item[-5])
                    tmp = np.loadtxt(path_+item, delimiter='\t', dtype='str')
                    clist.append(tmp)
                    tcidx = []
                    if (len(tmp) != 1) or (tmp != 'none'):
                        for j in range(len(tmp)):
                            idx_ = np.where(tmp[j]==dat_gname)[0][0]
                            tcidx.append(idx_)
                        cidx.append(tcidx)
                    else:
                        cidx.append('none')
        fclist.append(clist)
        fcidx.append(cidx)
    
    #ordering
    orderidx = np.array(orderidx).astype(np.int)-1
    ffclist, ffcidx = [], []
    for i in range(len(orderidx)):
        idx = np.where(orderidx==i)[0][0]
        ffclist.append(fclist[0][idx])
        ffcidx.append(fcidx[0][idx])
    fclist, fcidx = [ffclist], [ffcidx]

    #plot
    dlist = np.array([dat_value])
    for cl in range(len(clists)):
        tfcl, tfci = fclist[cl], fcidx[cl]
        plt.figure(figsize=(13,7))
        sname = path_+'new_6tv_ver3_'+clists[cl]+'_module.png'

        for i in range(2):
            for j in range(3):
                pidx = (3*i)+(j+1)
                plt.subplot(2,3,pidx)
                plt.subplots_adjust(wspace=0.2, hspace=0.3)
                if tfci[pidx-1]=='none': #tv
                    plt.title('Pattern'+str(pidx)+'_no member')
                    plt.ylim(-0.1, 1.1)
                    plt.plot(tv[:, pidx-1], c='red', linewidth=3)
                    plt.xticks([0,100], ('WT','AD'))
                else: #genx
                    plt.title('Pattern'+str(pidx)+' ('+str(len(tfcl[pidx-1]))+')', fontsize=16)
                    for m in range(len(tfcl[pidx-1])):
                        plt.plot(dlist[cl][:, tfci[pidx-1][m]], c='gray', alpha=0.3)
                        plt.plot(tv[:, pidx-1], c='red', linewidth=3)
                        plt.ylim(-0.1, 1.1)
                        plt.yticks(fontsize=17)
                    if ((pidx==1) or (pidx==4)) :
                        plt.ylabel('Rescaled RLD', fontsize=18)
                    plt.xticks([0,100], ('WT','AD'), fontsize=18)
