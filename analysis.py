import numpy as np
import scipy.stats
from model import WGAN_GP
import matplotlib.pyplot as plt
import pandas as pd
import os

"""
wgangp = WGAN_GP(xtr, n_te)
t_vars = tf.trainable_variables()
d_vars = [var for var in t_vars if var.name.startswith('discriminator')]
g_vars = [var for var in t_vars if var.name.startswith('generator')]
"""
##need to set path for plot(webgestalt results)
wgpath = 'rld1_gse104775_genx_247m_75k_125k/6TV_webgestalt/thr0.95_base0.9_60_10_10_50/'
up_flist = np.array(['WP_7M_up', 'KEGG_7M_up', 'GOBP_7M_up'])
dn_flist = np.array(['WP_7M_down', 'KEGG_7M_down', 'GOBP_7M_down'])

up_clist = np.array([[wgpath+'WP/7M_WP_type1',wgpath+'WP/7M_WP_type2',wgpath+'WP/7M_WP_type3','rld1_gse104775_genx_247m_75K_125K/Webgestalt/WikiPathway/7M_WP_up'],                  
                  [wgpath+'KEGG/7M_KEGG_type1', wgpath+'KEGG/7M_KEGG_type2', wgpath+'KEGG/7M_KEGG_type3', 'rld1_gse104775_genx_247m_75K_125K/Webgestalt/KEGG/7M_KEGG_up'],
                  [wgpath+'GOBP/7M_GOBP_type1', wgpath+'GOBP/7M_GOBP_type2', wgpath+'GOBP/7M_GOBP_type3', 'rld1_gse104775_genx_247m_75K_125K/Webgestalt/GO_BP/7M_GO_BP_up']])
dn_clist = np.array([[wgpath+'WP/7M_WP_type4',wgpath+'WP/7M_WP_type5',wgpath+'WP/7M_WP_type6','rld1_gse104775_genx_247m_75K_125K/Webgestalt/WikiPathway/7M_WP_down'],
                  [wgpath+'KEGG/7M_KEGG_type4', wgpath+'KEGG/7M_KEGG_type5', wgpath+'KEGG/7M_KEGG_type6', 'rld1_gse104775_genx_247m_75K_125K/Webgestalt/KEGG/7M_KEGG_down'],
                  [wgpath+'GOBP/7M_GOBP_type4', wgpath+'GOBP/7M_GOBP_type5', wgpath+'GOBP/7M_GOBP_type6', 'rld1_gse104775_genx_247m_75K_125K/Webgestalt/GO_BP/7M_GO_BP_down']])

def calculate_corr(dat1, dat2):
    #len(dat.shape)==3 : generated data
    #len(dat.shape)==2 : real data
    if len(dat1.shape)==3:
        if len(dat2.shape)==3: #gen vs. gen
            corr = np.zeros((dat1.shape[0], dat1.shape[1], dat2.shape[1]))
            for t in range(dat1.shape[0]):
                for n1 in range(dat1.shape[1]):
                    for n2 in range(dat2.shape[1]):
                        tmp1, tmp2 = dat1[t,n1], dat2[2,n2]
                        corr[t,n1,n2] = scipy.stats.pearsonr(tmp1,tmp2)[0]
        elif len(dat2.shape)==2: # gen vs. real
            corr = np.zeros((dat1.shape[0], dat1.shape[1], dat2.shape[0]))
            for t in range(dat1.shape[0]):
                for n1 in range(dat1.shape[1]):
                    for n2 in range(dat2.shape[0]):
                        tmp1, tmp2 = dat1[t,n1], dat2[n2]
                        corr[t,n1,n2] = scipy.stats.pearsonr(tmp1,tmp2)[0]
    elif len(dat1.shape)==2: # real vs. real
        corr = np.zeros((dat1.shape[0], dat2.shape[0]))
        for n1 in range(dat1.shape[0]):
            for n2 in range(dat2.shape[0]):
                tmp1, tmp2 = dat1[n1], dat2[n2]
                corr[n1,n2] = scipy.stats.pearsonr(tmp1,tmp2)[0]
    return corr

def generate_augx(re_rld, eps, wgangp, n_noise):
    wgangp.load_weight(eps)
    tz, tgenx = wgangp.generate_samples(n_noise)
    
    gen_z = []
    for nr in range(len(re_rld)):
        corr = np.zeros((n_noise,))
        for na in range(n_noise):
            corr[na] = scipy.stats.pearsonr(re_rld[nr], tgenx[na])[0]
        maxidx = np.argsort(corr)[-10:] #top10
        avgz = np.average(tz[maxidx], axis=0)
        #tmp_genx = wgangp.generate_samples(avgz.reshape((1,100)))
        #genaugx.append(tmp
        gen_z.append(avgz)
    genx = wgangp.generate_samples(gen_z)
    return genx

def latent_interpolation(path, wgangp, re_rld, g_vars, clist, ep_ary, n_noise=10000):
    """
    z:(epochs,types(WT/AD),age(2/4/7M),#augsamples,zdim)
    avgz:(epochs,ages(2/4/7M),zdim)
    genx:(epochs,ages(2/4/7M),#augsamples,interpolate_stpes,p)
    gencorr:(epochs,types(WT/AD),age(2/4/7M),#augsamples)
    fw:(epochs,fw.shape)
    """
    sname = path+'generated_from_75k_to_125k.npz'
    WTAD_z, WTAD_avgz, WTAD_genx, WTAD_gencorr, fw = [], [], [], [], []
    for ep in range(len(ep_ary)):
        #restore weights
        tpath = path+'ep'+str(ep_ary[ep])
        wgangp.load_weight(tpath)
        
        #generate_x
        tz, tgenx = wgangp.generate_samples(n_noise)
        gens, gene = [], []
        gen_scorr, gen_ecorr = [], []
        
        for cond in range(len(clist)):
            wt, ad = clist[cond][0], clist[cond][1]
            
            tmps, tmpe = [], []
            tmps_corr, tmpe_corr = [], []
            for sp in range(len(wt)):
                scorr_, ecorr_ = np.zeros((n_noise,)), np.zeros((n_noise,))
                for n in range(n_noise):
                    scorr_[n], ecorr_[n] = scipy.stats.pearsonr(re_rld[wt[sp]], tgenx[n])[0], scipy.stats.pearsonr(re_rld[ad[sp]], tgenx[n])[0]
                smaxidx, emaxidx = np.argsort(scorr_)[-10:], np.argsort(ecorr_)[-10:]
                tmps.append(np.average(tz[smaxidx], axis=0))
                tmpe.append(np.average(tz[emaxidx], axis=0))
                tmps_corr.append(np.average(scorr_[smaxidx]))
                tmpe_corr.append(np.average(ecorr_[emaxidx]))
            gens.append(tmps)
            gene.append(tmpe)
            gen_scorr.append(tmps_corr)
            gen_ecorr.append(tmpe_corr)
        gens, gene = np.array(gens), np.array(gene)
        gen_scorr, gen_ecorr = np.array(gen_scorr), np.array(gen_ecorr)
        delta_z = np.average(gene, axis=1)-np.average(gens, axis=1)
        inter = np.linspace(0., 1., num=101)
        WTAD_z.append([gens, gene])
        WTAD_avgz.append(delta_z)
        WTAD_gencorr.append([gen_scorr, gen_ecorr])
        
        #generate_x by intepolating latent vectors
        genx = []
        for ag in range(gens.shape[0]): #age
            tmpgenx = []
            for augsp in range(gens.shape[1]): #aug_samples
                ttmpgenx = []
                for z_ in inter: #interpolate
                    tmpz = gens[ag][augsp] + (delta_z[ag]*z_)
                    ttmpgenx.append(tmpz)
                ttmpgenx = wgangp.generate_samples(ttmpgenx)
                tmpgenx.append(ttmpgenx)
            genx.append(tmpgenx)
        genx = np.array(genx)
        WTAD_genx.append(genx)
        
        #generate_weights
        gen_fw = wgangp.get_weights(g_vars[4])
        fw.append(gen_fw)
        
    WTAD_z, WTAD_avgz, WTAD_genx, WTAD_gencorr, fw = np.array(WTAD_z), np.array(WTAD_avgz), np.array(WTAD_genx), np.array(WTAD_gencorr), np.array(fw)
    
    np.savez(sname, WTAD_z=WTAD_z, WTAD_avgz=WTAD_avgz, WTAD_genx=WTAD_genx, WTAD_gencorr=WTAD_gencorr, gen_fw=fw)
    print ('Done, ', sname)
    return WTAD_z, WTAD_avgz, WTAD_genx, WTAD_gencorr, fw

def extract_specific_cond_glists(corr_matrix, heatmap, threshold, ex_glists):
    #heatmap to ary
    gh_col = np.array(heatmap.data2d.columns)
    reorder = []
    for row in range(gh_col.shape[0]):
        tr = (corr_matrix[gh_col[row]])
        trc = (tr[gh_col])
        reorder.append(trc)
    reorder = np.array(reorder)
    reorder_glist = ex_glists[gh_col]
    #print (reorder.shape, reorder_glist.shape)
    plt.matshow(reorder)

    fmat = np.concatenate((reorder_glist.reshape((1, len(reorder_glist))), reorder), axis=0)
    np.fill_diagonal(fmat[1:], 0)
    #print (fmat.shape)
    
    #extract based on thr
    thr = threshold
    node1, node2, weight = [], [], []
    for i in range(len(fmat)-1):
        w = fmat[1:,i].astype(np.float)
        idx_ = np.where(w>=thr)[0]
        if len(idx_) != 0:
            node1.append(np.full(len(idx_), fmat[0,i]))
            node2.append(fmat[0][idx_])
            weight.append(w[idx_])

    node1, node2, weight = np.concatenate(node1), np.concatenate(node2), np.concatenate(weight)
    #print (node1.shape, node2.shape, weight.shape)

    df = pd.DataFrame({'node1':node1, 'node2':node2, 'weight':weight}, index=np.arange(len(node1)))
    df1 = pd.DataFrame(np.sort(df[['node1','node2']], axis=1))
    df = df[~df1.duplicated()]
    #print (df.shape)

    #save pd to csv
    df.to_csv('weights/75K_to_125K_thr'+str(thr)+'_weight_network.csv', sep=',', na_rep='NAN')
    print('Save, weights/weight_network.csv')

    #then, plot gene wetwork using cytoscape and this csv(input)
    
def parsing(cont):
    stextfile = cont.split('"')
    tgeneset, geneset, set_description = [], [], []
    size, ovlap, expect, ratio, pvalue, fdr = [], [], [], [], [], []
    for i in range(len(stextfile)):
        if (stextfile[i]=='description'):
            if (stextfile[i+1]==':'):
                tgeneset.append(stextfile[i-2]+'\":[{\"userId')
                geneset.append(stextfile[i-2])
                set_description.append(stextfile[i+2])
                size.append(stextfile[i+9][:-1].split(':')[1])
                ovlap.append(stextfile[i+11][:-1].split(':')[1])
                expect.append(stextfile[i+13][:-1].split(':')[1])
                ratio.append(stextfile[i+15][:-1].split(':')[1])
                pvalue.append(stextfile[i+17][:-1].split(':')[1])
                fdr.append(stextfile[i+19][:-1].split(':')[1])
    for i in range(len(stextfile)): #final
        if stextfile[i] == 'clusters':
            tgeneset.append(stextfile[i])
            geneset.append('eol')
            set_description.append('eol')
            size.append('eol')
            ovlap.append('eol')
            expect.append('eol')
            ratio.append('eol')
            pvalue.append('eol')
            fdr.append('eol')
    return np.vstack((tgeneset, geneset, set_description, size, ovlap, expect, ratio, pvalue, fdr))

def get_data_from_html(bpath, hlist):
    textfiles = []
    results = []
    for nl in range(len(hlist)):
        source=open(bpath+hlist[nl], 'r', encoding='utf-8')
        textfile = source.read()
        
        tres = parsing(textfile)
        textfiles.append(textfile)
        results.append(tres)
    return textfiles, results

def get_egidx(ttext, tdat, dat, datd, fidx_, glists):
    sidx, eidx = ttext.index(tdat[fidx_]), ttext.index(tdat[fidx_+1])
    etext = ttext[sidx:eidx]
    s_etext = etext.split('"')

    #get eid and match gname
    egidx = []
    for i in range(len(s_etext)):
        if s_etext[i]=='userId':
            #print (s_etext[i+2])
            idx_ = (np.where(s_etext[i+2]==glists)[0][0])
            egidx.append(idx_)
    egidx = np.array(egidx)
    #print (dat[fidx_], ',', datd[fidx_], ',' ,egidx.shape)
    return dat[fidx_], datd[fidx_], egidx

def get_extlist(lists):
    wp, kg, go = [], [], []
    for nl in range(len(lists)):
        if lists[nl][:2]=='WP':
            wp.append(lists[nl])
        elif lists[nl][:2]=='mm':
            kg.append(lists[nl])
        elif lists[nl][:2]=='GO':
            go.append(lists[nl])
    extlist = np.array([wp, kg, go])
    return extlist

def get_vcont(cond):
    vcont = []
    for i in range(len(cond)):
        tpath = cond[i]+'/'
        ilist = os.listdir(tpath)
        
        for item in ilist:
            if item.find('Report') is not -1:
                source = open(tpath+item, 'r', encoding='utf-8')
                textfile = source.read()

                #parsing
                tmp = parsing(textfile).T
                vcont.append(tmp)
    vcont = np.array(vcont)
    return vcont

def get_ninfo(target, list_, ucnt, dat):
    tresult = []
    
    idx_ = np.where(target==list_)[0][0]
    for p in range(len(ucnt[idx_])):
        ttresult = []
        tmp = ucnt[idx_][p]
        tmpcont = dat[tmp]
        tfidx = np.where(tmpcont[:,1]==target)[0]
        if len(tfidx)==0:
            print ('ERROR!',i)
        else:
            tfcont = dat[tmp][tfidx[0]] 
            ttresult.append(tmp) #idx
            ttresult.append(target) #geneset
            ttresult.append(tfcont[2]) #desc
            ttresult.append(tfcont[6]) #fdr
            ttresult.append(tfcont[8]) #er
        tresult.append(np.array(ttresult))
    return np.array(tresult)

def get_plotinfo(ud, lists, extlist):
    #ud-up:0, dn:1
    if ud==0:
        clists=up_clist
    elif ud==1:
        clists=dn_clist
        
    ext_results = []
    for case in range(len(clists)):
        textres = []
        vcont = get_vcont(clists[case])

        #concat cont except null
        fcont = []
        for i in range(len(vcont)):
            if len(vcont[i])!=0:
                fcont.append(vcont[i])
        fcont = np.concatenate(fcont, axis=0)

        #get union
        ulist = np.unique(fcont[:,1])

        #cnt about ulist
        ucnt = []
        for i in range(len(ulist)):
            tcnt = []
            for j in range(len(vcont)):
                if len(vcont[j])!=0:
                    idx_ = np.where(ulist[i]==vcont[j][:,1])[0]
                    if len(idx_)!=0:
                        tcnt.append(j)
            ucnt.append(tcnt)
        ucnt = np.array(ucnt)
        for i in range(len(extlist[case])):
            textres.append(get_ninfo(extlist[case][i], ulist, ucnt, vcont))
        ext_results.append(textres)
    ext_results = np.array(ext_results)
    
    #trimming
    ucase = []
    for i in range(len(ext_results)): #case
        tucase = []
        for j in range(len(ext_results[i])): #gs
            ufdr, uer = np.zeros((4,))-1, np.zeros((4,))-1
            for k in range(len(ext_results[i][j])): #pt
                #print (k, ext_results[i][j][k])
                ufdr[int(ext_results[i][j][k][0])]=ext_results[i][j][k][4]
                uer[int(ext_results[i][j][k][0])]=ext_results[i][j][k][3]
            tucase.append([ext_results[i][j][k][1], ext_results[i][j][k][2], ufdr, uer])
        ucase.append(tucase)
    #print (len(ucase))
    
    mat_fdr, mat_er, gs, desc = [], [], [], []
    for i in range(len(ucase)):
        tmp = ucase[i]
        for g in range(len(tmp)):
            gs.append(tmp[g][0])
            desc.append(tmp[g][1])
            mat_fdr.append(tmp[g][2])
            mat_er.append(tmp[g][3])
    gs, desc, mat_fdr, mat_er = np.array(gs), np.array(desc), np.array(mat_fdr), np.array(mat_er)

    #reseq : P1 P2 P3 UP -> UP P1 P2 P3
    reidx = [3,0,1,2]
    mat_fdr, mat_er = mat_fdr[:,reidx], mat_er[:,reidx]

    #alignment
    reidx = []
    for i in range(len(lists)):
        uidx_ = np.where(lists[i]==gs)[0][0]
        reidx.append(uidx_)

    gs, desc, mat_fdr, mat_er = gs[reidx], desc[reidx], mat_fdr[reidx], mat_er[reidx]
    #print (gs.shape, desc.shape, mat_fdr.shape, mat_er.shape)
    
    #transform and inf to zero
    mat_fdr, mat_er = -np.log10(mat_fdr), np.log2(mat_er)
    mat_fdr[mat_fdr==0]=0 #fdr=1
    mat_fdr[np.isnan(mat_fdr)]=0 #no
    mat_fdr[mat_fdr==np.inf]=20 #fdr=0
    mat_er[np.isnan(mat_er)]=0
    
    #title
    ttitle = []
    for i in range(len(gs)):
        t = desc[i]+'-'+gs[i]
        ttitle.append(t)
        
    return ttitle, mat_fdr, mat_er


