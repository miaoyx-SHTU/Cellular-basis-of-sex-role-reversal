#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc 
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
from samap.utils import save_samap
from samap.utils import load_samap
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pickle
import holoviews as hv
import anndata


# In[8]:


human_hy_15k = sc.read_h5ad('human_hy_15k.h5ad')
mouse_placenta_uterus_HY = sc.read_h5ad('mouse_placenta_uterus_HY.h5ad')
seahorse_symbol = sc.read_h5ad('seahorse_symbol.h5ad')


# In[11]:


human_hy_15k.X


# In[7]:


fn1 = 'human_hy_15k.h5ad'
fn2 = 'mouse_placenta_uterus_HY.h5ad'
fn3 = 'seahorse_symbol.h5ad'

filenames = {'hs':fn1,'mm':fn2,'hc':fn3}
keys = {'hs':'ident','mm':'ident','hc':'ident'}


# In[ ]:


sm = SAMAP(
        filenames,
        f_maps = 'example_data/maps/',
        keys=keys,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
    )

sm.run(pairwise=True,ncpus=20)


# In[ ]:


save_samap(sm, 'hyaline.SAMap_alignment.h5ad')


# In[27]:


fn2 = 'mouse_placenta_uterus_HY.h5ad'
fn3 = 'seahorse_symbol.h5ad'

filenames = {'mm':fn2,'hc':fn3}
keys = {'mm':'ident','hc':'ident'}
sm = SAMAP(
        filenames,
        f_maps = './ref_cds/maps/',
        keys=keys,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
            )


# In[18]:


help(SAMAP)


# In[42]:


raw_seahorse_symbol = seahorse_symbol.raw.to_adata()


# In[64]:


raw_var = pd.DataFrame(index=raw_seahorse_symbol.var['_index'].to_list())
raw_var['features'] = raw_seahorse_symbol.var['_index'].to_list()


# In[68]:


raw_write =   anndata.AnnData(X=seahorse_symbol.raw.X,obs=seahorse_symbol.obs,var=raw_var)


# In[69]:


raw_write.write('seahorse_symbol_raw.h5ad')


# In[14]:


raw_human_hy_15k = human_hy_15k.raw.to_adata()

raw_write =   anndata.AnnData(X=human_hy_15k.raw.X,obs=human_hy_15k.obs,var=raw_human_hy_15k.var)

raw_write.write('human_hy_15k_raw.h5ad')


# ## seahorse $ mouse

# In[ ]:


fn1 = 'seahorse_symbol_raw.h5ad'
fn2 = 'mouse_placenta_uterus_HY.h5ad'

filenames = {'hc':fn1,'mm':fn2}
keys = {'hc':'ident','mm':'ident'}
sm = SAMAP(
        filenames,
        f_maps = './ref_cds/maps/',
        keys=keys,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
            )


# In[ ]:


sm.run(pairwise=True,ncpus=20)


# In[ ]:


save_samap(sm, 'seahorse_mouse_SAMap_alignment')


# In[56]:


#https://github.com/pandas-dev/pandas/issues/56825
with open('seahorse_mouse_SAMap_alignment.pkl','rb') as f:
    #sm = pickle.load(f)
    sm = pd.read_pickle(f) 


# In[57]:


keys = {'hc':'ident','mm':'ident'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
MappingTable.to_csv('seahorse_mouse_MappingTable_0.txt', sep='\t')


# In[58]:


#sankey plot
#https://github.com/flatironinstitute/CaImAn/issues/1140
a=sankey_plot(MappingTable, align_thr=0.2, species_order = ['mm','hc'])
hv.extension('matplotlib')
#hv.save(a, filename="seahorse_mouse_sankey_0.15.html")
hv.save(a, filename='seahorse_mouse_sankey_0.2.pdf', fmt='pdf')

#github.com/holoviz/holoviews/issues/1819
renderer = hv.renderer('bokeh')

# Using renderer save
renderer.save(a, 'seahorse_mouse_sankey_0.2.html')


# In[59]:


select_list = ['mm_M_EPC','mm_M_TGC','mm_M_pTGC','hc_S_TGC(cxcl14)','hc_S_TGC(muc5)','hc_S_TGC(nucb2)','hc_S_decidual']
select_MappingTable=MappingTable.loc[select_list,select_list]


# In[60]:


select_MappingTable


# In[61]:


#sankey plot
#https://github.com/flatironinstitute/CaImAn/issues/1140
a=sankey_plot(select_MappingTable, align_thr=0.2, species_order = ['mm','hc'])
hv.extension('matplotlib')
#hv.save(a, filename="seahorse_mouse_sankey_0.15.html")
hv.save(a, filename='seahorse_mouse_select_sankey_0.2.pdf', fmt='pdf')

#github.com/holoviz/holoviews/issues/1819
renderer = hv.renderer('bokeh')

# Using renderer save
renderer.save(a, 'seahorse_mouse_select_sankey_0.2.html')


# In[53]:


#Calculating table of enriched gene pairs for cell type mappings with an alignment score above 0.4:
gpf = GenePairFinder(sm,keys=keys)
gene_pairs = gpf.find_all(align_thr=0.1)
gene_pairs.to_csv('seahorse_mouse_SAMap_gene_pairs.txt',  sep='\t')


# ## seahorse & human

# In[16]:


fn1 = 'seahorse_symbol_raw.h5ad'
fn2 = 'human_hy_15k_raw.h5ad'

filenames = {'hc':fn1,'hs':fn2}
keys = {'hc':'ident','hs':'ident'}
sm = SAMAP(
        filenames,
        f_maps = './ref_cds/maps/',
        keys=keys,
        save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
            )


# In[17]:


sm.run(pairwise=True,ncpus=20)


# In[18]:


save_samap(sm, 'seahorse_human_SAMap_alignment')


# In[62]:


#https://github.com/pandas-dev/pandas/issues/56825
with open('seahorse_human_SAMap_alignment.pkl','rb') as f:
    #sm = pickle.load(f)
    sm = pd.read_pickle(f) 


# In[63]:


sm


# In[64]:


keys = {'hc':'ident','hs':'ident'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
MappingTable.to_csv('seahorse_human_MappingTable_0.txt', sep='\t')


# In[65]:


#sankey plot
#https://github.com/flatironinstitute/CaImAn/issues/1140
a=sankey_plot(MappingTable, align_thr=0.2, species_order = ['hs','hc'])
hv.extension('matplotlib')
#hv.save(a, filename="seahorse_mouse_sankey_0.15.html")
hv.save(a, filename='seahorse_human_sankey_0.2.pdf', fmt='pdf')

#github.com/holoviz/holoviews/issues/1819
renderer = hv.renderer('bokeh')

# Using renderer save
renderer.save(a, 'seahorse_human_sankey_0.2.html')


# In[66]:


select_list = ['hs_H_EVT','hs_H_SCT','hs_H_VCT','hc_S_TGC(cxcl14)','hc_S_TGC(muc5)','hc_S_TGC(nucb2)','hc_S_decidual']
select_MappingTable=MappingTable.loc[select_list,select_list]


# In[67]:


select_MappingTable


# In[68]:


a=sankey_plot(select_MappingTable, align_thr=0.2, species_order = ['hs','hc'])
hv.extension('matplotlib')
hv.save(a, filename="seahorse_human_select_sankey_0.2.html")
hv.save(a, filename='seahorse_human_select_sankey_0.2.pdf', fmt='pdf')


# In[69]:


#github.com/holoviz/holoviews/issues/1819
renderer = hv.renderer('bokeh')

# Using renderer save
renderer.save(a, 'seahorse_human_select_sankey_0.2.html')


# In[50]:


#Calculating table of enriched gene pairs for cell type mappings with an alignment score above 0.4:
gpf = GenePairFinder(sm,keys=keys)
gene_pairs = gpf.find_all(align_thr=0.1)
gene_pairs.to_csv('seahorse_human_SAMap_gene_pairs.txt',  sep='\t')


# In[107]:


seahorse_human_SAMap_gene_pairs = pd.read_csv('seahorse_human_SAMap_gene_pairs.txt',  sep='\t',index_col=0)
seahorse_mouse_SAMap_gene_pairs = pd.read_csv('seahorse_mouse_SAMap_gene_pairs.txt',  sep='\t',index_col=0)


# In[81]:


seahorse_human_SAMap_gene_pairs


# In[108]:


seahorse_mouse_SAMap_gene_pairs


# ## epi

# In[84]:


human_choose_type = 'hc_S_TGC(cxcl14);hs_H_EVT'
mouse_choose_type = 'hc_S_TGC(cxcl14);mm_M_FB'


# In[121]:


seahorse_human_SAMap_gene_pairs.loc[:,['hc_S_TGC(cxcl14);hs_H_EVT']]


# In[122]:


#human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,[i.split(';')[0] ==  choose_type for i in seahorse_human_SAMap_gene_pairs.columns]]
human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,['hc_S_TGC(cxcl14);hs_H_EVT']]
#human_choose_gene = human_choose_gene.iloc[:,[0,3]]
human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene]


# In[123]:


pd.DataFrame(human_choose_gene).to_csv('hc_S_TGC(cxcl14)_human.csv')


# In[124]:


#mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,[i.split(';')[0] ==  choose_type for i in seahorse_mouse_SAMap_gene_pairs.columns]]
mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,['hc_S_TGC(cxcl14);mm_M_FB']]
#mouse_choose_gene = mouse_choose_gene.iloc[:,[0]]
#mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene]


# In[126]:


pd.DataFrame(mouse_choose_gene).to_csv('hc_S_TGC(cxcl14)_mouse.csv')


# In[ ]:





# In[188]:


human_choose_type = 'hc_S_TGC(muc5);hs_H_EVT'
mouse_choose_type = 'hc_S_TGC(muc5);mm_M_EC'


# In[189]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
#human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_TGC(muc5)_human.csv')


# In[190]:


mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,mouse_choose_type]
#mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene.dropna()]
pd.DataFrame(mouse_choose_gene).to_csv('hc_S_TGC(muc5)_mouse.csv')


# In[ ]:





# In[191]:


human_choose_type = 'hc_S_TGC(nucb2);hs_H_VCT'
#mouse_choose_type = 'hc_S_TGC(muc5);mm_M_EC'


# In[192]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
#human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_TGC(nucb2)_human.csv')


# In[ ]:





# In[193]:


human_choose_type =[ 'hc_S_decidual;hs_H_SCT','hc_S_decidual;hs_H_VCT']
mouse_choose_type = ['hc_S_decidual;mm_M_EPC','hc_S_decidual;mm_M_pTGC']


# In[194]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_decidual_human.csv')


# In[195]:


mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,mouse_choose_type]
mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene.dropna()]
pd.DataFrame(mouse_choose_gene).to_csv('hc_S_decidual_mouse.csv')


# ## immu

# In[196]:


human_choose_type = 'hc_S_immune;hs_H_immune'
mouse_choose_type = 'hc_S_immune;mm_M_immune'


# In[197]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
#human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_immune_human.csv')

mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,mouse_choose_type]
#mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene.dropna()]
pd.DataFrame(mouse_choose_gene).to_csv('hc_S_immune_mouse.csv')


# In[198]:


human_choose_type = 'hc_S_FB;hs_H_FB'
mouse_choose_type = 'hc_S_FB;mm_M_FB'


# In[199]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
#human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_FB_human.csv')

mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,mouse_choose_type]
#mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene.dropna()]
pd.DataFrame(mouse_choose_gene).to_csv('hc_S_FB_mouse.csv')


# In[200]:


human_choose_type = 'hc_S_endo;hs_H_endo'
mouse_choose_type = 'hc_S_endo;mm_M_endo'


# In[201]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
#human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_endo_human.csv')

mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,mouse_choose_type]
#mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene.dropna()]
pd.DataFrame(mouse_choose_gene).to_csv('hc_S_endo_mouse.csv')


# In[202]:


human_choose_type =[ 'hc_S_basal;hs_H_VCT','hc_S_progenitor;hs_H_EVT']
mouse_choose_type = ['hc_S_basal;mm_M_ES','hc_S_progenitor;mm_M_FB']


# In[203]:


human_choose_gene = seahorse_human_SAMap_gene_pairs.loc[:,human_choose_type]
human_choose_gene = human_choose_gene.melt()['value'].dropna()
human_choose_gene = [i.split(';')[1][3:] for i in human_choose_gene.dropna()]
pd.DataFrame(human_choose_gene).to_csv('hc_S_basal&progenitor_human.csv')

mouse_choose_gene = seahorse_mouse_SAMap_gene_pairs.loc[:,mouse_choose_type]
mouse_choose_gene = mouse_choose_gene.melt()['value'].dropna()
mouse_choose_gene = [i.split(';')[1][3:] for i in mouse_choose_gene.dropna()]
pd.DataFrame(mouse_choose_gene).to_csv('hc_S_basal&progenitor_mouse.csv')


# ## read

# In[2]:


#https://github.com/pandas-dev/pandas/issues/56825
with open('seahorse_mouse_SAMap_alignment.pkl','rb') as f:
    #sm = pickle.load(f)
    sm = pd.read_pickle(f) 


# In[5]:


smp = sm.samap


# In[6]:


smp.adata.varp["homology_graph_reweighted"]


# In[8]:


adata=smp.adata.copy()


# In[13]:


sc.pl.umap(adata,color=['species'])


# In[ ]:




