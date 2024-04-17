import os
import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import anndata as ad 
import re 
random.seed(10)
from scipy import sparse
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import pairwise_distances
from collections import defaultdict
import anndata as ad
from copy import copy
import pandas as pd
from path import Path
from sagenet.utils import glasso
import numpy as np


# data_path = 'data_tidy/00_anndata_objects_ult'

import torch
if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)

# atlas_visium = sc.read_10x_mtx('data_raw/liver_atlas/rawData_mouseStStVisium/countTable_mouseStStVisium')


VISIUM_DATA_DIR    = Path('data_tidy/composed_anndata_objects_2023_12_01')

ANNOTS_DIR        = Path('data_raw/2023-11-24_AllAnnosFinal')
NODULE_ANNOT_DIR  = ANNOTS_DIR / 'Nodules'
mapping_reps   = pd.read_csv(( NODULE_ANNOT_DIR / 'mapping_replicates.csv'), sep=';')


markers_all = {
      "Fibroblasts"  : ["Col1a1", "Col3a1", "Col1a2", "Col5a1", "Dpt", "Gas6", "Thbs1", "Col6a3", 'Sparc', 'Mmp2', 'Tagln', 'Thbs1', 'Postn'],
      "Erythrocytes" : ["Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1", "Snca", "Hemgn", "Gata1", "Rhd"],
      "Kupffer cells": ["Clec4f", "Csf1r", "Marco", "Cd163", "C1qc", "C1qb", "C1qa", "Ctss", "Ctsc"],
      "Neutrophils"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
      "Mast cells"   : ["Cpa3", "Cma1", "Tpsb2", "Mcpt4", 'Tpsab1'],
    #   "NK cells"     : ["Itgax", "S100a4", "Gzmb", "Gzma", "Lat", "Cd27", "Rgs1"],
      "Platelets"    : ["Pf4", "Itga2b", "Ppbp"],
      "B cells"      : ["Igkc", "Jchain"],
    "Periportal HCC": ["Sds", "Sdsl", "Hal", 'Bex1', 'Bex2', 
                    'Upp2', 'Aspg', 'Serpina6', 'Hamp2', 
                    'Ptger3', 'Bdh2'],
    "Pericentral HCC": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2", "Gstm3", "Gstm2", "Axin2"],
    "Cholangiocytes": ["Krt19", "Cldn7", "Krt7", "Mapk13", 'Epcam', "Gp2", 
                       "Slc15a2", 'Kcnk1', "Gp2", "Ezr", 
                       "Sox9", 'Spp1', 'Ppp1r1b', 'Car9', 'Tspan8'],
    "Other": ['Serpina6', 'Vil1', 'Car3', 'Upp2'],
    "Histone": ["Hist1h3d","Hist1h2ag", "Hist1h1c", "Hist1h1e", "Hist1h4h", "Hist4h4", "Hist1h3c"]
}


matrix_dict = defaultdict(dict)
for cell_type in markers_all.keys():
    for marker in markers_all[cell_type]:
        matrix_dict[marker][cell_type] = 0
        
        
        
for cell_type in markers_all.keys():
    for marker in markers_all[cell_type]:
        matrix_dict[marker][cell_type] = 1
        
        
matrix_df = pd.DataFrame.from_dict(matrix_dict, orient="index")
matrix_df = matrix_df.reindex(sorted(matrix_df.index), axis=0).fillna(0)


    
adatas     = {}
genes_all  = {}
for sample_f in os.listdir(VISIUM_DATA_DIR):
    if sample_f != 'ML_II_B_2Cyt.h5ad':
        mapping_df   = mapping_reps.loc[mapping_reps['sample_name'] == sample,['old_annotation', 'new_annotation']]
        mapping_dict = mapping_df.set_index('old_annotation').to_dict()['new_annotation']
        print(sample)
        adatas[sample] = sc.read(VISIUM_DATA_DIR / (sample + '.h5ad'))
        adata          = adatas[sample]
        adata.obs['histo_annotation_num'] = [mapping_dict.get(item, 'NaN') for item in adata.obs['histo_annotation']]
        sample = sample_f.strip('\\.h5ad')
        adatas[sample] = sc.read_h5ad(VISIUM_DATA_DIR / sample_f)
        adata          = adatas[sample]
        adata.obs['sample'] = sample
        sc.pp.filter_genes(adata, max_counts=1e06)
        # sc.pp.normalize_total(adata, target_sum=1e4)
        adata_copy = adata.copy()
        sc.pp.normalize_total(adata_copy, target_sum=1e4)
        sc.pp.log1p(adata_copy)
        sc.pp.highly_variable_genes(adata_copy, n_top_genes=1000, flavor='seurat_v3')
        hvgs = adata.var_names[adata_copy.var['highly_variable']]   
        genes          = list((((set(matrix_df.index) ) | set(hvgs)) )& set(adata.var_names))
        adata          = adata[:, genes]
        # sc.pp.filter_genes(adata, max_counts=1e6)
        # sc.pp.filter_genes(adata, min_counts=1e2)
        adatas[sample] = adata
        genes_all[sample]  = adata.var_names

  
adata_joint = ad.concat(adatas)
genes = adata_joint.var_names

for sample in adatas.keys():
    adata = adatas[sample][:, genes]
    le = LabelEncoder()
    adata.X = adata.X.toarray()
    adata.obs['nodule_encoded'] = le.fit_transform(adata.obs['histo_annotation_num'])
    glasso(adata)
    if not 'Cyt' in sample:
        adata.X = adata.X.toarray()
        sg_obj.add_ref(adata, comm_columns=['nodule_encoded'], tag=sample, epochs=100, verbose = True)
    
if not 'Cyt' in sample:
	sg_obj.add_ref(adata, comm_columns=['nodule_encoded'], tag=sample, epochs=100, verbose = True)


os.makedirs('objects/sagenet_model_all_sections_nodules_25_03_2024')
sg_obj.save_model_as_folder('objects/sagenet_model_all_sections_nodules_25_03_2024')
# adata_q = ad.concat([adata_dic[i] for i in ['ML_II_A_2', 'ML_I']])
adata_q = ad.concat(adata_dic)
sc.pp.combat(adata_q, key='section')

# sg_obj.load_model_as_folder('objects/sagenet_model_all_sections_nodules_17_08_2023')

sg_obj.map_query(adata_q, save_prob=True)



adata_q.write('data_tidy/01_sagenet_integration/integrated_query_nodules_23_11_2023.h5ad')


adata_q = sc.read('data_tidy/01_sagenet_integration/integrated_query_nodules_23_11_2023.h5ad')


import numpy as np

adata_q.obsm['prob_combined'] = np.column_stack((
		np.nan_to_num(adata_q.obsm['prob_ML_III_B_nodule_encoded'],0), 
		np.nan_to_num(adata_q.obsm['prob_ML_I_nodule_encoded'],0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_A_2_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_A_1_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_III_A_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_B_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_I_2_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_C_nodule_encoded'], 0)
		# np.nan_to_num(adata_q.obsm['prob_ML_III_A_2Cyt_nodule_encoded'], 0),
		# np.nan_to_num(adata_q.obsm['prob_ML_II_B_2Cyt_nodule_encoded'], 0),
		# np.nan_to_num(adata_q.obsm['prob_ML_II_A_3Cyt_nodule_encoded'], 0),
  		# np.nan_to_num(adata_q.obsm['prob_ML_II_B_3Cyt_nodule_encoded'], 0)

	)) 



base_dists = np.zeros((adata_q.obsm['prob_combined'].shape[0], adata_q.obsm['prob_combined'].shape[0]))

del adata_q.obsm['prob_combined']
prob_list = ['prob_ML_III_A_nodule_encoded', 'prob_ML_III_B_nodule_encoded', 'prob_ML_II_A_1_nodule_encoded', 'prob_ML_II_A_2_nodule_encoded', 'prob_ML_II_B_nodule_encoded', 'prob_ML_II_C_nodule_encoded', 'prob_ML_I_2_nodule_encoded', 'prob_ML_I_nodule_encoded']
for prob in prob_list:
	print(prob)
	pd = pairwise_distances(adata_q.obsm[prob])
	del adata_q.obsm[prob]
	pd /= np.linalg.norm(pd, 'fro')
	base_dists += pd
 
adata_q.obsp['sagenet_dist'] = base_dists
adata_q.write('data_tidy/01_sagenet_integration/integrated_query_nodules_umapped_23_11_2023.h5ad')

	


# import umap
# my_model = umap.UMAP(metric='precomputed')
# gg= my_model.fit(adata_q.obsp['sagenet_dist'])


sc.tl.leiden(adata_q, key_added = "leiden_raw")
# sc.tl.umap(dist_adata)

k=50
adata_q.obsp["distances"], adata_q.obsp["connectivities"] = sc.neighbors._compute_connectivities_umap(
    *sc.neighbors._get_indices_distances_from_dense_matrix(adata_q.obsp['sagenet_dist'], k),
    adata_q.obsp['sagenet_dist'].shape[0],
    k,
)
adata_q.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": None}}
del adata_q.obsp['sagenet_dist']

sc.tl.leiden(adata_q, 0.1, key_added = "leiden_sagenet_low")
sc.tl.leiden(adata_q, 0.5, key_added = "leiden_sagenet_medium")
sc.tl.leiden(adata_q, 1, key_added = "leiden_sagenet_high")
sc.tl.umap(adata_q)
sc.pl.umap(adata_q, color='leiden_sagenet_high', save='leiden_high.pdf')
sc.pl.umap(adata_q, color='section', save='_section.pdf')


# adata_q.obsm['X_umap_sagenet'] = gg.embedding_

del adata_q.uns['neighbors']
sc.pp.neighbors(adata_q)
# adata_q_subset = copy(adata_q[:,plasmid_genes])
# sc.pp.neighbors(adata_q, use_rep='X')
sc.tl.umap(adata_q)
sc.pl.umap(adata_q, color='section', save='_integrated_nodules.pdf')
sc.pl.umap(adata_q, color=plasmid_genes, save='_integrated_nodules_genes.pdf')



adata_q.write('data_tidy/01_sagenet_integration/integrated_query_nodules_umapped.h5ad')



del adata_q.obsp


k=50
adata_q.obsp["distances"], adata_q.obsp["connectivities"] = sc.neighbors._compute_connectivities_umap(
    *sc.neighbors._get_indices_distances_from_dense_matrix(base_dists, k),
    base_dists.shape[0],
    k,
)
adata_q.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": None}}
sc.tl.umap(adata_q)
sc.pl.umap(adata_q, color='section', save='_integrated_all_additive_normalized.pdf')
sc.pl.umap(adata_q, color=plasmid_genes, save='_integrated_all_additive_normalized_genes.pdf')


from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform


pairwise_distances(adata_q.obsm['prob_ML_III_B_nodule_encoded'])





sc.tl.tsne(adata_q, use_rep='prob_combined')
sc.pl.tsne(adata_q, color='section', save='_integrated_all_uncorrected_normalized.pdf')
sc.pl.tsne(adata_q, color=plasmid_genes, save='_integrated_all_uncorrected_normalized_genes.pdf')


from copy import copy
adata_q_subset = copy(adata_q[:,plasmid_genes])
sc.pp.neighbors(adata_q_subset, use_rep='X')
sc.tl.umap(adata_q_subset)
sc.pl.umap(adata_q_subset, color='section', save='_raw_all_subset.pdf')
sc.pl.umap(adata_q_subset, color=plasmid_genes, save='_raw_all_subset_genes.pdf')






vis_clus(
    spe = samples[[1]],
    clustervar = "nodule",
    colors = .palette_all, 
    size=0.5,
    alpha=0.5
) + theme(legend.position='none')



