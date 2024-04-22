import cell2module
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
# import cell2location
import scvi
# import txnburst
import hashlib
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns
import cell2module
from pathlib import Path
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder  
from cell2module.tl._initialisation import (
    find_waypoint_gene_clusters, compute_w_initial_waypoint,
    align_plot_stability, compute_pcs_knn_umap, knn_building, _subset_cells
)
from collections import defaultdict
import anndata as ad
from copy import copy
from cell2module.models import Cell2ModuleModel
from numpy.linalg import inv

results_folder = './results/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/ref_2023_12_08'
run_name = f'{results_folder}/run_map_2023_12_08'

VISIUM_DATA_DIR    = Path('data_tidy/composed_anndata_objects_2023_12_01')

markers_all = {
    "Fibroblasts"  : ["Col1a1", "Col3a1", "Col1a2", "Col5a1", "Dpt", "Gas6", "Thbs1", "Col6a3", 'Sparc', 'Mmp2', 'Tagln', 'Thbs1', 'Postn'],
    "Erythrocytes" : ["Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1", "Snca", "Hemgn", "Gata1", "Rhd"],
    "Kupffer cells": ["Clec4f", "Csf1r", "Marco", "Cd163", "C1qc", "C1qb", "C1qa", "Ctss", "Ctsc"],
    "Neutrophils"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
    "Mast cells"   : ["Cpa3", "Cma1", "Tpsb2", "Mcpt4", 'Tpsab1'],
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
        sample = sample_f.strip('\\.h5ad')
        adatas[sample] = sc.read_h5ad(VISIUM_DATA_DIR / sample_f)
        adata          = adatas[sample]
        adata.obs['sample'] = sample
        sc.pp.filter_genes(adata, max_counts=1e06)
        adata_copy = adata.copy()
        sc.pp.normalize_total(adata_copy, target_sum=1e4)
        sc.pp.log1p(adata_copy)
        sc.pp.highly_variable_genes(adata_copy, n_top_genes=500, flavor='seurat_v3')
        hvgs = adata.var_names[adata_copy.var['highly_variable']]   
        genes          = list((((set(matrix_df.index) ) | set(hvgs)) )& set(adata.var_names))
        adata          = adata[:, genes]
        sc.pp.filter_genes(adata, max_counts=1e6)
        sc.pp.filter_genes(adata, min_counts=1e2)
        adatas[sample] = adata
        genes_all[sample]  = adata.var_names
    

adata_joint = ad.concat(adatas)

tech_category_key=None

adata_joint.uns['mod'] = dict()
adata_joint.layers['counts'] = adata_joint.X

adata_joint_pc = compute_pcs_knn_umap(
    adata_joint,
    plot_category_keys='sample',
    scale_max_value=10, n_comps=20, n_neighbors=15,
)

adata_joint_g, n_factors = find_waypoint_gene_clusters(
    adata_joint_pc,
    k='aver_norm',
    n_factors=40,
    margin_of_error=1,
    n_neighbors=20,
    labels_key=None,
    label_filter=None,
    verbose=False,
)


adata_joint = compute_w_initial_waypoint(
    adata_joint, adata_joint_g, n_factors,
    scale=True,
    # tech_category_key=tech_category_key,
    use_x=False, layer='counts',
    knn_smoothing=True,
)


adata_joint.uns['mod_init'] = adata_joint.uns['mod_init'].copy()


init_vals = {
    'cell_modules_w_cf': adata_joint.uns['mod_init']['initial_values']['w_init']['cell_factors_w_cf'],
}
init_vals = {k: v.values.astype('float32') for k, v in init_vals.items()}
n_factors = init_vals[list(init_vals.keys())[0]].shape[1]


adata_test = copy(adata_joint)
adata_test.obs_names = adata_test.obs_names + "_" + adata_test.obs['sample'].astype(str)

init_vals = {
    'cell_modules_w_cf': adata_test.uns['mod_init']['initial_values']['w_init']['cell_factors_w_cf'],
}
init_vals = {k: v.values.astype('float32') for k, v in init_vals.items()}
n_factors = init_vals[list(init_vals.keys())[0]].shape[1]

cell2module.models.Cell2ModuleModel.setup_anndata(
    adata=adata_joint,
    layer='counts',
    batch_key='sample',
    labels_key=None,
    variance_categories="sample",
)

import torch
model_kwargs = {
    "rna_model": True,
    "chromatin_model": False,
    "use_non_linear_decoder": False,
    "bayesian_decoder": False,
    "n_hidden": 1024,
    "n_layers": 2,
    "use_orthogonality_constraint": False,

    "detection_hyp_prior": {"alpha": 20.0, "mean_alpha": 1.0, "mean_beta": 1.0},

    "amortised": False,
    "encoder_mode": "multiple",
    "encoder_kwargs": {'dropout_rate': 0.1,
                    'n_hidden': {
                "single": 1024,
                "detection_y_c": 10,
                "detection_chr_y_c": 10,
                "factors_per_cell": 3,
                "cell_modules_w_cf": 1024,
                "cell_modules_w_cf_amortised_prior": 1024,
                "cell_type_modules_w_cz_prior": 1024,
            },
                    'use_batch_norm': False, 'use_layer_norm': True,
                    'n_layers': 2, 'activation_fn': torch.nn.ELU,
    },
}

mod = Cell2ModuleModel(
    adata_joint,
    n_factors=n_factors,
    init_vals=init_vals,
    **model_kwargs,
)
# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(
    max_epochs=10000, batch_size=None,
    lr=0.002,
    use_aggressive_training=False,
    use_gpu=True,
)



factor_names_keys = {
    "g_fg": "factor_names",
    "cell_modules_w_cf": "factor_names",
}
adata_test= mod.export_posterior(
    adata_test, use_quantiles=True,
    add_to_varm=['means', "q50"],
    export_varm_variables=["g_fg"],
    export_obsm_variables=["cell_modules_w_cf"],
    sample_kwargs={'batch_size': 2000, 'use_gpu': True, 'exclude_vars': ['data_rna', 'data_chromatin']},
    factor_names_keys=factor_names_keys,
)

adata_file = f"{ref_run_name}/sp_test.h5ad"
adata_test.write(adata_file)

adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf_q_high'] = \
    1 / np.quantile(adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf'], 0.9995, axis=0).flatten()


r = {
        'x': adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf_q_high'], 
        'y': adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf_q_high'], 
        'x_cut': 1, 
        'y_cut': 1
     }

adata_test.uns['mod']['factor_selection'] = (r['x'] < r['x_cut']) & (r['y'] < r['y_cut'])
adata_test.uns['mod']['selected_factor_names'] = list(np.array(adata_test.uns['mod']['factor_names'])[adata_test.uns['mod']['factor_selection']])




cell_loadings = adata_test.obsm['q50_cell_modules_w_cf'][
    [
        f'q50_cell_modules_w_cf_{i}'
        for i in adata_test.uns['mod']['selected_factor_names']
    ]
]

df = adata_test.varm['q50_g_fg'].iloc[:,adata_test.uns['mod']['factor_selection']]
df = (df.T / df.T.sum()).T


genes = list(set(matrix_df.index) & set(adata_test.varm['q50_g_fg'].index))


type_factors = np.matmul(matrix_df.loc[genes,:].T, df.loc[genes,:])


a1 = adata_test.uns['mod']['selected_factor_names']


a2 = type_factors.idxmax(0).tolist()

sorted_values = np.sort(type_factors.T, axis=1)

# Keep only the top 5 values in each row
top_values = sorted_values[:, -1]

# Create a mask to identify elements that are part of the top 5
mask = np.isin(type_factors.T, top_values)
cell_loadings.columns = list(map('::'.join, zip(a1, a2)))


cell_loadings = cell_loadings.loc[:, (type_factors.loc[:,adata_test.uns['mod']['factor_selection']].sum(0) > 2).tolist()]

factor_columns = [col for col in adata_test.obs.columns if 'factor' not in col]

# Create a new DataFrame with only the selected columns
adata_test.obs = adata_test.obs[factor_columns]

adata_test.obs[cell_loadings.columns] = cell_loadings
adata_test.obsm['selected_cell_loadings'] = cell_loadings
adata_test.obsm['celltype_loadings'] = np.matmul(cell_loadings, type_factors.T)
adata_test.obsm['celltype_loadings'].columns = type_factors.index
adata_test.obsm['celltype_loadings_masked'] = np.matmul(cell_loadings, mask)
adata_test.obsm['celltype_loadings_masked'].columns = type_factors.index
adata_test.obs.loc[:, type_factors.index] = adata_test.obsm['celltype_loadings_masked']


adata_test.write('data_tidy/joint_factors_2024_03_21_masked.h5ad')

df.columns = list(map('::'.join, zip(a1, a2)))
df = df.loc[genes,:]
top_indices = {}
for column in df.columns:
    top_indices[column] = df[column].nlargest(20).index.tolist()

pd.DataFrame.from_dict(top_indices).to_csv('data_tidy/top_genes_per_factor_2024_03_21_new.csv')