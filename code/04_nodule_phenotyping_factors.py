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
# GENO_DIR        = Path('data_tidy/geno_data')

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


# markers_all = {
#     "hepatocyte": ["Alb", "Ttr", "Apoa1", "Serpina1c", "Fabp1", "Echs1", "Glul", "Acly", "Asl", "Cyp2e1", "Cyp2f2", "Ass1", "Mup3", "Pck1", "G6pc", "Apoa1", "Ass1", "G6pc", "Mup3"],
#     "hybrid_hepatocyte": ["Meg3", "Igfbp3", "Trp53inp1"],
#     "midzonal_hepatocyte": ["Oat", "Cyp2e1", "Cyp1a2", "Alb"],
#     "periportal_hepatocyte": ["Sds", "Sdsl", "Hal", "Apoa1", "Apoa2", "Apoa5", "Apob", "Alb", "Cyp2f2"],
#     "hepatic_stellate_cell": ["Ecm1", "Colec11", "Vipr1", "Hgf", "Rgs5", "Lrat", "Ngfr", "Reln", "Pth1r", "Col1a1", "Acta2", "Col3a1", "Dcn"],
#     "hepatocellular_carcinoma": ["Trf", "Serpina1a", "Orm1", "Hnf4a", "Fbp1", "Mat1a", "Sult1a1", "Apoa4", "Apob", "Cyp1a2", "Cyp3a11", "Cyp3a25", "Ugt1a1", "Fmo5", "Scd1", "Afp", "Gpc3"],
#     "stellate_activated": ["Ecm1", "Ccl2", "Col3a1", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"],
#     "stellate_activated_MYCi": ["Ecm1", "Trp53inp1", "Sfpq", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"],
#     "stellate_fibrotic": ["Ecm1", "Lrat", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"],
#     "stellate_quiescent": ["Ecm1", "Tagln", "Lrat", "Pdgfra", "Pdgfrb", "Rgs5", "Col14a1", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"],
#     "Endothelial": [
#         "Ptprb", "Vwf", "Klf4", "Sdc1", "Eng", "Ehd3", "Cd300lg", "Ramp2", "Fam167b", "Pecam1", "Oit3", "Esam",
#         "Sema6a", "Kdr", "Adam23", "Fcgr2b", "Gpr182", "Lrrc8a", "Tek", "Flt1", "Flt4", "Nrp1", "Nrp2"
#     ],
#     "B Cell": [
#         "Cd79a", "Cd79b", "Cd74", "Cd19", "Fcmr", "Jchain", "Mzb1", "Igkc", "Cd22", "Ebf1", "Pax5", "Cxcl10",
#         "Vpreb3", "Vpreb1"
#     ],
#     "Erythroblast": [
#         "Hbb-bs", "Hba-a2", "Hba-a1", "Klf1", "Alad", "Blvrb", "Hmbs", "Rhd", "Gata1"
#     ],
#     "Kupffer Cell": [
#         "Adgre1", "Emr1", "Clec4f", "Cd68", "Irf7", "C1qa", "C1qb", "C1qc", "Csf1r"
#     ],
#     "Macrophage": [
#         "Ccl9", "Cd14", "Marco", "Adgre1", "Cd68", "Csf1r", "Ptprc", "Cd52", "H2-Aa"
#     ],
#     "Monocytes": [
#         "Itgam", "Ly6g", "Ly6c1", "Ccr2"
#     ],
#     "Natural Killer": [
#         "Zap70", "Il2rb", "Nkg7", "Cxcr6", "Itga1"
#     ],
#     "T Cell": [
#         "Ccl5", "Cd7", "Cd3g", "Trbc1", "Trbc2", "Thy1", "Lat", "Cd3d", "Cdse", "Nkg7", "Ptprc", "Cd68", "Cd52",
#         "Cxcr6", "Cd163l1", "Itgae"
#     ],
#     "Plasmacytoid DCs": [
#         "Bcl11a", "Runx2", "Ccr9", "Siglech", "Spib", "Irf8"
#     ],
#     "Conventional DCs": [
#         "Irf5", "Irf8", "Ccl22", "Il12b", "Ccr7", "Id2", "Id3", "Siglech", "Ly6d", "Bst2", "Itgax", "Cd80", "Cd83"
#     ],
#     "Mesenchymal Cells": [
#         "Ecm1", "Colec11", "Vim", "Col1a2", "Mest", "Mmp2", "Pdgfra", "Ncam1", "Lhx2", "Fn1", "Cdh2", "Sparc",
#         "Tagln", "Tpm2", "Acta2", "Cnn1", "Pln", "Myh11"
#     ],
#     "Megakaryocytes": [
#         "Pf4", "Ppbp", "Cd9", "Itga2b", "Plek", "Cxcr4", "Itgb3"
#     ],
#     "Neutrophils": ['Sepx1', 'Retnlg', 'C1qa', 'Xcr1', 'Lcn2', 'Lgals3', 'Csf1', 'Ly6g', 'Csf3r', 'Slpi', 'Timd4',
#         'Cd209a', 'Cebpe', 'Cxcl2', 'Ccl4', 'Elane', 'Siglech', 'Ngp', 'Jun', 'Junb', 'S100a4', 'F13a1', 'Mpo', 
#         'S100a8', 'Cst7', 'Ccl3', 'Zfp36', 'Cd5l', 'Vsig4', 'Clec4f', 'Apoe', 'Gda', 'C1qc', 'Itgam', 'Chil3', 
#         'S100a9', 'Folr2', 'Fos', 'Adgre1', 'C1qb', 'Fcnb', 'Ltf'],
#     "Cholangiocytes": [
#         "Epcam", "Krt19", "Spp1", "Hnf1b", "Prom1", "St14", "Foxj1", "Cftr", "Sctr", "Car2", "Cd44", "Sox9", "Krt7"
#     ],
#     "Fibroblast-like": [
#         "Gsn", "Clec3b", "Dpt", "Cd34", "Mfap4", "Entpd2", "Fbln2", "Col15a1", "Ccnb2", "Apoa2", "C3", "Cpa3",
#         "Hdc", "Srgn", "Tpsb2", "Poln", "Pga5", "Fosb", "Cma1", "Tagln", "Tpm2", "Acta2", "Cnn1", "Pln", "Myh11"
#     ],
#     "Mast Cells": [
#         "Cpa3", "Hdc", "Srgn", "Tpsb2", "Poln", "Pga5", "Fosb", "Cma1"
#     ],
#     "Epithelial-like": [
#         "Krt8", "Krt19", "Mmp7", "Krt18", "Muc1", "Krt23", "Epcam", "Cldn6", "Cldn7", "Cdh1", "Utf1", "Epacm", 
#         "Pou5f1", "Dnmt3b"
#     ]
# }

# markers_all = {
#     "hepatocyte": ["Alb", "Ttr", "Apoa1", "Serpina1c", "Fabp1"],
#     "periportal_hepatocyte": ["Oat", "Sds", "Sdsl", "Hal", "Apoa2", "Apoa5","Cyp2f2"],
#     "hepatocellular_carcinoma": ["Trf", "Orm1", "Hnf4a", "Fbp1", "Mat1a"],
#     "Endothelial": [
#         "Ptprb", "Vwf", "Klf4", "Sdc1", "Eng", "Nrp1", "Nrp2"
#     ],
#     "B Cell": [
#         "Cd79a", "Cd79b", "Cd74", "Cd19", "Fcmr", "Jchain", "Cd22", "Ebf1", "Pax5", "Cxcl10"
#     ],
#     "Erythroblast": [
#         "Hbb-bs", "Hba-a2", "Hba-a1", "Klf1"
#     ],
#     "Kupffer Cell": [
#         "Adgre1", "Emr1", "Clec4f", "Cd68", "Irf7", "C1qa", "C1qb", "C1qc", "Csf1r"
#     ],
#     "Macrophage": [
#         "Ccl9", "Cd14", "Marco", "H2-Aa"
#     ],
#     "Monocytes": [
#         "Itgam", "Ly6g", "Ly6c1", "Ccr2"
#     ],
#     "Natural Killer": [
#         "Zap70", "Il2rb", "Nkg7", "Cxcr6", "Itga1"
#     ],
#     "T Cell": [
#         "Ccl5", "Cd7", "Cd3g", "Trbc1", "Trbc2", "Itgae", "Cdse", "Nkg7", "Ptprc", "Cd68", "Cd52"
#     ],
#     "Conventional DCs": [
#         "Irf5", "Irf8", "Ccl22", "Il12b", "Ccr7", "Id2", "Id3", "Siglech", "Ly6d", "Bst2", "Itgax", "Cd80", "Cd83"
#     ],
#     "Mesenchymal Cells": [
#         "Ecm1", "Colec11", "Vim", "Col1a2", "Mest", "Mmp2", "Pdgfra", "Ncam1", "Lhx2", "Fn1", "Cdh2", "Sparc",
#         "Tagln", "Tpm2", "Acta2", "Cnn1", "Pln", "Myh11"
#     ],
#     "Megakaryocytes": [
#         "Pf4", "Ppbp", "Cd9", "Itga2b", "Plek", "Cxcr4", "Itgb3"
#     ],
#     "Neutrophils": ['Elane', 'Siglech', 'Ngp', 'Mpo', 'S100a8', 'Ly6G', 'Lft', 'mmp8', 'S100a9'],
#     "Cholangiocytes": [
#         "Epcam", "Krt19", "Spp1", "Cd44", "Sox9", "Krt7"
#     ],
#     "Fibroblast-like": [
#         "Gsn", "Clec3b", "Dpt", "Cd34", "C3", 'Col1a1', 'lgfbp5', 'Col3a1', "Entpd2", "Fbln2", "Col15a1", "Ccnb2"
#     ],
#     "Mast Cells": [
#         "Cpa3", "Hdc", "Srgn", "Tpsb2", "Poln", "Pga5", "Fosb", "Cma1"
#     ],
#     "Epithelial-like": [
#         "Krt8", "Krt18", "Muc1", "Krt23", "Cldn6", "Cldn7", "Cdh1", "Utf1", "Pou5f1"
#     ]
# }

matrix_dict = defaultdict(dict)
for cell_type in markers_all.keys():
    for marker in markers_all[cell_type]:
        matrix_dict[marker][cell_type] = 0
        
        
        
for cell_type in markers_all.keys():
    for marker in markers_all[cell_type]:
        matrix_dict[marker][cell_type] = 1
        
        
matrix_df = pd.DataFrame.from_dict(matrix_dict, orient="index")
matrix_df = matrix_df.reindex(sorted(matrix_df.index), axis=0).fillna(0)

# ORGANS = ['Connective tissue', 'Immune system', 'Liver', 'Blood']
# MARKERS_PDB_F    = Path('data_raw/PanglaoDB_markers.tsv')
# markers_df2 = pd.read_csv(MARKERS_PDB_F, sep='\t')
# markers_df2 = markers_df2[(markers_df2['species'].isin(['Mm Hs', 'Mm'])) & \
#                 (markers_df2['sensitivity_mouse'] > .8)  &\
#                 (markers_df2['specificity_mouse'] < .2) &
#                 (markers_df2['organ'].isin(ORGANS))
#             ]
# markers_df2['gene_id'] = markers_df2['official gene symbol'].apply(lambda x: x.capitalize())
# markers_df2 = markers_df2[['gene_id', 'cell type']]
# markers_df2['cell type'] = markers_df2['cell type'].str.replace(' ', '_')
# markers_df2 = markers_df2.rename({'cell type': 'cell_type'})
# matrix_df2 = pd.get_dummies(markers_df2.set_index('gene_id'), columns=['cell type'], prefix='cell type') * 1
# matrix_df2.columns = matrix_df2.columns.str.replace('cell type_', '')
# matrix_df2 = matrix_df2[~matrix_df2.index.duplicated(keep='first')]


    
adatas     = {}
genes_all  = {}
for sample_f in os.listdir(VISIUM_DATA_DIR):
    if sample_f != 'ML_II_B_2Cyt.h5ad':
        sample = sample_f.strip('\\.h5ad')
        adatas[sample] = sc.read_h5ad(VISIUM_DATA_DIR / sample_f)
        adata          = adatas[sample]
        adata.obs['sample'] = sample
        sc.pp.filter_genes(adata, max_counts=1e06)
        # sc.pp.normalize_total(adata, target_sum=1e4)
        adata_copy = adata.copy()
        sc.pp.normalize_total(adata_copy, target_sum=1e4)
        sc.pp.log1p(adata_copy)
        sc.pp.highly_variable_genes(adata_copy, n_top_genes=500, flavor='seurat_v3')
        hvgs = adata.var_names[adata_copy.var['highly_variable']]   
        genes          = list((((set(matrix_df.index) ) | set(hvgs)) )& set(adata.var_names))
        # genes          = list((((set(matrix_df.index) )))& set(adata.var_names))
        # genes          = list((set(matrix_df.index) & set(adata.var_names)))
        # genes = hvgs
        adata          = adata[:, genes]
        sc.pp.filter_genes(adata, max_counts=1e6)
        sc.pp.filter_genes(adata, min_counts=1e2)
        # genes          = list(set(matrix_df.index) & set(adata.var_names))
        # adata.varm['celltype_mask'] = matrix_df.loc[genes,:]
        # adata.var['cell_type'] = adata.varm['celltype_mask'].loc[genes,:].idxmax(1)
        adatas[sample] = adata
        genes_all[sample]  = adata.var_names
    

adata_joint = ad.concat(adatas)

tech_category_key=None
# plot_category_keys=['nodule']

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
# waypoints = ['Alb', 'C1qb', 'Col1a1', 'Cd74', 'Hbb-bs', 'Adgre1', 'H2-Aa', 'Marco', 'Mpo', 'Vwf','Nkg7', 'Cd7', 'Acta2', 'S100a8', 'Krt19', 'Cma1', 'Fosb', 'Krt23', 'Sds', 'Dpt']
# adata_joint_g.obs['is_waypoint_size'] = 1
# adata_joint_g.obs['is_waypoint'] = False
# adata_joint_g.obs.loc[waypoints, 'is_waypoint'] = True
# adata_joint_g.obs.loc[waypoints, 'is_waypoint_size'] = 10
# n_factors = 20

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
    # 10X reaction / sample / batch
    batch_key='sample',
    labels_key=None,
    # multiplicative technical effects (platform, 3' vs 5', donor effect)
    # categorical_covariate_keys=['Donor'],
    # unexplained variance categories (normally batch)
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
    # plan_kwargs={
    #    "n_aggressive_epochs": 1000,
    #    "n_aggressive_steps": 10,
    # },
    #ignore_warnings=True,
    use_gpu=True,
    # use_aggressive_training=True
)


# Save model
mod.save(f"{ref_run_name}", overwrite=True)
# Cell2ModuleModel.load(f"{ref_run_name}", adata_joint)

# Save anndata object with results
adata_file = f"{ref_run_name}/sp.h5ad"
adata_joint.write(adata_file)
# adata_joint = sc.read(adata_file)

factor_names_keys = {
    "g_fg": "factor_names",
    "cell_modules_w_cf": "factor_names",
}
adata_test= mod.export_posterior(
    adata_test, use_quantiles=True,
    # choose quantiles
    # add_to_varm=["means", "stds"],
    # add_to_varm=["q05","q50", "q95"],
    add_to_varm=['means', "q50"],
    # choose variables to export
    export_varm_variables=["g_fg"],
    export_obsm_variables=["cell_modules_w_cf"],
    sample_kwargs={'batch_size': 2000, 'use_gpu': True, 'exclude_vars': ['data_rna', 'data_chromatin']},
    factor_names_keys=factor_names_keys,
)

adata_file = f"{ref_run_name}/sp_test.h5ad"
adata_test.write(adata_file)

adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf_q_high'] = \
    1 / np.quantile(adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf'], 0.9995, axis=0).flatten()

def factor_expressed_plot(x, y, x_cut=4, y_cut=15,
                          x_lab='1 / 99.9% quantile of cell loadings',
                          y_lab='1 / 99.9% quantile of gene loadings',
                          invert_selection=False):
  # Expression shape and rate across cells
  plt.scatter(x, y);
  plt.xlabel(x_lab)
  plt.ylabel(y_lab)
  #plt.vlines(x_cut, 0, y_cut)
  #plt.hlines(y_cut, 0, x_cut)
  low_lab = high_lab = 'active'
  if not invert_selection:
    high_lab = f'not {high_lab}'
  else:
    low_lab = f'not {low_lab}'
  #plt.text(x_cut - 0.5 * x_cut, y_cut - 0.5 * y_cut, low_lab)
  #plt.text(x_cut + 0.1 * x_cut, y_cut + 0.1 * y_cut, high_lab)
  return {'x': x, 'y': y, 'x_cut': x_cut, 'y_cut': y_cut}


cutoff = 1.0
r = factor_expressed_plot(x_cut = cutoff, y_cut = cutoff,
                               x=adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf_q_high'],
                               y=adata_test.uns['mod']['post_sample_q50']['cell_modules_w_cf_q_high'],
                              )

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
# type_factors = np.matmul(matrix_df.loc[adata_test.varm['q50_g_fg'].index.tolist(),:].T, adata_test.varm['q50_g_fg'])

# type_factors = np.matmul(matrix_df.loc[adata_test.varm['q50_g_fg'].index.tolist(),:].T, df)


type_factors = np.matmul(matrix_df.loc[genes,:].T, df.loc[genes,:])


a1 = adata_test.uns['mod']['selected_factor_names']


a2 = type_factors.idxmax(0).tolist()

sorted_values = np.sort(type_factors.T, axis=1)

# Keep only the top 5 values in each row
top_values = sorted_values[:, -1]

# Create a mask to identify elements that are part of the top 5
mask = np.isin(type_factors.T, top_values)
# a2 = ['Mix', 'Cholangiocytes', 'Endothelial', 'Mix', 
#       'HCC', 'Neutrophils/Macrophages', 
#       'Mix', 'Megakaryocytes', 
#       'Mix', 'Lymphocytes', 'Periportal (Sds, Sdsl)', 
#       'Endothelial', 'Erythroblast', 'Fibro/Mes', 'Mix', 
#       'Periportal (Cyp2f2, Pck1, Hal)', 'Neutrophils (S100a8, S100a9)', 
#       'Fibro/Mes (Col*)', 'Mix']

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
  