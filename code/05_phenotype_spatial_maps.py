import scanpy as sc
import squidpy as sq
import pandas as pd 
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder  
from scipy import stats
# from matplotlib.colors import LinearSegmentedColormap, shiftedColorMap

sc.set_figure_params(scanpy=True, fontsize=20)

clustered_nodules = pd.read_csv('output/clustred_nodules_q_.6_2024_03_22.csv')
factors_obs       = sc.read('data_tidy/joint_factors_2024_03_21_new.h5ad').obs

VISIUM_DATA_DIR = Path('data_tidy/composed_anndata_objects_2024_01_11')

ANNOTS_DIR        = Path('data_raw/2023-11-24_AllAnnosFinal')
NODULE_ANNOT_DIR  = ANNOTS_DIR / 'Nodules'
mapping_reps   = pd.read_csv(( NODULE_ANNOT_DIR / 'mapping_replicates.csv'), sep=';')


all_rois = ['ML_I_2', 'ML_II_B', 'ML_II_C', 'ML_III_A', 'ML_III_B', 'ML_II_A_1'] + ['ML_II_B_3Cyt', 'ML_II_A_3Cyt', 'ML_III_A_2Cyt', 'ML_I',
       'ML_II_A_2']

type_genes = {
    "Periportal HCC": ["Sds", "Sdsl", "Hal", 'Bex1', 'Bex2', 
                    'Upp2', 'Aspg', 'Serpina6', 'Hamp2', 
                    'Ptger3', 'Bdh2'],
    "Pericentral HCC": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2", "Gstm3", "Gstm2", "Axin2"],
    "Cholangiocytes": ["Krt19", "Cldn7", "Krt7", "Mapk13", 'Epcam', "Gp2", 
                       "Slc15a2", 'Kcnk1', "Gp2", "Ezr", 
                       "Sox9", 'Spp1', 'Ppp1r1b', 'Car9', 'Tspan8', 'Mup3', 'Serpina1c', 'Lgals2'],
    'Other': ['Hamp2', 'Hamp', 'Car3', 'Upp2'],
    "Histone": ["Hist1h3d","Hist1h2ag", "Hist1h1c", "Hist1h1e", "Hist1h4h", "Hist4h4", "Hist1h3c"]
}

    
type_bins = [s + '_bin' for s in type_genes.keys()]

flattened_tumor = [value for sublist in tumor_genes.values() for value in sublist]

def compute_quantiles(data):
    data = np.array(data, dtype=np.float32)
    sorted_data = np.sort(data)
    n = len(data)
    quantiles = np.zeros(n, dtype=np.float32)
    for i, value in enumerate(data):
        pos = np.where(sorted_data == value)[0][0]
        quantile = (pos + 1) / n
        quantiles[i] = quantile
    return quantiles

def normalize_array(arr):
    # Calculate the quantiles
    q25 = np.quantile(arr, 0.25, axis=0)
    q99 = np.quantile(arr, 0.99, axis=0)
    normalized_arr = (arr - q25) / (q99 - q25)
    normalized_arr = np.clip(normalized_arr, 0, 1)
    return normalized_arr

adatas = {}
adatas_geno = {}
for sample in all_rois:
    sample = sample.replace(".h5ad", "")
    mapping_df   = mapping_reps.loc[mapping_reps['sample_name'] == sample,['old_annotation', 'new_annotation']]
    mapping_dict = mapping_df.set_index('old_annotation').to_dict()['new_annotation']
    print(sample)
    adatas[sample] = sc.read(VISIUM_DATA_DIR / (sample + '.h5ad'))
    adata          = adatas[sample]
    adata.obs['number'] = [mapping_dict.get(item, 'NaN') for item in adata.obs['histo_annotation']]
    cl_nodules = clustered_nodules.loc[(clustered_nodules[['sample']]==sample).values,:]
    adata.obs['barcode'] = adata.obs_names
    new_obs = pd.merge(adata.obs, cl_nodules, on='number')
    new_obs.index = new_obs.loc[:,'barcode']
    for type in type_bins:
        adata.obs.loc[:, type] = 0
        adata.obs.loc[new_obs['barcode'].values, type] = new_obs.loc[:,type]
    adata.obsm['spatial'] = adata.obsm['spatial'].astype('float') 
    adata.uns['spatial'] = adata.uns['spatial']
    adata.obs['library_id'] = lib_id =  list(adata.uns['spatial'].keys())[0]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    factors_obs_sample = factors_obs.loc[(factors_obs[['sample']]==sample).values,:]
    factors_obs_sample.index = factors_obs_sample.index.str.split('_' + sample).str[0]
    for p in type_genes:  
        adata.obs[f'{p}_aggregate'] = normalize_array(compute_quantiles(np.array(adata[:,type_genes[p]].X.mean(1))[:,0]))
        adata.obs[f'{p}_factor']    = normalize_array(factors_obs_sample.loc[adata.obs_names, p]) # / factors_obs_sample.loc[adata.obs_names, p].max()
        adata.obs[f'{p}_binarised'] = adata.obs[f'{p}_bin']
        del adata.obs[f'{p}_bin']
    adatas[sample] = adata

type_dict = {}
cmaps = {} 
for type in type_genes:
    type_dict[type]  = []
    type_dict[type] += [f'{type}_aggregate', f'{type}_binarised', f'{type}_factor']
    type_dict[type] += type_genes[type]
    cmaps[f'{type}_aggregate'] = 'magma'
    cmaps[f'{type}_binarised'] = 'viridis'
    cmaps[f'{type}_factor']    = 'plasma'
    for g in type_genes[type]:
        cmaps[g] = 'coolwarm'

sc.set_figure_params(dpi=20)
for type in type_genes:
    plt.rcParams['pdf.fonttype'] = 42
    fig = plt.figure(layout='constrained', figsize=((5 * len(type_dict[type]), 5 * len(adatas))))
    subfigs = fig.subfigures(len(adatas), len(type_dict[type]), wspace=0.0)
    for si, s in enumerate(all_rois):
        for i, k in enumerate(type_dict[type]):
            ax = subfigs[si, i].subplots(1, 1, sharey=True)
            sc.pl.spatial(adatas[s], img_key='lowres', color=k,
                cmap=cmaps[k],vmin=0, vmax=1, ax=ax, show=False, colorbar_loc=None)
            # sc.pl.spatial(adata_dict[k], img_key="lowres", color='histo_decision', scale_factor=1, ax=ax, show=False)
            ax.set_title('')
            if si == 0:
                ax.set_title(k, fontsize=30)
            ax.set_rasterization_zorder(3)
            if si > 5:
                ax.invert_xaxis()
                ax.invert_yaxis()
            ax.axis('off')
    fig.savefig(f'output/phenotype_spatial_plots/{type}.pdf', dpi=50)