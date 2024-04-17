import scanpy as sc
import squidpy as sq
import pandas as pd 
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder  


RAW_VISIUM_DIR    = Path('data_raw/raw_visium_files')
ANNOTS_DIR        = Path('data_raw/2023-11-24_AllAnnosFinal')
NODULE_ANNOT_DIR  = ANNOTS_DIR / 'Nodules'
NODULE_ANNOT_DIR  = ANNOTS_DIR / 'Nodules'
PLASMID_ANNOT_DIR = ANNOTS_DIR / 'Plasmids'
IHC_ANNOT_DIR     = ANNOTS_DIR / 'IHC'
TYPE_ANNOT_DIR    = ANNOTS_DIR / 'Types'
SAVE_DIR          = Path('data_tidy/composed_anndata_objects_2024_01_30')
os.makedirs(SAVE_DIR)

mapping_reps   = pd.read_csv(( NODULE_ANNOT_DIR / 'mapping_replicates.csv'), sep=';')

plasmid_annots = pd.read_csv(PLASMID_ANNOT_DIR / 'Plasmid_DECODEv1.csv').loc[:,['List', 'ID']]
plasmid_annots.rename(columns={'ID': 'gene_ids', 'List': 'Plasmid'}, inplace=True)
olft_annots = pd.read_csv(PLASMID_ANNOT_DIR / 'OLFTASVMN_DECODEv1.csv').loc[:,['List', 'ID']]
olft_annots.rename(columns={'ID': 'gene_ids', 'List': 'OLFT'}, inplace=True)
olft_annots = olft_annots.iloc[:1398,:4]
bc_annots   = pd.read_csv(PLASMID_ANNOT_DIR / 'BC_DECODEv1.csv').loc[:,['List', 'ID']]
bc_annots.rename(columns={'ID': 'gene_ids', 'List': 'BC'}, inplace=True)


adatas = {}
nodule_annot = {}
for f in os.listdir(RAW_VISIUM_DIR):
    print(f)
    if 'ML' in f:
        # read file_names
        f_name         = f.rstrip(".h5ad")
        adatas[f_name] = sq.read.visium(RAW_VISIUM_DIR / f / 'outs')
        # nodule histopathologica annotation
        mapping_df   = mapping_reps.loc[mapping_reps['sample_name'] == f_name,['old_annotation', 'new_annotation']]
        mapping_dict = mapping_df.set_index('old_annotation').to_dict()['new_annotation']
        nodule_annot = pd.read_csv(NODULE_ANNOT_DIR / (f_name + '.csv'))
        assert np.all(nodule_annot['Barcode'] == adatas[f_name].obs_names), 'nodule barcodes are not equal to anndata barcodes'
        nodule_annot.fillna('normal_tissue', inplace=True)
        adatas[f_name].obs['histo_annotation'] = nodule_annot.iloc[:,1].values
        # adatas[f_name].obs['histo_annotation'] = adatas[f_name].obs['histo_annotation'].apply(lambda x: x.capitalize() if x != 'normal_tissue' else x)
        adatas[f_name].obs['histo_annotation'] =np.array([(lambda s: s[0].upper() + s[1:] if s != 'normal_tissue' else s)(s) for s in adatas[f_name].obs['histo_annotation'].values])
        adatas[f_name].obs['histo_annotation'] = adatas[f_name].obs['histo_annotation'].replace({'Nx1mb':'Nx1MB',  'Nx2mb':'Nx2MB', 'Nx3mb':'Nx3MB'})
        # adatas[f_name].obs['histo_annotation'] = adatas[f_name].obs['histo_annotation'].str.upper() 
        adatas[f_name].obs['histo_annotation_num'] = [mapping_dict.get(item, 'NaN') for item in adatas[f_name].obs['histo_annotation']]
        # IHC annotations
        for stain in ['GFP', 'GS', 'RFP']:
            f_ihc = IHC_ANNOT_DIR / (f_name + '_' + stain + '.csv')
            if os.path.isfile(f_ihc):
                stain_annot = pd.read_csv(f_ihc)
                assert np.all(stain_annot['Barcode'] == adatas[f_name].obs_names), 'staining barcodes are not equal to anndata barcodes'
                stain_annot = stain_annot.iloc[:,1].str.replace(r'%s' % stain, '')
                stain_annot = stain_annot.str.replace(r'%s' % stain.lower(), '')
                stain_annot = stain_annot.replace({'lo': 'low'})
                assert (stain_annot.isna() | stain_annot.isin(['low', 'moderate', 'hi'])).all()
                adatas[f_name].obs[stain + '_stain'] = stain_annot.values
        # types
        for type in ['CCC', 'SDS']:
            f_type = TYPE_ANNOT_DIR / (f_name + '_' + type + '.csv')
            if os.path.isfile(f_type):
                type_annot = pd.read_csv(f_type)
                assert np.all(type_annot['Barcode'] == adatas[f_name].obs_names), 'type barcodes are not equal to anndata barcodes'
                type_annot = type_annot.iloc[:,1]
                assert (type_annot.isna() | type_annot.isin([type.lower()])).all(), 'type'
                type_annot = type_annot.fillna(False)
                type_annot = type_annot.replace(type.lower(), True)
                adatas[f_name].obs[type + '_type'] = type_annot.values
        # gene annots
        adatas[f_name].var_names_make_unique()
        a = adatas[f_name].copy()
        a.var['symbol'] = a.var_names
        a.var_names = a.var.loc[:,'gene_ids']
        annots       = pd.merge(adatas[f_name].var, plasmid_annots , on='gene_ids', how='outer')
        annots       = pd.merge(olft_annots, annots, on='gene_ids', how='outer')
        annots       = pd.merge(annots, bc_annots, on='gene_ids', how='outer')
        annots       = annots[~annots['feature_types'].isna()]
        annots.index = a.var.loc[annots.gene_ids.values,'symbol']
        adatas[f_name].var= annots.loc[adatas[f_name].var_names,:]
        # store the object
        if f == 'ML_II_B_3Cyt':
            nan_coords = np.argwhere(np.isnan(adatas[f_name].obsm['spatial']))
            nan_coords = np.unique(nan_coords[:,0])
            idx = ~adatas[f_name].obs.index.isin(adatas[f_name].obs.index[nan_coords])
            adatas[f_name] = adatas[f_name][idx,:]
        # adatas[f_name].obsm['spatial'] = adatas[f_name].obs[['array_col', 'array_row']]
        sq.gr.spatial_neighbors(adatas[f_name], n_rings=2, coord_type="grid") 
        ohe = OneHotEncoder()
        minor_one = ohe.fit_transform(adatas[f_name].obs[['histo_annotation_num']])
        adatas[f_name].obsm['nodules_ohe'] = minor_one
        minor_dist = adatas[f_name].obsp['spatial_connectivities'] @ minor_one 
        inside_spots = minor_dist.multiply(minor_one)
        outside_spots = (minor_dist.multiply((1-minor_one.toarray()))).toarray()
        max_values = np.max(inside_spots.toarray(), axis=0)
        max_indices = np.where(inside_spots.toarray() >= 15)
        labels = np.full_like(inside_spots.toarray(), 'outside', dtype='<U7')  # Assuming a string length of 7 is enough
        labels[inside_spots.toarray() != 0] = 'edge'
        labels[max_indices] = 'core'
        labels[np.where(outside_spots > 0)] = 'closure'
        adatas[f_name].uns['nodules'] = n_all = ohe.categories_[0].astype('str')
        adatas[f_name].obsm['nodule_compositions'] = np.array(pd.DataFrame(labels, index=adatas[f_name].obs_names, columns=n_all))
        adatas[f_name].obsm['nodules_core']    = pd.DataFrame((adatas[f_name].obsm['nodule_compositions'] == 'core') * 1, columns=n_all, index=adatas[f_name].obs_names)
        adatas[f_name].obsm['nodules_edge']    =  pd.DataFrame((adatas[f_name].obsm['nodule_compositions'] == 'edge') * 1, columns=n_all, index=adatas[f_name].obs_names)
        adatas[f_name].obsm['nodules_inside']  =  pd.DataFrame(((adatas[f_name].obsm['nodules_core'] + adatas[f_name].obsm['nodules_edge']) > 0) * 1, columns=n_all, index=adatas[f_name].obs_names)
        is_closure = (adatas[f_name].obsm['nodule_compositions'] == 'closure') * 1
        to_include = 1 - (np.array(adatas[f_name].obsm['nodules_inside'])[:,:-1].sum(1) > 0) * 1
        adatas[f_name].obsm['nodules_closure'] = pd.DataFrame(np.multiply(is_closure, to_include[:, np.newaxis]), columns=n_all, index=adatas[f_name].obs_names) 
        # adatas[f_name].write(SAVE_DIR / (f_name + '.h5ad'))

genes = list(set.intersection(*[set(adatas[f].var_names) for f in adatas.keys()]))
for f_name in adatas.keys():
    adatas[f_name] = adatas[f_name][:,genes]
    adatas[f_name].write(SAVE_DIR / (f_name + '.h5ad'))
    


