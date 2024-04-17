import scanpy as sc
import squidpy as sq
import pandas as pd 
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder  
import pickle
from collections import defaultdict
import math



VISIUM_DATA_DIR = Path('data_tidy/composed_anndata_objects_2024_01_30')
GENO_DIR        = Path('data_tidy/geno_data')

GENE_SIG_F  = Path('data_tidy/2024-01-05_Modules_GeneSignatures.xlsx')
GENE_SEL_F  = Path('data_tidy/2024-01-17_MarkerGenes_to_include.xlsx')
GENO_PROB_F = GENO_DIR / 'plasmid_probs-12-samples-final.pkl'


# SAVE_DIR    = Path('data_tidy/geno_nodule_2024_01_15')
# os.makedirs(SAVE_DIR)
# OVERWRITE = False

with open(GENO_PROB_F, 'rb') as f:
    geno_data = pickle.load(f)
    
def calculate_z_scores(classes_matrix, features_matrix):
    class_means = np.array([np.mean(features_matrix[classes_matrix[:, i] == 1], axis=0) for i in range(classes_matrix.shape[1])])
    class_std_devs = np.array([np.std(features_matrix[classes_matrix[:, i] == 1], axis=0) for i in range(classes_matrix.shape[1])])
    z_scores = np.nan_to_num(class_means / 1, 0)
    return z_scores

gene_sig = pd.read_excel(GENE_SIG_F)
gene_sig['Gene'] = gene_sig['Gene'].str.capitalize()
gene_sig = gene_sig.set_index('Gene')
gene_sel = pd.read_excel(GENE_SEL_F)
gene_sel['MarkerGenes'] = gene_sel['MarkerGenes'].str.capitalize()
gene_sel['CoreMarker'] = gene_sel['CoreMarker'].str.capitalize()
gene_core = list(set(gene_sel['CoreMarker'].dropna().values))
pheno     = gene_sel.loc[~gene_sel['MarkerGenes'].isnull(), ['MarkerGenes', 'Phenotype']]
pheno = pheno.drop_duplicates('MarkerGenes', keep='first')
pheno = pheno.set_index('MarkerGenes')

gene_sel  = list(set(gene_sel['MarkerGenes'].dropna().values))

ANNOTS_DIR        = Path('data_raw/2023-11-24_AllAnnosFinal')
NODULE_ANNOT_DIR  = ANNOTS_DIR / 'Nodules'
mapping_reps   = pd.read_csv(( NODULE_ANNOT_DIR / 'mapping_replicates.csv'), sep=';')


# tme_genes = {
#     "Fibroblasts"  : ["Col1a1", "Col1a2", "Col3a1", "Igfbp5", "Pdgfrb", 'Dpt'],
#     "Erythrocytes" : ["Hbb−bs", "Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1"],
#     "Kupffer cells": ["Clec4f", "Csf1r", "Marco", "Cd163"],
#     "Neutrophils"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
#     "Mast cells"   : ["Cpa3", "Cma1", "Tpsb2", "Tpsab1"],
#     "NK cells"     : ["Itgax", "S100a4", "Gzmb", "Gzma", "Lat", "Cd27", "Rgs1"],
#     "Platelets"    : ["Pf4", "Itga2b", "Ppbp"],
#     "B cells"      : ["Cd79b", "Jchain", "Ighd", "Cd79a", "Pou2af1", "Cd19", 'Igkc'],
#     "T cells"      : ["Trbc1", "Zap70", "Cd4", "Cd8a", 'Trbc2', 'Ccl5'],
#     # "Monocytes"    : ["Cd163", "Mrc1", "Irf7", "Ptprc"],
#     "Endothelial"  : ["Fam167b", "Pecam1", "Tek"]
# }

# tumor_genes = {
#     # "Other"          : ["Serpina1c", "Fabp1", "Echs1", "Acly", "Asl", "Mup3", "G6pc", "Ass1", "G6pc"],
#     "Periportal HCC" : ["Sds", "Sdsl", "Hal"],
#     "Pericentral HCC": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2", "Gstm1", 'Nrn1', 'Axin2'],
#     "Cholangiocytes" : ["Krt19", "Cldn7", "Krt7", "Epcam", 'Mapk13']
# }

tme_genes = {
      "Fibroblasts"  : ["Col1a1", "Col3a1", "Col1a2", "Col5a1", "Dpt"],
      "Erythrocytes" : ["Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1"],
      "Kupffer cells": ["Clec4f", "Csf1r", "Marco", "Cd163"],
      "Neutrophils"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
      "Mast cells"   : ["Cpa3", "Cma1", "Tpsb2", "Mcpt4", 'Tpsab1'],
    #   "NK cells"     : ["Itgax", "S100a4", "Gzmb", "Gzma", "Lat", "Cd27", "Rgs1"],
      "Platelets"    : ["Pf4", "Itga2b", "Ppbp"],
      "B cells"      : ["Igkc", "Jchain"],
      "add_on"       : ["C1qc", "C1qb", "C1qa", "Ctss", "Ctsc", "Gas6", "Thbs1", "Col6a3", "Snca", "Hemgn", "Gata1", "Rhd", 'Sparc', 'Mmp2', 'Tagln', 'Thbs1', 'Postn'],
    
    #   "T cells"      : ["Trbc1", "Zap70", "Cd4", "Cd8a", 'Trbc2', 'Ccl5']
    # "Monocytes"    : ["Cd163", "Mrc1", "Irf7", "Ptprc"],
    #   "Endothelial"  : ["Vwf", "Klf4", "Eng", "Fam167b", "Pecam1", "Tek"]
}


tumor_genes = {
    "Periportal HCC": ["Sds", "Sdsl", "Hal"],
    "Pericentral HCC": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2"],
    "Cholangiocytes": ["Krt19", "Cldn7", "Krt7", "Mapk13", 'Epcam'],
    "add_on": ["Gp2", "Slc15a2", 'Kcnk1', "Elf3", "Gp2", "Ezr", "Ehf", "Gstm3", "Gstm2", "Axin2", 
               "Sox9", 'Spp1', 'Ccl9', 'Cyp2f2', 'Ly6d', 'Asl', 'Mup3', 'Serpina1c', 'Rida',
               'Ppp1r1b', 'Tmem229a', 'Bex1', 'Bex2', 'Upp2', 'Aspg', 'Serpina6', 'Hamp2', 
               'Ptger3', 'Bdh2', 'Car9', 'Tspan8', 'Pck1', 'Vil1', 'Car3', 'Lgals2'],
    "Histone": ["Hist1h3d","Hist1h2ag", "Hist1h1c", "Hist1h1e", "Hist1h4h", "Hist4h4", "Hist1h3c"]
    # ,
    # "Other": ["Upp2", "Vil1", 'Bpifb1', 'Serpina1e', 'Apoa4', 'Car3', 'Serpina6', 'Lgals2']
}



flattened_tumor = [value for sublist in tumor_genes.values() for value in sublist]

keys_6ROI = ['ML_III_A', 'ML_III_B', 'ML_II_A_1','ML_II_B', 'ML_II_C', 'ML_I_2']


adatas = {}
adatas_geno = {}
# for sample in list(set(geno_data['samples'])):
for sample in keys_6ROI:
    mapping_df   = mapping_reps.loc[mapping_reps['sample_name'] == sample,['old_annotation', 'new_annotation']]
    mapping_dict = mapping_df.set_index('old_annotation').to_dict()['new_annotation']
    print(sample)
    adatas[sample] = sc.read(VISIUM_DATA_DIR / (sample + '.h5ad'))
    adata          = adatas[sample]
    adata.obs['histo_annotation_num'] = [mapping_dict.get(item, 'NaN') for item in adata.obs['histo_annotation']]
    # genes          = list(set(matrix_df.index) & set(adata.var_names))
    # adata          = adata[:, genes]
    adata          = adata[:, adata.var['BC'].isna()]
    # sc.pp.filter_genes(adata, min_counts=)
    sc.pp.filter_genes(adata, max_counts=1e06)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # genes          = list(set(matrix_df.index) & set(adata.var_names))
    sc.pp.highly_variable_genes(adata, n_top_genes=15000, flavor='seurat_v3')
    hvgs = adata.var_names[adata.var['highly_variable']]   
    genes = list((set(hvgs) | set(gene_core)|  set(tme_genes) | set(tumor_genes) | set(add_on_genes))& set(adata.var_names))
    adata = adata[:,genes]
    # sc.pp.normalize_total(adata, target_sum=1e4)
    sig_genes = list(set(adata.var_names) & set(gene_sig.index))
    gene_sig_copy = gene_sig.loc[sig_genes,:]
    out_genes = list(set(adata.var_names) - set(gene_sig_copy.index))
    gene_sel = list(set(adata.var_names) & set(gene_sel))
    gene_core = list(set(adata.var_names) & set(gene_core))
    sig_mask = pd.get_dummies(gene_sig_copy['Module'].astype(str)).groupby(level=0).max().sort_index()
    sig_mask = pd.concat([sig_mask, pd.DataFrame(0, index=out_genes, columns=sig_mask.columns)])
    # assert adata.obs['histo_annotation']
    ind            = np.where(np.array(geno_data['samples']) == sample)[0]
    nodules        = np.array(geno_data['names'])[ind]
    plasmid_mask   = pd.DataFrame(geno_data['probs'][ind,:].reshape([len(ind), 8]), columns=geno_data['plasmids'])
    to_add_nodules = list(set(adata.obs['histo_annotation_num']) - set(nodules) - set(['NaN']))
    to_add_nodules = [x for x in to_add_nodules if not isinstance(x, float) or not math.isnan(x)]
    plasmid_normal = pd.DataFrame(np.zeros([len(to_add_nodules), 8]), columns=geno_data['plasmids'])
    plasmid_mask = pd.concat([plasmid_mask, plasmid_normal])
    nodules        = list(nodules) + to_add_nodules
    plasmid_mask.index = nodules
    inside_scores    = np.array(adata.obsm['nodules_inside'].loc[:,nodules])
    z_scores         = np.log(calculate_z_scores(inside_scores, adata.X.toarray()) + 1)
    geno_inside_exp  = pd.DataFrame(z_scores, index=nodules)
    geno_pheno_adata = sc.AnnData(geno_inside_exp, obsm={'plasmid_mask': plasmid_mask}, var=adata[:,genes].var, varm={'sig_mask': sig_mask.loc[genes,:]})
    geno_pheno_adata.var['sig'] = 'NaN'
    geno_pheno_adata.var.loc[sig_genes,'sig'] = geno_pheno_adata.varm['sig_mask'].loc[sig_genes,:].idxmax(1).values
    geno_pheno_adata.var['sel'] = 'F'
    geno_pheno_adata.var.loc[gene_sel,'sel'] = 'T'
    geno_pheno_adata.var['core'] = 'F'
    geno_pheno_adata.var.loc[gene_core,'core'] = 'T'
    geno_pheno_adata.var['pheno'] = 'NaN'
    geno_pheno_adata.var.loc[gene_sel,'pheno'] = pheno.loc[gene_sel,'Phenotype'].values
    geno_pheno_adata.var['tme'] = 'NaN'
    for tme, marker_genes in tme_genes.items():
        geno_pheno_adata.var.loc[geno_pheno_adata.var_names.intersection(marker_genes), 'tme'] = tme
    # geno_pheno_adata.var['add_on'] = 'NaN'
    # for add_on, marker_genes in add_on_genes.items():
    #     geno_pheno_adata.var.loc[geno_pheno_adata.var_names.intersection(marker_genes), 'add_on'] = add_on
    geno_pheno_adata.var['tumor'] = 'NaN'
    for tumor, marker_genes in tumor_genes.items():
        geno_pheno_adata.var.loc[geno_pheno_adata.var_names.intersection(marker_genes), 'tumor'] = tumor
    for ct in tumor_genes.keys():
        adata_sub = geno_pheno_adata.copy()[:, tumor_genes[ct]]
        norm_X = adata_sub.X.toarray()  / adata_sub.X.toarray().sum(axis=0)
        # norm_X = stats.zscore(norm_X, 0)
        geno_pheno_adata.obs[ct] = norm_X.mean(1)
    geno_pheno_adata.obs['dom_type'] = geno_pheno_adata.obs.loc[:, tumor_genes.keys()].idxmax(axis=1)
    geno_pheno_adata.obs['conf_type'] = geno_pheno_adata.obs.loc[:, tumor_genes.keys()].max(axis=1)
    geno_pheno_adata.obs.loc[geno_pheno_adata.obs['conf_type'] < .01,'dom_type']= 'other'
    # core
    core_scores   = np.array(adata.obsm['nodules_core'].loc[:,nodules])
    z_scores      = np.log(calculate_z_scores(core_scores, adata.X.toarray()) + 1)
    geno_pheno_adata.layers['core'] = z_scores
    # edge
    edge_scores   = np.array(adata.obsm['nodules_edge'].loc[:,nodules])
    z_scores      = np.log(calculate_z_scores(edge_scores, adata.X.toarray()) + 1)
    geno_edge_exp = z_scores
    geno_pheno_adata.layers['edge'] = z_scores
    # closure
    closure_scores   = np.array(adata.obsm['nodules_closure'].loc[:,nodules])
    z_scores         = np.log(calculate_z_scores(closure_scores, adata.X.toarray()) + 1)
    geno_closure_exp = z_scores
    geno_pheno_adata.layers['closure'] = z_scores
    # geno_pheno_adata.obs['expected_number'] = probs.sum(1)
    geno_pheno_adata.obs['CCC_score'] = calculate_z_scores(inside_scores, adata.obs['CCC_type']) 
    geno_pheno_adata.obs['CCC_score'].fillna(0, inplace=True)
    geno_pheno_adata.obs['SDS_score'] = calculate_z_scores(inside_scores, adata.obs['SDS_type']) 
    geno_pheno_adata.obs['SDS_score'].fillna(0, inplace=True)
    for stain in ['GFP_stain', 'GS_stain', 'RFP_stain']:
        if stain in adata.obs.columns:
            adata.obs[stain] = adata.obs[stain].astype('object')
            # Replace values
            adata.obs[stain] = adata.obs[stain].replace({np.nan: 0, 'low': 1, 'moderate': 2, 'hi': 3})
            z_scores         = calculate_z_scores(inside_scores, adata.obs[stain]) 
            geno_pheno_adata.obs['inside_'+stain]  = (z_scores > 0 ) * 1
            geno_pheno_adata.obs['inside_'+stain].fillna(0, inplace=True)
            z_scores         = calculate_z_scores(core_scores, adata.obs[stain]) 
            geno_pheno_adata.obs['core_'+stain]  = (z_scores > 0 ) * 1
            geno_pheno_adata.obs['core_'+stain].fillna(0, inplace=True)
            z_scores         = calculate_z_scores(edge_scores, adata.obs[stain]) 
            geno_pheno_adata.obs['edge_'+stain]  = (z_scores > 0 ) * 1
            geno_pheno_adata.obs['edge_'+stain].fillna(0, inplace=True)
            z_scores         = calculate_z_scores(closure_scores, adata.obs[stain]) 
            geno_pheno_adata.obs['closure_'+stain]  = (z_scores > 0 ) * 1
            geno_pheno_adata.obs['closure_'+stain].fillna(0, inplace=True)
    adatas_geno[sample] = geno_pheno_adata
    # nodules_all = [mapping_dict.get(x, 'NaN') for x in nodules_all]
    adatas_geno[sample].obs_names = [sample + '::' + str(s) if not (isinstance(s, float) and math.isnan(s)) else sample + '::' + 'NaN' for s in nodules]
        # Convert back to Categorical
        # adata.obs[stain] = pd.Categorical(adata.obs[stain], categories=[0, 1, 2, 3], ordered=False)
    # geno_pheno_adata.write(SAVE_DIR / ('geno_pheno_' + sample + '.h5ad'))
    
genes = list(set.intersection(*[set(adatas_geno[f].var_names) for f in adatas_geno.keys()]))
for f_name in adatas_geno.keys():
    adatas_geno[f_name] = adatas_geno[f_name][:,genes]
    
import anndata as ad
adatas_joint = ad.concat(adatas_geno, label='batch')
sc.pp.filter_genes(adatas_joint, min_counts=10)
sc.pp.neighbors(adatas_joint)
sc.tl.umap(adatas_joint)
sc.pl.umap(adatas_joint, color='batch', save='_batch_effects.pdf')
sc.pp.combat(adatas_joint, key='batch')
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, save='_batch_effects_removed.pdf')

adatas_joint.var = adatas_geno['ML_I_2'].var
adatas_joint.varm = adatas_geno['ML_I_2'].varm
adatas_joint.write('data_tidy/new_geno_nodule_2024_02_19.h5ad')


add = adatas_joint[:, flattened_tumor].copy()
sc.pp.neighbors(add)
sc.tl.leiden(adatas_joint)





