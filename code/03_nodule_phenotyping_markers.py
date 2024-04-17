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
    
def calculate_means(classes_matrix, features_matrix):
    class_means = np.array([np.mean(features_matrix[classes_matrix[:, i] == 1], axis=0) for i in range(classes_matrix.shape[1])])
    z_scores = np.nan_to_num(class_means / 1, 0)
    return z_scores

def calculate_quantiles(classes_matrix, features_matrix):
    class_means = np.array([np.quantile(features_matrix[classes_matrix[:, i] == 1], q=.95, axis=0) for i in range(classes_matrix.shape[1])])
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


# B cells —> B cell
# erythroblast —> erythroblast
# fibroblast —> fibroblast
# kupffer cell/macrophage —> kupffer cell/macrophage
# mast cell —> mast cell
# neutrophil —> neutrophil
# platelet —> platelet
# histone-enriched —> histone-enriched
# Midlobular/Other —> Hamp2/Upp2-enriched
# central —> central
# portal —> portal


tme_genes = {
      "fibroblast"  : ["Col1a1", "Col3a1", "Col1a2", "Col5a1"],
      "erythroblast" : ["Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1"],
      "kupffer cell/macrophage": ["Clec4f", "Csf1r", "Marco", "Cd163"],
      "neutrophil"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
      "mast cell"   : ["Cpa3", "Cma1", "Tpsb2", "Mcpt4", 'Tpsab1'],
    #   "NK cells"     : ["Itgax", "S100a4", "Gzmb", "Gzma", "Lat", "Cd27", "Rgs1"],
      "platelet"    : ["Pf4", "Itga2b", "Ppbp"],
      "B cell"      : ["Igkc", "Jchain"],
      "tme_add_on"       : [ "C1qa", "Ctss", "Ctsc", "Gas6", "Thbs1", "Hemgn", "Gata1", "Rhd", 'Sparc', 'Mmp2', 'Thbs1', 'Ighj1'],
        
    #   "T cells"      : ["Trbc1", "Zap70", "Cd4", "Cd8a", 'Trbc2', 'Ccl5']
    # "Monocytes"    : ["Cd163", "Mrc1", "Irf7", "Ptprc"],
    #   "Endothelial"  : ["Vwf", "Klf4", "Eng", "Fam167b", "Pecam1", "Tek"]
}


tumor_genes = {
    "portal": ["Sds", "Sdsl", "Hal"],
    "central": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2"],
    "cholangiocytic": ["Krt19", "Cldn7", "Krt7", "Mapk13", 'Epcam'],
    "tumor_add_on": ["Gp2", "Slc15a2", "Gp2", "Ezr", "Ehf", "Gstm3", "Gstm2", 
                'Upp2', 'Aspg', 'Hamp2', 
               'Bdh2', 'Tspan8', 'Vil1', 'Car3'],
    'Hamp2/Upp2-enriched': ['Hamp2', 'Hamp', 'Car3', 'Upp2'],
    "histone-enriched": ["Hist1h3d", "Hist1h1c", "Hist1h1e", "Hist1h4h", "Hist1h3c"]
    # "Hamp2/Upp2-enriched": ["Upp2", "Vil1", 'Bpifb1', 'Serpina1e', 'Apoa4', 'Car3', 'Serpina6', 'Lgals2']
}



tme_genes_ext = {
      "fibroblast"  : ["Col1a1", "Col3a1", "Col1a2", "Col5a1", "Dpt", "Gas6", "Thbs1", "Col6a3", 'Sparc', 'Mmp2', 'Tagln', 'Thbs1', 'Postn'],
      "erythroblast" : ["Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1", "Snca", "Hemgn", "Gata1", "Rhd"],
      "kupffer cell/macrophage": ["Clec4f", "Csf1r", "Marco", "Cd163", "C1qc", "C1qb", "C1qa", "Ctss", "Ctsc"],
      "neutrophil"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
      "mast cell"   : ["Cpa3", "Cma1", "Tpsb2", "Mcpt4", 'Tpsab1'],
    #   "NK cells"     : ["Itgax", "S100a4", "Gzmb", "Gzma", "Lat", "Cd27", "Rgs1"],
      "platelet"    : ["Pf4", "Itga2b", "Ppbp"],
      "B cell"      : ["Igkc", "Jchain", 'Ighj1']
    
    #   "T cells"      : ["Trbc1", "Zap70", "Cd4", "Cd8a", 'Trbc2', 'Ccl5']
    # "Monocytes"    : ["Cd163", "Mrc1", "Irf7", "Ptprc"],
    #   "Endothelial"  : ["Vwf", "Klf4", "Eng", "Fam167b", "Pecam1", "Tek"]
}


tumor_genes_ext = {
    "portal": ["Sds", "Sdsl", "Hal", 'Bex1', 'Bex2', 
                    'Upp2', 'Aspg', 'Serpina6', 'Hamp2', 
                    'Ptger3', 'Bdh2'],
    "central": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2", "Gstm3", "Gstm2", "Axin2"],
    "cholangiocytic": ["Krt19", "Cldn7", "Krt7", "Mapk13", 'Epcam', "Gp2", 
                       "Slc15a2", 'Kcnk1', "Gp2", "Ezr", 
                       "Sox9", 'Spp1', 'Ppp1r1b', 'Car9', 'Tspan8', 'Mup3', 'Serpina1c', 'Lgals2'],
    "Hamp2/Upp2-enriched": ['Hamp2', 'Hamp', 'Car3', 'Upp2'],
    "histone-enriched": ["Hist1h3d","Hist1h2ag", "Hist1h1c", "Hist1h1e", "Hist1h4h", "Hist4h4", "Hist1h3c"]
    # ,
    # "Hamp2/Upp2-enriched": ["Upp2", "Vil1", 'Bpifb1', 'Serpina1e', 'Apoa4', 'Car3', 'Serpina6', 'Lgals2']
}



flattened_tumor = [value for sublist in tumor_genes.values() for value in sublist]

# keys_6ROI = ['ML_III_A', 'ML_III_B', 'ML_II_A_1','ML_II_B', 'ML_II_C', 'ML_I_2']
all_rois = ['ML_I_2', 'ML_II_B', 'ML_II_C', 'ML_III_A', 'ML_III_B', 'ML_II_A_1'] + ['ML_II_B_3Cyt', 'ML_II_A_3Cyt', 'ML_III_A_2Cyt', 'ML_I',
       'ML_II_A_2']

adatas = {}
adatas_geno = {}
# for sample in list(set(geno_data['samples']) - set(['ML_II_B_2Cyt'])):
for sample in all_rois:
    mapping_df   = mapping_reps.loc[mapping_reps['sample_name'] == sample,['old_annotation', 'new_annotation']]
    mapping_dict = mapping_df.set_index('old_annotation').to_dict()['new_annotation']
    print(sample)
    adatas[sample] = sc.read(VISIUM_DATA_DIR / (sample + '.h5ad'))
    adata          = adatas[sample]
    adata.obs['histo_annotation_num'] = [mapping_dict.get(item, 'NaN') for item in adata.obs['histo_annotation']]
    # inner_nodules = [mapping_dict.get(item, 'NaN') for item in adata.uns['nodules'][0]]
    inner_nodules = adata.uns['nodules']
    adata.obsm['nodules_inside'] = pd.DataFrame(adata.obsm['nodules_inside'], columns = inner_nodules, index = adata.obs_names)
    adata.obsm['nodules_closure'] = pd.DataFrame(adata.obsm['nodules_closure'], columns = inner_nodules, index = adata.obs_names)
    # genes          = list(set(matrix_df.index) & set(adata.var_names))
    # adata          = adata[:, genes]
    adata          = adata[:, adata.var['BC'].isna()]
    # sc.pp.filter_genes(adata, min_counts=)
    # sc.pp.filter_genes(adata, max_counts=1e06)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # genes          = list(set(matrix_df.index) & set(adata.var_names))
    sc.pp.highly_variable_genes(adata, n_top_genes=15000, flavor='seurat_v3')
    hvgs = adata.var_names[adata.var['highly_variable']]   
    genes = list((set(hvgs) | set(gene_core)|  set(tme_genes) | set(tumor_genes))& set(adata.var_names))
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
    plasmid_mask   = pd.DataFrame(geno_data['probs'][1, ind,:].reshape([len(ind), 8]), columns=geno_data['plasmids'])
    to_add_nodules = list(set(adata.obs['histo_annotation_num']) - set(nodules) - set(['NaN']))
    to_add_nodules = [x for x in to_add_nodules if not isinstance(x, float) or not math.isnan(x)]
    plasmid_normal = pd.DataFrame(np.zeros([len(to_add_nodules), 8]), columns=geno_data['plasmids'])
    plasmid_mask = pd.concat([plasmid_mask, plasmid_normal])
    nodules        = list(nodules) + to_add_nodules
    plasmid_mask.index = nodules
    inside_scores    = np.array(adata.obsm['nodules_inside'].loc[:,nodules])
    nodules          = list(np.array(nodules)[inside_scores.sum(0) != 0])
    inside_scores    = np.array(adata.obsm['nodules_inside'].loc[:,nodules])
    # closure_scores      = np.array(adata.obsm['nodules_closure'].loc[:,nodules])
    # inside_scores    = inside_scores + closure_scores
    z_scores         = np.log(calculate_means(inside_scores, adata.X.toarray()) + 1)
    geno_inside_exp  = pd.DataFrame(z_scores, index=nodules)
    plasmid_mask = plasmid_mask.loc[nodules,:]
    geno_pheno_adata = sc.AnnData(geno_inside_exp, obsm={'plasmid_mask': plasmid_mask}, var=adata[:,genes].var, varm={'sig_mask': sig_mask.loc[genes,:]})
    z_scores         = np.log(calculate_quantiles(inside_scores, adata.X.toarray()) + 1)
    geno_pheno_adata.layers['quantiles'] = z_scores
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
    for ct in tumor_genes_ext.keys():
        adata_sub = geno_pheno_adata.copy()[:, list(set(geno_pheno_adata.var_names) & set(tumor_genes_ext[ct]))]
        norm_X = adata_sub.X.toarray()  / adata_sub.X.toarray().sum(axis=0)
        # norm_X = stats.zscore(norm_X, 0)
        geno_pheno_adata.obs[ct + '_exp'] = norm_X.mean(1)
    # geno_pheno_adata.obs['dom_type'] = geno_pheno_adata.obs.loc[:, tumor_genes_ext.keys()].idxmax(axis=1)
    # geno_pheno_adata.obs['conf_type'] = geno_pheno_adata.obs.loc[:, tumor_genes_ext.keys()].max(axis=1)
    for ct in tme_genes_ext.keys():
        adata_sub = geno_pheno_adata.copy()[:, list(set(geno_pheno_adata.var_names) & set(tme_genes_ext[ct]))]
        norm_X = adata_sub.X.toarray()  / adata_sub.X.toarray().sum(axis=0)
        # norm_X = stats.zscore(norm_X, 0)
        geno_pheno_adata.obs[ct + '_exp'] = norm_X.mean(1)
    # geno_pheno_adata.obs['dom_tme'] = geno_pheno_adata.obs.loc[:, tme_genes_ext.keys()].idxmax(axis=1)
    # geno_pheno_adata.obs['conf_tme'] = geno_pheno_adata.obs.loc[:, tme_genes_ext.keys()].max(axis=1)
    # geno_pheno_adata.obs.loc[geno_pheno_adata.obs['conf_tme'] < .01,'dom_tme']= 'Hamp2/Upp2-enriched'
    # core
    # core_scores   = np.array(adata.obsm['nodules_core'].loc[:,nodules])
    # z_scores      = np.log(calculate_z_scores(core_scores, adata.X.toarray()) + 1)
    # geno_pheno_adata.layers['core'] = z_scores
    # # edge
    # edge_scores   = np.array(adata.obsm['nodules_edge'].loc[:,nodules])
    # z_scores      = np.log(calculate_z_scores(edge_scores, adata.X.toarray()) + 1)
    # geno_edge_exp = z_scores
    # geno_pheno_adata.layers['edge'] = z_scores
    # # closure
    # closure_scores   = np.array(adata.obsm['nodules_closure'].loc[:,nodules])
    # z_scores         = np.log(calculate_z_scores(closure_scores, adata.X.toarray()) + 1)
    # geno_closure_exp = z_scores
    # geno_pheno_adata.layers['closure'] = z_scores
    # # geno_pheno_adata.obs['expected_number'] = probs.sum(1)
    # geno_pheno_adata.obs['CCC_score'] = calculate_z_scores(inside_scores, adata.obs['CCC_type']) 
    # geno_pheno_adata.obs['CCC_score'].fillna(0, inplace=True)
    # geno_pheno_adata.obs['SDS_score'] = calculate_z_scores(inside_scores, adata.obs['SDS_type']) 
    # geno_pheno_adata.obs['SDS_score'].fillna(0, inplace=True)
    # for stain in ['GFP_stain', 'GS_stain', 'RFP_stain']:
    #     if stain in adata.obs.columns:
    #         adata.obs[stain] = adata.obs[stain].astype('object')
    #         # Replace values
    #         adata.obs[stain] = adata.obs[stain].replace({np.nan: 0, 'low': 1, 'moderate': 2, 'hi': 3})
    #         z_scores         = calculate_z_scores(inside_scores, adata.obs[stain]) 
    #         geno_pheno_adata.obs['inside_'+stain]  = (z_scores > 0 ) * 1
    #         geno_pheno_adata.obs['inside_'+stain].fillna(0, inplace=True)
    #         z_scores         = calculate_z_scores(core_scores, adata.obs[stain]) 
    #         geno_pheno_adata.obs['core_'+stain]  = (z_scores > 0 ) * 1
    #         geno_pheno_adata.obs['core_'+stain].fillna(0, inplace=True)
    #         z_scores         = calculate_z_scores(edge_scores, adata.obs[stain]) 
    #         geno_pheno_adata.obs['edge_'+stain]  = (z_scores > 0 ) * 1
    #         geno_pheno_adata.obs['edge_'+stain].fillna(0, inplace=True)
    #         z_scores         = calculate_z_scores(closure_scores, adata.obs[stain]) 
    #         geno_pheno_adata.obs['closure_'+stain]  = (z_scores > 0 ) * 1
    #         geno_pheno_adata.obs['closure_'+stain].fillna(0, inplace=True)
    adatas_geno[sample] = geno_pheno_adata
    adatas_geno[sample].obs.loc[:,'sample'] = sample
    adatas_geno[sample].obs.loc[:,'number'] = nodules
    # nodules_all = [mapping_dict.get(x, 'NaN') for x in nodules_all]
    adatas_geno[sample].obs_names = [sample + '::' + str(s) if not (isinstance(s, float) and math.isnan(s)) else sample + '::' + 'NaN' for s in nodules]
        # Convert back to Categorical
        # adata.obs[stain] = pd.Categorical(adata.obs[stain], categories=[0, 1, 2, 3], ordered=False)
    # geno_pheno_adata.write(SAVE_DIR / ('geno_pheno_' + sample + '.h5ad'))

def normalize_array(arr):
    # Calculate the quantiles
    q25 = np.quantile(arr, 0.25, axis=0)
    q99 = np.quantile(arr, 0.99, axis=0)
    # Clip values below 0.25 quantile to 0
    # arr[arr < q25] = 0
    # # Clip values above 0.99 quantile to 1
    # arr[arr >= q99] = 1
    # Normalize between 0 and 1
    normalized_arr = (arr - q25) / (q99 - q25)
    # Clip values to ensure they are between 0 and 1
    normalized_arr = np.clip(normalized_arr, 0, 1)
    return normalized_arr

genes = list(set.intersection(*[set(adatas_geno[f].var_names) for f in adatas_geno.keys()]))
for f_name in adatas_geno.keys():
    adatas_geno[f_name] = adatas_geno[f_name][:,genes]
    
import anndata as ad
adatas_joint = ad.concat(adatas_geno, label='batch')


adatas_joint.layers['norm_tme'] = normalize_array(np.power(adatas_joint.layers['quantiles'],1))
for gene_set in tme_genes.keys():
    gene_s = list(set(tme_genes[gene_set]) & set(adatas_joint.var_names))
    ind = adatas_joint.obs_names[np.where(np.sum(adatas_joint[:,gene_s].layers['norm_tme']> 0.6, 1) > 1)[0]]
    adatas_joint.obs[[(gene_set + '_bin')]] = 0
    adatas_joint.obs.loc[ind, (gene_set+ '_bin')] = 1
    


adatas_joint.obs['dom_tme'] = adatas_joint.obs.loc[:,[s + '_bin' for s in tme_genes.keys()]].idxmax(axis=1)
    

adatas_joint.layers['norm_tumor'] = normalize_array(np.power(adatas_joint.X,1))
for gene_set in tumor_genes.keys():
    if gene_set != 'tumor_add_on':
        gene_s = list(set(tumor_genes[gene_set]) & set(adatas_joint.var_names))
        # if gene_set != 'portal':
        ind = adatas_joint.obs_names[np.where(np.sum(adatas_joint[:,gene_s].layers['norm_tumor']> 0.5, 1) > 1)[0]]
        # else:
        #     ind = adatas_joint.obs_names[np.where(np.sum(adatas_joint[:,gene_s].layers['norm_tumor']> 0.6, 1) > 1)[0]]
        adatas_joint.obs[[(gene_set + '_bin')]] = 0
        adatas_joint.obs.loc[ind, (gene_set+ '_bin')] = 1
        adatas_joint.obs.loc[:, gene_set] = np.mean(adatas_joint[:,gene_s].layers['norm_tumor'], 1)
        adatas_joint.obs.loc[:, gene_set] = adatas_joint.obs.loc[:, gene_set] * adatas_joint.obs.loc[:, (gene_set+ '_bin')]


tumors = ['central', 'cholangiocytic', 'portal','Hamp2/Upp2-enriched', 'histone-enriched']
adatas_joint.obs['dom_tumor'] = adatas_joint.obs.loc[:, [s + '_bin' for s in tumors]].idxmax(axis=1)
# adatas_joint.obs.loc[adatas_joint.obs.loc[:, [s + '_bin' for s in tumors]].sum(1) == 0, 'dom_tumor'] = 'Unclassified'
adatas_joint.obs['dom_tumor'] = adatas_joint.obs['dom_tumor'].str.replace('_bin', '', regex=True)

tumors = ['central', 'cholangiocytic', 'portal','Hamp2/Upp2-enriched', 'histone-enriched']
adatas_joint.obs.loc[adatas_joint.obs['dom_tumor'] != 'histone-enriched', 'dom_tumor'] = adatas_joint.obs.loc[adatas_joint.obs['dom_tumor'] != 'histone-enriched', [s  for s in list(set(tumors) - set(['histone-enriched']))]].idxmax(axis=1)
adatas_joint.obs.loc[adatas_joint.obs.loc[:, [s  for s in tumors]].sum(1) == 0, 'dom_tumor'] = 'Unclassified'
adatas_joint.obs['dom_tumor'] = adatas_joint.obs['dom_tumor'].str.replace('_bin', '', regex=True)

adatas_joint.obs.loc[adatas_joint.obs['number'].str.contains("Normal"), 'dom_tumor'] = 'Normal'

adatas_joint.obs.to_csv('output/clustred_nodules_q_.6_2024_04_03.csv')
adatas_joint.var = adatas_geno['ML_I_2'].var
adatas_joint.varm = adatas_geno['ML_I_2'].varm
adatas_joint.write('data_tidy/new_geno_nodule_2024_04_03.h5ad')
