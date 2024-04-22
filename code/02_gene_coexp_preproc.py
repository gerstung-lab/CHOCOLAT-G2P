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
GENO_PROB_F = GENO_DIR / 'genotypes-6-samples.pkl'


SAVE_DIR    = Path('data_tidy/geno_nodule_2024_01_15')
os.makedirs(SAVE_DIR)
OVERWRITE = False

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



tme_genes = {
      "Fibroblasts"  : ["Col1a1", "Col3a1", "Col1a2", "Col5a1"],
      "Erythrocytes" : ["Hbb-bt", "Hba-a2", "Alas2", "Tmcc2", "Slc4a1"],
      "Kupffer cells": ["Clec4f", "Csf1r", "Marco", "Cd163"],
      "Neutrophils"  : ["S100a9", "S100a8", "Ngp", "Ltf", "Camp", "Elane", "Ctsg", "Mpo"],
      "Mast cells"   : ["Cpa3", "Cma1", "Tpsb2", "Mcpt4", 'Tpsab1'],
    #   "NK cells"     : ["Itgax", "S100a4", "Gzmb", "Gzma", "Lat", "Cd27", "Rgs1"],
      "Platelets"    : ["Pf4", "Itga2b", "Ppbp"],
      "B cells"      : ["Igkc", "Jchain"],
      "tme_add_on"       : [ "C1qa", "Ctss", "Ctsc", "Gas6", "Thbs1", "Hemgn", "Gata1", "Rhd", 'Sparc', 'Mmp2', 'Thbs1'],
}


tumor_genes = {
    "Periportal HCC": ["Sds", "Sdsl", "Hal"],
    "Pericentral HCC": ["Oat", "Gulo", "Cyp2e1", "Cyp1a2"],
    "Cholangiocytes": ["Krt19", "Cldn7", "Krt7", "Mapk13", 'Epcam'],
    "tumor_add_on": ["Gp2", "Slc15a2", "Gp2", "Ezr", "Ehf", "Gstm3", "Gstm2", 
                'Upp2', 'Aspg', 'Hamp2', 
               'Bdh2', 'Tspan8', 'Vil1', 'Car3'],
    'Other': ['Hamp2', 'Vil1', 'Car3', 'Upp2'],
    "Histone": ["Hist1h3d", "Hist1h1c", "Hist1h1e", "Hist1h4h", "Hist1h3c"]
    # "Other": ["Upp2", "Vil1", 'Bpifb1', 'Serpina1e', 'Apoa4', 'Car3', 'Serpina6', 'Lgals2']
}


adatas = {}
for sample in os.listdir(VISIUM_DATA_DIR):
    sample = sample.rstrip(".h5ad")
    if sample != 'ML_II_B_2Cyt':
        mapping_df   = mapping_reps.loc[mapping_reps['sample_name'] == sample,['old_annotation', 'new_annotation']]
        mapping_dict = mapping_df.set_index('old_annotation').to_dict()['new_annotation']
        print(sample)
        adatas[sample] = sc.read(VISIUM_DATA_DIR / (sample + '.h5ad'))
        adata          = adatas[sample]
        adata.obs['histo_annotation_num'] = [mapping_dict.get(item, 'NaN') for item in adata.obs['histo_annotation']]
        sc.pp.filter_genes(adata, max_counts=1e06)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # genes          = list(set(matrix_df.index) & set(adata.var_names))
        sc.pp.highly_variable_genes(adata, n_top_genes=10000, flavor='seurat_v3')
        hvgs = adata.var_names[adata.var['highly_variable']]   
        genes = list((set(hvgs) | (set(gene_core) & set(adata.var_names)) |  {x for v in tme_genes.values() for x in v} | {x for v in tumor_genes.values() for x in v} | set(adata.var_names[~adata.var['BC'].isna()]) | (set(adata.var_names) & set(gene_sig.index))) & set(adata.var_names))
        adata = adata[:,genes]
        sc.pp.normalize_total(adata, target_sum=1e4)
        sig_genes = list(set(adata.var_names) & set(gene_sig.index))
        gene_sig_copy = gene_sig.loc[sig_genes,:]
        out_genes = list(set(adata.var_names) - set(gene_sig_copy.index))
        gene_sel = list(set(adata.var_names) & set(gene_sel))
        gene_core = list(set(adata.var_names) & set(gene_core))
        adata.var['tme'] = 'NaN'
        for tme, marker_genes in tme_genes.items():
            adata.var.loc[adata.var_names.intersection(marker_genes), 'tme'] = tme
        adata.var['tumor'] = 'NaN'
        for tumor, marker_genes in tumor_genes.items():
            adata.var.loc[adata.var_names.intersection(marker_genes), 'tumor'] = tumor
        adatas[sample] = adata
    
genes = list(set.intersection(*[set(adatas[f].var_names) for f in adatas.keys()]))
for f_name in adatas.keys():
    adatas[f_name] = adatas[f_name][:,genes]
    
import anndata as ad
adatas_joint = ad.concat(adatas, label='batch')

adatas_joint.var = adatas['ML_I_2'].var.loc[genes,:]
adatas_joint.varm = adatas['ML_I_2'][:,genes].varm


adatas_joint.obs['in_tissue'] = adatas_joint.obs['in_tissue'].astype('int')
adatas_joint.obs['array_row'] = adatas_joint.obs['array_row'].astype('int')
adatas_joint.obs['array_col'] = adatas_joint.obs['array_col'].astype('int')
adatas_joint.obsm['spatial'] = adatas_joint.obsm['spatial'].astype('float')
adatas_joint.write('data_tidy/joint_coexp_2024_02_29.h5ad')