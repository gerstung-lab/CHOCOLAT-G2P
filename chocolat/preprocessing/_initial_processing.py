import numpy as np

def annotate_normal_cancer(adata, th=800, row_th=None, col_th=None):
    # Step 1: Identify cancer spots
    cancer_spots = np.array([not x.upper().startswith('NORM') for x in adata.obs['histo_annotation']])

    # Step 2: Calculate distances
    coordinates = adata.obsm['spatial']
    distances = np.sqrt((coordinates[:, np.newaxis] - coordinates)**2).sum(axis=2)

    # Step 3: Filter non-cancer spots
    non_cancer_spots = ~cancer_spots

    # Step 4: Determine distance criteria
    min_distance_to_cancer = distances[non_cancer_spots][:, cancer_spots].min(axis=1)
    far_from_cancer = min_distance_to_cancer >= th
    
    # Step 5: Put annotations into new column
    adata.obs['histo_decision'] = 'excluded'
    adata.obs.loc[non_cancer_spots, 'histo_decision'] = ['normal' if x else 'excluded' for x in far_from_cancer]
    if row_th:
        adata.obs.loc[adata.obs['array_row'].astype(int) < row_th, 'histo_decision'] = 'excluded'
    if col_th:
        adata.obs.loc[adata.obs['array_col'].astype(int) < col_th, 'histo_decision'] = 'excluded'
    adata.obs.loc[cancer_spots, 'histo_decision'] = 'cancer'
