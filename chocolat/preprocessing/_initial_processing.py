import numpy as np

def annotate_normal_cancer(adata, th=800, row_th=None, col_th=None, histo_column='histo_annotation'):
    """
    Annotate a spatial data object with 'normal', 'cancer', or 'excluded' based on histological
    annotations and spatial distances.

    Parameters
    ----------
    adata : AnnData
        An annotated data matrix with `.obs` containing histological annotations in 
        `histo_column` and `.obsm['spatial']` containing spatial coordinates.
    th : int, optional
        Threshold distance to define whether non-cancer spots are considered far from 
        cancer spots. Default is 800.
    row_th : int, optional
        Threshold for excluding spots based on their 'array_row' value. If None, this
        filter is not applied.
    col_th : int, optional
        Threshold for excluding spots based on their 'array_col' value. If None, this
        filter is not applied.
    histo_column : str, optional
        The column name in `adata.obs` that contains histological annotations. If 
        annotation start with norm, it will be treated as normal tissue. 
        Everything else will be treated as cancer. Default is 'histo_annotation'.

    Returns
    -------
    None
        Modifies `adata.obs` in-place by adding a 'histo_decision' column with the 
        annotations 'normal', 'cancer', or 'excluded'.
    """
    cancer_spots = np.array([not x.upper().startswith('NORM') for x in adata.obs[histo_column]])

    # Calculate distances
    coordinates = adata.obsm['spatial']
    distances = np.sqrt((coordinates[:, np.newaxis] - coordinates)**2).sum(axis=2)

    # Filter non-cancer spots
    non_cancer_spots = ~cancer_spots

    # Determine distance criteria
    min_distance_to_cancer = distances[non_cancer_spots][:, cancer_spots].min(axis=1)
    far_from_cancer = min_distance_to_cancer >= th
    
    # Put annotations into new column
    adata.obs['histo_decision'] = 'excluded'
    adata.obs.loc[non_cancer_spots, 'histo_decision'] = ['normal' if x else 'excluded' for x in far_from_cancer]
    if row_th:
        adata.obs.loc[adata.obs['array_row'].astype(int) < row_th, 'histo_decision'] = 'excluded'
    if col_th:
        adata.obs.loc[adata.obs['array_col'].astype(int) < col_th, 'histo_decision'] = 'excluded'
    adata.obs.loc[cancer_spots, 'histo_decision'] = 'cancer'
