import matplotlib.pyplot as plt
import anndata
import scanpy as sc

def plot_annotations(adata, numerical_colum, exclude_mask=None, ax=None, border=40, lims=None, color_cycler=None):
    ad = adata.copy()
    
    if exclude_mask is not None:
        ad = ad[~exclude_mask]

    ## correct for plotting
    key = list(ad.uns['spatial'].keys())[0]
    ad = ad.copy()
    ad.obsm['spatial'] = ad.obsm['spatial'] * ad.uns['spatial'][key]['scalefactors']['tissue_lowres_scalef']
    ad.uns['spatial'][key]['scalefactors']['spot_diameter_fullres'] *= ad.uns['spatial'][key]['scalefactors']['tissue_lowres_scalef']
    if lims is None:
        y, x, _ = ad.uns['spatial'][key]['images']['lowres'].shape
    else:
        y, x = lims
    sc.pl.spatial(ad, img_key="lowres", color=numerical_colum, legend_loc='on data', scale_factor=1, ax=ax,
                  show=False,palette=color_cycler)
    
    if ax is None:
        ax = plt.gca()
    ax.set_xlim(border, x-border)
    ax.set_ylim(border, y-border)
    ax.invert_yaxis()
    # ax.set_title(k)
    ax.set_xlabel('')
    ax.set_ylabel('')
    # ax.set_rasterization_zorder(1)
