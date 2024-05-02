import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import numpy as np

def plot_annotations(adata, numerical_colum, exclude_mask=None, ax=None, border=40, lims=None, color_cycler=None):
    """
    Plot spatial data region annotations of from an AnnData object, optionally excluding certain data points,
    adjusting image borders, and customising color scales.

    Parameters
    ----------
    adata : AnnData
        An annotated data matrix with spatial information under `.obsm['spatial']` and scale 
        factors under `.uns['spatial'][key]['scalefactors']`.
    numerical_column : str
        The name of the column in `adata.obs` with region ids, used to colour the data points.
    exclude_mask : array-like, optional
        Boolean array where True indicates observations to exclude from the plot. If None, all
        observations are excluded.
    ax : matplotlib.axes.Axes, optional
        The axes object where the plot will be drawn. If None, creates new plot.
    border : int, optional
        Margin size to apply around the plot within the axes. Default is 40.
    lims : tuple of int, optional
        A tuple (y, x) used to manually set the limits for the plot dimensions. If None, the
        limits are derived from the data.
    color_cycler : Cycler, optional
        Matplotlib cycler object for customising color cycles in the plot.

    Returns
    -------
    None
        Modifies the given or current axes object with the spatial data plot.
    """
    ad = adata.copy()
    
    if exclude_mask is not None:
        ad = ad[~exclude_mask]
        
    if ax is None:
        fig = plt.figure(figsize=((5, 5)))
        subfigs = fig.subfigures(1, 1)
        ax = subfigs.subplots(1, 1, sharey=True)

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

def plot_spatial_pertrubation_probabilities(adata_dict_filtered, samples_dict, region_ids_tensor, plasmids, keys=None):
    """
    Plot spatial perturbation probabilities in array view, for specified samples and plasmids.

    Parameters
    ----------
    adata_dict_filtered : dict of AnnData
        A dictionary of filtered AnnData objects indexed by sample identifiers, containing 
        spatial gene expression data.
    samples_dict : dict
        A dictionary containing samples from the inferred posterior distribution.
    region_ids_tensor : dict of tensors
        A dictionary of tensors encoding backwards mapping of regioins into spots.
    plasmids : list of str
        A list of plasmid names for which the perturbation probabilities are plotted.
    keys : list of str, optional
        A list of sample keys to include in the plot. If None, all keys from `samples_dict` are used.

    Returns
    -------
    None
        Generates a multi-panel plot, each panel showing spatial perturbation probability for a given plasmid 
        and sample.

    """

    if keys is None:
        keys = list(samples_dict.keys())
    
    fig = plt.figure(layout='constrained', figsize=((3 * len(plasmids), 3 * len(keys))))
    subfigs = fig.subfigures(len(keys), len(plasmids), wspace=0.0)

    # Convert the numpy array of Axes objects to a list

    for si, s in enumerate(keys):
        dims = samples_dict[s]['region_params'].mean(0).shape
        noise = np.ones([1, dims[1], dims[2]])
        noise[:,:,0] = 0.995
        noise[:,:,0] = 0.001
        region_params = np.concatenate([samples_dict[s]['region_params'].mean(0),noise])

        for i, n in enumerate(plasmids):
            adata_dict_filtered[s].obs[n] = region_params[region_ids_tensor[s].cpu().numpy(),:][:,i][:,1:].sum(-1)
        # fig, ax = plt.subplots(1,8, layout='constrained', figsize=(5*8,5))

        ad = adata_dict_filtered[s]
        ad = ad[~((ad.obs['histo_annotation_num'] == 'NaN'))]

        for i in range(len(plasmids)):
            ax = subfigs[si,i].subplots(1, 1, sharey=True)
            sc.pl.spatial(ad, img_key=None, color=plasmids[i], vmin=0, vmax=1,
                          cmap='RdBu_r', ax=ax, spot_size=220, show=False, colorbar_loc=None)
            # ax.set_rasterization_zorder(3)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # Remove x-axis ticks and labels
            ax.xaxis.set_ticks_position('none') 
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticks_position('none') 
            ax.yaxis.set_ticklabels([])
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')

            if i == 0:
                ax.set_ylabel(s)
            if si == 0:
                ax.set_title(plasmids[i])