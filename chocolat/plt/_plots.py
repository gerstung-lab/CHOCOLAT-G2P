import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import numpy as np
from chocolat.preprocessing import compute_quantiles

### Order the reporters
reporters_order = [
    "Olfr103", "Tas2r102", "Vmn1r1", "Olfr1018", "Tas2r118", "Vmn1r174", "Olfr1", 
    "Olfr1000", "Tas2r103", "Vmn1r178", "Olfr1019", "Tas2r119", "Vmn1r175", 
    "Olfr1002", "Tas2r104", "Vmn1r12", 
    "Olfr1012", "Tas2r109", "Vmn1r169",
    "Olfr1006", "Tas2r105", "Vmn1r139", 
    "Olfr1008", "Tas2r106", "Vmn1r157", 
    "Olfr1009", "Tas2r107", "Vmn1r167", "Olfr107", "Olfr1014", "Tas2r113", "Vmn1r171", 
    "Olfr1013", "Tas2r110", "Vmn1r170", "Olfr1015", "Tas2r114", "Vmn1r172", 
]

### Correspondence between colours and plasmids
colours = {
    "MYC-plasmid": "#56b4e9",
    "coVEGFA-plasmid": "#cc79a7",
    "NICD-plasmid": "#e69f00", 
    "shRen-plasmid": "#f0e442", 
    "shtrp53-plasmid": "#d55e00",
    "shPTEN-plasmid" : "#009e73", 
    "shMLL3-plasmid": "#0072b2", 
    "mtCtnnb1": "#000000",  
}

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
            if len(keys) > 1 and len(plasmids) > 1:
                ax = subfigs[si, i].subplots(1, 1, sharey=True)
            elif len(keys) > 1:
                ax = subfigs[si].subplots(1, 1, sharey=True)
            elif len(plasmids) > 1:
                ax = subfigs[i].subplots(1, 1, sharey=True)
            else:
                ax = subfigs[0].subplots(1, 1, sharey=True)

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
                
def plot_spatial_barcodes_average(adata_dict, plasmids_df, plasmids, keys=None):
    """
    Plot the average expression of spatial barcodes in array view for specified samples and plasmids.

    Parameters
    ----------
    adata_dict : dict of AnnData
        A dictionary of AnnData objects indexed by sample identifiers, containing spatial gene expression data.
    plasmids_df : DataFrame
        A DataFrame containing the plasmid-to-gene mapping, with columns 'List' and 'Name'.
    plasmids : list of str
        A list of plasmid names for which the average barcode expression is plotted.
    keys : list of str, optional
        A list of sample keys to include in the plot. If None, all keys from `adata_dict` are used.

    Returns
    -------
    None
        Generates a multi-panel plot, each panel showing the average barcode expression for a given plasmid and sample.
    """
    
    if keys is None:
        keys = list(samples_dict.keys())
    
    fig = plt.figure(layout='constrained', figsize=((3 * len(plasmids), 3 * len(keys))))
    subfigs = fig.subfigures(len(keys), len(plasmids), wspace=0.0)
        
    for si, s in enumerate(keys):
        
        ad = adata_dict[s].copy()
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)

        for i, p in enumerate(plasmids):
            ad.obs[f'expression_{p}'] = compute_quantiles(np.array(ad[:,plasmids_df[plasmids_df['List']== p].Name].X.mean(1))[:,0])

            if len(keys) > 1 and len(plasmids) > 1:
                ax = subfigs[si, i].subplots(1, 1, sharey=True)
            elif len(keys) > 1:
                ax = subfigs[si].subplots(1, 1, sharey=True)
            elif len(plasmids) > 1:
                ax = subfigs[i].subplots(1, 1, sharey=True)
            else:
                ax = subfigs[0].subplots(1, 1, sharey=True)

            
            sc.pl.spatial(ad, img_key='lowres', color=f'expression_{p}',
              cmap='coolwarm',vmin=0, vmax=1, ax=ax, show=False, colorbar_loc=None)
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

def plot_average_barcode_expression(adata_dict, plasmids_df, plasmids, keys=None, ax=None):
    
    reporter_genes_list = []
    reporter_raw_counts = []
    plasmids_df_ordered = plasmids_df.set_index('Name', drop=False).loc[reporters_order]
    for pl in plasmids:
        # sorted_names = sorted(plasmids_df_ordered[plasmids_df_ordered['List'] == pl]['Name'].values,
        #                       key=lambda x: next((i for i, prefix in enumerate(order) if x.startswith(prefix)), len(order)))
        sorted_names = plasmids_df_ordered[plasmids_df_ordered['List'] == pl]['Name'].values
        reporter_genes_list.append(sorted_names)

        norm_adatas = []
        for k in keys:
            adata = adata_dict[k].copy()
            sc.pp.normalize_total(adata, 1e4)
            sc.pp.log1p(adata)
            norm_adatas.append(adata[:,sorted_names].X.todense())

        raw_expression_data = np.concatenate(norm_adatas, axis=0).mean(0)

        reporter_raw_counts.append(np.array(raw_expression_data))

    x_axis = np.ones(plasmids_df.shape[0])
    bc = 0
    x_gap=3
    for i in range(7):
        n_elements = reporter_raw_counts[i].shape[-1]
        x_axis[bc+n_elements] = x_gap
        bc+=n_elements

    x_axis = np.cumsum(x_axis)

    ccolors = [[colours[pl]]*reporter_raw_counts[i].shape[-1] for i, pl in enumerate(plasmids)]

    ys = [x for y in reporter_raw_counts for x in y[0]]
    cs = [x for y in ccolors for x in y]
    xl = [x for y in reporter_genes_list for x in y]
    xs = x_axis

    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(10,3))
        
    ax.plot(np.array([[xe, xe] for xe in xs]).T, np.array([[0, ye] for ye in ys]).T, c='k', alpha=0.7, lw=1, linestyle='--',zorder=0);
    ax.scatter(xs, ys, c = cs)
    ax.set_xticks(xs);
    ax.set_xticklabels(xl, rotation=90)
    ax.set_ylim(0,np.max(ys) + 0.05 * np.max(ys))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Average Expression')
    ax.set_title(f'Average barcode expression aggregated per {", ".join(keys)}')
