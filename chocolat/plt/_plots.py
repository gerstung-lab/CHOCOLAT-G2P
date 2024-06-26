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
    "intersection": "white"
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
    """
    Plots the average barcode expression for a list of plasmids across multiple samples.

    Parameters
    ----------
    adata_dict : dict
        Dictionary of AnnData objects, indexed by sample keys.
    plasmids_df : pandas.DataFrame
        DataFrame containing information about the plasmids.
    plasmids : list of str
        List of plasmid names to be plotted.
    keys : list of str, optional
        List of sample keys to be included in the plot. If None, all keys in `adata_dict` are used.
    ax : matplotlib.axes.Axes, optional
        Axis object to draw the plot onto, otherwise uses the current axis.

    Returns
    -------
    None
    """

    reporter_genes_list = []
    reporter_raw_counts = []
    plasmids_df_ordered = plasmids_df.set_index('Name', drop=False).loc[reporters_order]
    for pl in plasmids:
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
    
def plot_glmm_weights(mean_marginal, variance_marginal, plasmids, feature_list, dim=None):
    """
    Plots the weights of the Generalised Linear Mixed Models (GLMM) for each feature and plasmid.
    Whiskers are 3σ 

    Parameters
    ----------
    mean_marginal : numpy.ndarray
        Mean of the marginal posterior distributions of the model parameters.
    variance_marginal : numpy.ndarray
        Variance of the marginal posterior distributions of the model parameters.
    plasmids : list of str
        List of plasmid names.
    feature_list : list of str
        List of feature names.
    dim : tuple of int, optional
        Dimensions of the subplot grid. Default is (1, len(feature_list)).

    Returns
    -------
    None
    """

    if dim is None:
        dim = (1, len(feature_list))
    plt.figure(figsize=(4*dim[1],4*dim[0]), dpi=100)

    for c, f in enumerate(feature_list):
        plt.subplot(dim[0],dim[1],c+1)
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)
        offset = 0  
        for i in range(len(plasmids) + 1):
            if c == 0:
                plt.text(-3, offset-0.2, (plasmids + ['intersection'])[i].split('-')[0])

            mean = mean_marginal[c,i]
            std = variance_marginal[c,i]**0.5
            colour = colours[(plasmids + ['intersection'])[i]]
            plt.scatter(mean, offset, c=colour, linewidths=1.5, edgecolors='k', s=80, zorder=5, )
            plt.plot([mean - 3*std, mean + 3*std], [offset, offset], lw=1.5, color='k', zorder=2)
            offset -= 2

        plt.axvline(0, c='k', linestyle='--', zorder=-10)
        plt.xlim(-3, 3)
        plt.title(f'{f} expr.')
        
def plot_seqential_pred_correlations(predicted_composition_df, correspondent_annotations, corresponding_samples, corresponding_column):
    """
    Plots the sequential correlations between predicted compositions and annotations for multiple samples.

    Parameters
    ----------
    predicted_composition_df : pandas.DataFrame
        DataFrame containing the predicted compositions.
    correspondent_annotations : pandas.DataFrame
        DataFrame containing the correspondence annotations.
    corresponding_samples : list of tuple of str
        List of tuples with sample names for the primary and replica.
    corresponding_column : list of str
        List of column names corresponding to the annotations.

    Returns
    -------
    None
    """
    xs = []
    ys = []

    for i in range(np.array(corresponding_samples).shape[0]):
        k1 = corresponding_samples[i][0]
        k2 = corresponding_samples[i][1]
        cc = corresponding_column[i]

        correspondence_df_sub = correspondent_annotations.set_index(['sample_name', 'new_annotation']).loc[k1][cc].dropna()

        if i == 0:
            k1 = corresponding_samples[0][1]
            k2 = corresponding_samples[i][0]

            ordered_inedex2 = correspondence_df_sub.index
            ordered_inedex2 = np.array(list(ordered_inedex2))[np.isin(ordered_inedex2, predicted_composition_df.loc[k2].index)]
            ordered_inedex1 = correspondence_df_sub.loc[ordered_inedex2].values.astype(int).astype(str)

        else:
            ordered_inedex1 = correspondence_df_sub.index
            ordered_inedex1 = np.array(list(ordered_inedex1))[np.isin(ordered_inedex1, predicted_composition_df.loc[k1].index)]
            ordered_inedex2 = correspondence_df_sub.loc[ordered_inedex1].values.astype(int).astype(str)

        xs.append(predicted_composition_df.loc[k1].loc[ordered_inedex1])
        ys.append(predicted_composition_df.loc[k2].loc[ordered_inedex2])

        print(k1, k2, len(ordered_inedex1))

    xs = pd.concat(xs)

    ys = pd.concat(ys)


    fig, axs = plt.subplots(1,8, figsize=(24,3))
    for i in range(8):
        x, y = xs[plasmids_ordered_list[i]].values, ys[plasmids_ordered_list[i]].values

        corr_coefficient, _ = pearsonr(x, y)
        slope, _, _, _ = np.linalg.lstsq(x[:, np.newaxis], y, rcond=None)
        # best_fit_line = slope * x

        axs[i].plot([0,1],[0,float(slope)], color='k', linestyle='--')
        axs[i].scatter(x,y, s=10, alpha=0.5, c='tomato')
        # axs[i].plot([xmin,xmax],[y,y],c='tomato', alpha=0.2, lw=0.2)
        # axs[i].plot([x,x],[ymin,ymax],c='tomato', alpha=0.2, lw=0.2)

        axs[i].text(0.1, 0.9, f'r = {np.round(corr_coefficient,2)}')
        axs[i].set_xlim(0,1)
        axs[i].set_ylim(0,1)
        axs[i].spines[['top', 'right']].set_visible(False)
        axs[i].set_title(plasmids_ordered_list[i])
        axs[i].set_ylabel('Replica')
        axs[i].set_xlabel('Primary')
        axs[i].set_xticks([0,0.5,1])
        axs[i].set_xticklabels([0,0.5,1])
        axs[i].set_yticks([0,0.5,1])
        axs[i].set_yticklabels([0,0.5,1])

    plt.tight_layout()
    
def plot_genotype_enrichment_OR(df_combined_per_nodule, ph, plasmid_list, rest_group=None, ax=None):
    """
    Plots the odds ratios of genotype enrichment for specified plasmids.

    Parameters
    ----------
    df_combined_per_nodule : pandas.DataFrame
        DataFrame containing combined data for each nodule.
    ph : str
        Phenotype column to be used for selecting nodules.
    plasmid_list : list of str
        List of plasmids to be included in the analysis.
    rest_group : str, optional
        Column name for the rest group of nodules. If None, the rest group is defined as all nodules
        not in the selected group. Default is None.
    ax : matplotlib.axes.Axes, optional
        Axis object to draw the plot onto, otherwise creates a new axis.

    Returns
    -------
    None
    """

    selected_nodules = df_combined_per_nodule[df_combined_per_nodule[ph] == 1].index
    if rest_group is None:
        rest_nodules = df_combined_per_nodule[~(df_combined_per_nodule[ph] == 1)].index
    else:
        rest_nodules = df_combined_per_nodule[df_combined_per_nodule[rest_group] == 1].index

    selected_plasmid_dat = df_combined_per_nodule.loc[selected_nodules, plasmid_list]
    rest_plasmid_dat = df_combined_per_nodule.loc[rest_nodules, plasmid_list]

    sim_genotypes_selected = np.array([np.random.binomial(1,np.repeat(selected_plasmid_dat.values[None,:,i],
                                                                      20000, axis=0)) for i in range(8)])
    sim_genotypes_full = np.array([np.random.binomial(1,np.repeat(rest_plasmid_dat.values[None,:,i],
                                                                      20000, axis=0)) for i in range(8)])
    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(3, 8))
    epsilon = 1e-10
    na, ne = np.sum(sim_genotypes_full[:,:,:] == 0, axis=2) + epsilon, np.sum(sim_genotypes_full[:,:,:] == 1, axis=2) + epsilon
    pa, pe = np.sum(sim_genotypes_selected[:,:,:] == 0, axis=2) + epsilon, np.sum(sim_genotypes_selected[:,:,:] == 1, axis=2) + epsilon

    log_or = (pe * na) / (pa * ne)
    # log_or = log_onr[:,~np.any(np.isinf(log_or), 0)]
    for i, c in enumerate([colours[x] for x in plasmid_list]):
        ax.scatter(i, np.percentile(log_or[:], 50, 1)[i], edgecolor='k', s=100, color=c, lw=2, zorder=5);
#         print((np.log10(log_or[i]) > 0).mean())
    ax.plot([np.arange(8), np.arange(8)], [np.percentile(log_or[:], 5, 1), np.percentile(log_or[:], 95, 1)], c='k');

    ax.axhline(1, c='k', linestyle='--')
    ax.spines[['right', 'top', 'bottom']].set_visible(False)
    ax.set_xlabel('Plasmid')
    ax.set_ylim(2**(-4),2**(4))
    ax.set_yscale('log', base=2)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_ylabel('Odds ratios')
    ax.set_title(ph + ' ' + str(df_combined_per_nodule[ph].sum()))
    # ax.invert_xaxis()

def plot_heatmap_cutouts(data_matrix, data_names, gene_name, l, r, ax=None):
    """
    Plots a heatmap cutout around a specified gene. Corresponds to the heatmap of GLM weights.

    Parameters
    ----------
    data_matrix : numpy.ndarray
        Matrix containing the data to be plotted.
    data_names : list of str
        List of gene names corresponding to the columns of the data matrix.
    gene_name : str
        Name of the gene around which the heatmap cutout is plotted.
    l : int
        Number of genes to include to the left of the specified gene.
    r : int
        Number of genes to include to the right of the specified gene.
    ax : matplotlib.axes.Axes, optional
        Axis object to draw the plot onto, otherwise creates a new axis.

    Returns
    -------
    None
    """

    fp = data_names.index(gene_name)
    dat_matrix_sub = data_matrix[:,fp-l:fp+r]
    gene_name_list_ordered_sub = data_names[fp-l:fp+r]

    if ax is None:
        fig, ax = plt.subplots(figsize=((r+l)*0.25,2))
    ax.imshow(dat_matrix_sub, aspect='auto', cmap='RdBu_r', vmax=1, vmin=-1, interpolation='none')
    ax.set_xticks(np.arange(len(gene_name_list_ordered_sub))) 
    ax.set_xticklabels(gene_name_list_ordered_sub, rotation=90)
    ax.set_aspect('equal')


