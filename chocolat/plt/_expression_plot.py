import numpy as np
import matplotlib.pyplot as plt
from alpha_shapes import Alpha_Shaper, plot_alpha_shape
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import warnings
import math

def get_rgb_function(cmap, min_value, max_value):
    r"""Generate a function to map continous values to RGB values using colormap between min_value & max_value."""

    if min_value > max_value:
        raise ValueError("Max_value should be greater or than min_value.")

    if min_value == max_value:
        warnings.warn(
            "Max_color is equal to min_color. It might be because of the data or bad parameter choice. "
            "If you are using plot_contours function try increasing max_color_quantile parameter and"
            "removing cell types with all zero values."
        )

        def func_equal(x):
            factor = 0 if max_value == 0 else 0.5
            return cmap(np.ones_like(x) * factor)

        return func_equal

    def func(x):
        return cmap((np.clip(x, min_value, max_value) - min_value) / (max_value - min_value))

    return func


def rgb_to_ryb(rgb):
    """
    Converts colours from RGB colorspace to RYB

    Parameters
    ----------

    rgb
        numpy array Nx3

    Returns
    -------
    Numpy array Nx3
    """
    rgb = np.array(rgb)
    if len(rgb.shape) == 1:
        rgb = rgb[np.newaxis, :]

    white = rgb.min(axis=1)
    black = (1 - rgb).min(axis=1)
    rgb = rgb - white[:, np.newaxis]

    yellow = rgb[:, :2].min(axis=1)
    ryb = np.zeros_like(rgb)
    ryb[:, 0] = rgb[:, 0] - yellow
    ryb[:, 1] = (yellow + rgb[:, 1]) / 2
    ryb[:, 2] = (rgb[:, 2] + rgb[:, 1] - yellow) / 2

    mask = ~(ryb == 0).all(axis=1)
    if mask.any():
        norm = ryb[mask].max(axis=1) / rgb[mask].max(axis=1)
        ryb[mask] = ryb[mask] / norm[:, np.newaxis]

    return ryb + black[:, np.newaxis]


def ryb_to_rgb(ryb):
    """
    Converts colours from RYB colorspace to RGB

    Parameters
    ----------

    ryb
        numpy array Nx3

    Returns
    -------
    Numpy array Nx3
    """
    ryb = np.array(ryb)
    if len(ryb.shape) == 1:
        ryb = ryb[np.newaxis, :]

    black = ryb.min(axis=1)
    white = (1 - ryb).min(axis=1)
    ryb = ryb - black[:, np.newaxis]

    green = ryb[:, 1:].min(axis=1)
    rgb = np.zeros_like(ryb)
    rgb[:, 0] = ryb[:, 0] + ryb[:, 1] - green
    rgb[:, 1] = green + ryb[:, 1]
    rgb[:, 2] = (ryb[:, 2] - green) * 2

    mask = ~(ryb == 0).all(axis=1)
    if mask.any():
        norm = rgb[mask].max(axis=1) / ryb[mask].max(axis=1)
        rgb[mask] = rgb[mask] / norm[:, np.newaxis]

    return rgb + white[:, np.newaxis]

def transform_row(row, pmin, pmax):
    q25, q99 = np.quantile(row, [pmin/100, pmax/100])
    return np.where(
        row > q99, 1,
        np.where(row < q25, 0, (row - q25) / (q99 - q25))
    )

def create_colormap(R, G, B, white_spacing=1):
    spacing = int(white_spacing * 2.55)

    N = 255
    M = 3

    alphas = np.concatenate([[0] * spacing * M, np.linspace(0, 1.0, (N - spacing) * M)])

    vals = np.ones((N * M, 4))
    for i, color in enumerate([R, G, B]):
        vals[:, i] = color / 255
    vals[:, 3] = alphas

    return ListedColormap(vals)

def hexagon_points(center, distance):
    points = [center]
    for i in range(6):
        angle_rad = math.pi / 3 * i
        x = center[0] + distance * math.cos(angle_rad)
        y = center[1] + distance * math.sin(angle_rad)
        points.append((x, y))
    
    return points

def filter_near_duplicates(points, threshold=10):
    filtered_points = []
    for point in points:
        is_duplicate = False
        for other_point in filtered_points:
            distance = np.sqrt((other_point[0] - point[0])**2 + (other_point[1] - point[1])**2)
            if distance < threshold:
                is_duplicate = True
                break
        if not is_duplicate:
            filtered_points.append(point)
    return np.array(filtered_points)

def kernel(w):
    return w**10

def plot_spatial_cell_signature(adata, phen_genes, outline=True, pmin = 25, pmax = 99, s=6, ax=None, cmaps=None, rasterised=False):
    
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", category=UserWarning, message="Some cells have zero counts")

    
    # Function is largely borrowed from cell2location method
    
    # Create linearly scaled colormaps
    if cmaps is None:
        white_spacing = pmin
        YellowCM = create_colormap(240, 228, 66) 
        RedCM = create_colormap(213, 94, 0) 
        BlueCM = create_colormap(86, 180, 233) 
        GreenCM = create_colormap(0, 158, 115) 
        GreyCM = create_colormap(200, 200, 200) 
        WhiteCM = create_colormap(50, 50, 50)  
        PurpleCM = create_colormap(90, 20, 165)  

        cmaps = [YellowCM, RedCM, BlueCM, GreenCM, PurpleCM, GreyCM, WhiteCM]

    alpha_scaling=1

    if outline:
        outline_list = []
        n_nodules = np.sum([x[0].upper() != 'N' for x in adata.obsm['nodules_inside'].columns])
        for i in range(n_nodules):
            try:
                nodule_points = adata.obsm['spatial'][adata.obsm['nodules_inside'].values[:,i] == 1]
                if len(nodule_points) != 0:

                    #double
                    new_coord_hex = []
                    for j in range(nodule_points.shape[0]):
                        new_coord_hex.append(np.array(hexagon_points(nodule_points[j], 198/2+16)))
                    new_coord_hex = filter_near_duplicates(np.array([x for s in new_coord_hex for x in s]))
                    shaper = Alpha_Shaper(new_coord_hex)
                    try:
                        alpha = 12
                        alpha_shape = shaper.get_shape(alpha=alpha)
                        try:
                            outline_list.append(alpha_shape.boundary.coords.xy)
                        except NotImplementedError:
                            alpha = 5
                            alpha_shape = shaper.get_shape(alpha=alpha)
                            outline_list.append(alpha_shape.boundary.coords.xy)
                    except AttributeError:
                        alpha = 2
                        alpha_shape = shaper.get_shape(alpha=alpha)
                        outline_list.append(alpha_shape.boundary.coords.xy)

                    # outline_list.append(outwords_points[hull].T)
                    # plt.plot(*outwords_points[hull].T, c='k')
                    # plt.scatter(*nodule_points.T, c='grey')
                else:
                    # print('empty ' + str(i))
                    pass
            except IndexError:
                print(i)
                
    for k, v in phen_genes.items():
        adata.obs[k] = adata[:,v].X.mean(1)
        vmin = np.percentile(adata.obs[k], pmin)
        vmax = np.percentile(adata.obs[k], pmax)
        adata.obs[k + '_normalised'] = transform_row(adata.obs[k], pmin, pmax)
        
        
    max_color_quantile = pmax / 100

    coords = adata.obsm['spatial']
    counts = adata.obs[[k + '_normalised' for k in phen_genes.keys()]].values

    c_ord = list(np.arange(0, counts.shape[1]))

    colors = np.zeros((*counts.shape, 4))
    weights = np.zeros(counts.shape)

    for c in c_ord:
        min_color_intensity = counts[:, c].min()
        max_color_intensity = np.min([np.quantile(counts[:, c], max_color_quantile), np.inf])

        rgb_function = get_rgb_function(cmap=cmaps[c], min_value=min_color_intensity, max_value=max_color_intensity)

        color = rgb_function(counts[:, c])
        color[:, 3] = color[:, 3] * alpha_scaling

        norm = mpl.colors.Normalize(vmin=min_color_intensity, vmax=max_color_intensity)

        colors[:, c] = color
        weights[:, c] = np.clip(counts[:, c] / (max_color_intensity + 1e-10), 0, 1)
        weights[:, c][counts[:, c] < min_color_intensity] = 0

    colors_ryb = np.zeros((*weights.shape, 3))

    for i in range(colors.shape[0]):
        colors_ryb[i] = rgb_to_ryb(colors[i, :, :3])

    kernel_weights = kernel(weights[:, :, np.newaxis])
    weighted_colors_ryb = (colors_ryb * kernel_weights).sum(axis=1) / kernel_weights.sum(axis=1)

    weighted_colors = np.zeros((weights.shape[0], 4))

    weighted_colors[:, :3] = ryb_to_rgb(weighted_colors_ryb)

    weighted_colors[:, 3] = colors[:, :, 3].max(axis=1)

    if ax is None:
        fig, ax = plt.sbuplots(figsize=(8,8))

    ax.scatter(x=coords[:, 0], y=coords[:, 1], c=weighted_colors, s=s**2, lw=0, zorder=-1)

    if outline:
        for i in range(len(outline_list)):
            ax.fill(*outline_list[i], fill=False, edgecolor='white', lw=1, zorder=3)
    ax.set_aspect('equal')
    ax.set_facecolor('black')
    ax.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if rasterised:
        ax.set_rasterization_zorder(0)
    ax.invert_yaxis()

