import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy.spatial.distance import pdist
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches     import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from scipy import stats
import copy
import networkx.algorithms.community as nx_comm

esfuerzo=pd.Series([7,9,7,9,8,10,12,9], index=[2015,2016,2017,2018,2019,2020,2021,2022])
Years=[2015,2016,2017,2018,2019,2020,2021,2022]
Years_str=[str(i) for i in Years]
Sites=['aznalcazar','bonares','convento de la luz','el pozo','esparragal','la cuña','la rocina','pinares de hinojos','pino del cuervo','urbanizaciones','villamanrique este','villamanrique sur']

def read_database():
    filename="./DATA/Ilist_byyear.csv" 
    
    Ilist=pd.read_csv(filename,index_col=0)
    Ilist.set_index(["Site_ID", "Year"],inplace=True)
    #Ilist["values"]=1

    return Ilist

#helpers for plotting weighted matrices
def _hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
 
def _weight_to_alpha(w, wmax, min_alpha=0.15, max_alpha=1.0):
    """Escala pesos crudos a opacidad (mas peso -> mas intenso)."""
    w = np.asarray(w, dtype=float)
    if wmax is None or wmax <= 0:
        return np.full_like(w, max_alpha)
    frac = np.clip(w / wmax, 0, 1)
    return min_alpha + frac * (max_alpha - min_alpha)

def _joint_order_binary(web_t1, web_t2):
    """
    UNICA fuente de verdad para el orden conjunto de especies (plantas y polinizadores),
    usada tanto por plot_betadiversity_Wmatrix_to_ax como por
    plot_networks_separately_matched_order -- asi las dos figuras SIEMPRE coinciden,
    porque comparten la misma llamada, no una copia duplicada del calculo.

    Orden: comunes primero (por conectividad dentro del bloque comun), luego exclusivas
    de t1, luego exclusivas de t2. Las especies exclusivas se ordenan por su
    conectividad TOTAL dentro de su propia red (incluyendo conexiones con especies
    comunes), no solo su conexion con otras especies tambien exclusivas.
    """
    web_t1_b = (web_t1 != 0).astype(int)
    web_t2_b = (web_t2 != 0).astype(int)

    pollis1 = set(web_t1.columns)
    plants1 = set(web_t1.index)
    pollis2 = set(web_t2.columns)
    plants2 = set(web_t2.index)

    common_pollis = pollis1.intersection(pollis2)
    common_plants = plants1.intersection(plants2)
    pollis_y1 = pollis1 - pollis2
    plants_y1 = plants1 - plants2
    pollis_y2 = pollis2 - pollis1
    plants_y2 = plants2 - plants1

    web_bin_full = web_t1_b.add(web_t2_b, fill_value=0)
    common_plants_ordered = list(
        web_bin_full.loc[sorted(common_plants), sorted(common_pollis)]
        .sum(axis=1).sort_values(ascending=False, kind="stable").index
    )
    common_pollis_ordered = list(
        web_bin_full.loc[sorted(common_plants), sorted(common_pollis)]
        .sum(axis=0).sort_values(ascending=False, kind="stable").index
    )

    # conectividad TOTAL dentro de su propia red (todas las columnas/filas,
    # no solo la submatriz periferia-periferia)
    plants_y1_ordered = list(
        web_t1_b.loc[sorted(plants_y1), :]
        .sum(axis=1).sort_values(ascending=False, kind="stable").index
    )
    pollis_y1_ordered = list(
        web_t1_b.loc[:, sorted(pollis_y1)]
        .sum(axis=0).sort_values(ascending=False, kind="stable").index
    )

    plants_y2_ordered = list(
        web_t2_b.loc[sorted(plants_y2), :]
        .sum(axis=1).sort_values(ascending=False, kind="stable").index
    )
    pollis_y2_ordered = list(
        web_t2_b.loc[:, sorted(pollis_y2)]
        .sum(axis=0).sort_values(ascending=False, kind="stable").index
    )

    ordered_plants = common_plants_ordered + plants_y1_ordered + plants_y2_ordered
    ordered_pollis = common_pollis_ordered + pollis_y1_ordered + pollis_y2_ordered

    return ordered_plants, ordered_pollis

def plot_network_matrix_to_ax(web, ax=None, order=None, cmap="Blues",
                               fmt="{:.0f}", fontsize=8, show_labels=True,
                               show_colorbar=False, vmax=None):
    """Plotea una unica red como matriz de calor, con el valor escrito en cada celda."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    web = web.fillna(0)

    if order is not None:
        ordered_rows, ordered_cols = order
    else:
        ordered_rows = list(web.sum(axis=1).sort_values(ascending=False).index)
        ordered_cols = list(web.sum(axis=0).sort_values(ascending=False).index)

    mat = web.loc[ordered_rows, ordered_cols]
    values = mat.values.astype(float)

    vmax_eff = vmax if vmax is not None else (np.nanmax(values) if values.size else 1)
    norm = Normalize(vmin=0, vmax=vmax_eff)
    cmap_obj = plt.get_cmap(cmap)

    display_values = np.where(values > 0, values, np.nan)
    cax = ax.imshow(display_values, cmap=cmap_obj, norm=norm)

    n_rows, n_cols = values.shape
    for i in range(n_rows):
        for j in range(n_cols):
            v = values[i, j]
            if v <= 0:
                continue
            rgba = cmap_obj(norm(v))
            luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            text_color = "white" if luminance < 0.5 else "black"
            ax.text(j, i, fmt.format(v), ha="center", va="center",
                    fontsize=fontsize, color=text_color)

    if show_labels:
        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(ordered_cols, rotation=90, fontsize=fontsize)
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(ordered_rows, fontsize=fontsize)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    if show_colorbar:
        plt.colorbar(ScalarMappable(norm=norm, cmap=cmap_obj), ax=ax, fraction=0.046, pad=0.04)

    return ax, (ordered_rows, ordered_cols)
#----------------
def plot_betadiversity_Wmatrix_to_ax(web_t1, web_t2, ax=None, order=None, overlap=None,
                                      weighted=True, min_alpha=0.15, vmax=None,
                                      marker_size=35, marker_color="#1a1a1a",
                                      names=False, fontsize=8):
    """
    Plotea la matriz plantas x polinizadores comparando dos redes (t1, t2).

    El ORDEN de especies se calcula con _joint_order_binary (unica fuente de verdad) --
    la MISMA que usa plot_networks_separately_matched_order, para que ambas figuras
    coincidan siempre.
    """
    if not weighted:
        web_t1 = (web_t1 != 0).astype(int)
        web_t2 = (web_t2 != 0).astype(int)

    # ------------------------------------------------------------------
    # 1) CLASIFICACION: usa las redes binarizadas (presencia/ausencia)
    # ------------------------------------------------------------------
    web_t1_b = (web_t1 != 0).astype(int)
    web_t2_b = (web_t2 != 0).astype(int)

    cat = web_t1_b + 2 * web_t2_b  # 0=nada, 1=solo t1, 2=solo t2, 3=permanente
    core_cat = cat

    pollis1 = set(web_t1.columns)
    plants1 = set(web_t1.index)
    pollis2 = set(web_t2.columns)
    plants2 = set(web_t2.index)

    common_pollis = pollis1.intersection(pollis2)
    common_plants = plants1.intersection(plants2)

    pollis_y1 = pollis1 - pollis2
    plants_y1 = plants1 - plants2
    pollis_y2 = pollis2 - pollis1
    plants_y2 = plants2 - plants1

    for p in plants_y1:
        for a in pollis_y1:
            core_cat.loc[p, a] = web_t1_b.loc[p, a]
        for a in common_pollis:
            core_cat.loc[p, a] = web_t1_b.loc[p, a]
    for a in pollis_y1:
        for p in common_plants:
            core_cat.loc[p, a] = web_t1_b.loc[p, a]

    for p in plants_y2:
        for a in pollis_y2:
            core_cat.loc[p, a] = 2 * web_t2_b.loc[p, a]
        for a in common_pollis:
            core_cat.loc[p, a] = 2 * web_t2_b.loc[p, a]
    for a in pollis_y2:
        for p in common_plants:
            core_cat.loc[p, a] = 2 * web_t2_b.loc[p, a]

    for a in pollis_y2:
        for p in plants_y1:
            core_cat.loc[p, a] = -1
    for a in pollis_y1:
        for p in plants_y2:
            core_cat.loc[p, a] = -1

    core_cat.replace(0, np.nan, inplace=True)

    # ------------------------------------------------------------------
    # 2) INTENSIDAD: matriz de pesos crudos, alineada con core_cat.
    # ------------------------------------------------------------------
    weight_web = web_t1.reindex(index=core_cat.index, columns=core_cat.columns).fillna(0) \
        + web_t2.reindex(index=core_cat.index, columns=core_cat.columns).fillna(0)

    diff_web = web_t2.reindex(index=core_cat.index, columns=core_cat.columns).fillna(0) \
        - web_t1.reindex(index=core_cat.index, columns=core_cat.columns).fillna(0)

    # ------------------------------------------------------------------
    # 3) ORDEN de especies -- delegado a _joint_order_binary (UNICA fuente de verdad)
    # ------------------------------------------------------------------
    ordered_plants, ordered_pollis = _joint_order_binary(web_t1, web_t2)

    if order:
        ordered_plants = order[0]
        ordered_pollis = order[1]

    c_grey = '#d7d7d6'
    c_blue = '#0072B2'
    c_yellow = '#E69F00'
    c_green = '#009455'

    if overlap is not None:
        over_ordered = overlap.reindex(index=ordered_plants, columns=ordered_pollis).replace(0, np.nan)
        cax = ax.matshow(over_ordered.loc[ordered_plants, ordered_pollis], cmap="brg", alpha=0.4)
    else:
        cat_ordered = core_cat.loc[ordered_plants, ordered_pollis]
        weight_ordered = weight_web.loc[ordered_plants, ordered_pollis]
        diff_ordered = diff_web.loc[ordered_plants, ordered_pollis]

        cat_vals = cat_ordered.values
        w_vals = weight_ordered.values
        diff_vals = diff_ordered.values

        n_rows, n_cols = cat_vals.shape
        rgba = np.zeros((n_rows, n_cols, 4))

        mask_green = cat_vals == 3
        mask_blue = cat_vals == 1
        mask_yellow = cat_vals == 2

        green_mag = np.abs(diff_vals[mask_green]) if mask_green.any() else np.array([])
        blue_mag = w_vals[mask_blue] if mask_blue.any() else np.array([])
        yellow_mag = w_vals[mask_yellow] if mask_yellow.any() else np.array([])

        if not weighted:
            alpha_fn = lambda mag: np.ones_like(mag, dtype=float)
        elif vmax is not None:
            shared_vmax = vmax
            alpha_fn = lambda mag: _weight_to_alpha(mag, shared_vmax, min_alpha=min_alpha)
        else:
            all_mags = np.concatenate([green_mag, blue_mag, yellow_mag])
            shared_vmax = np.nanmax(all_mags) if all_mags.size > 0 and np.nanmax(all_mags) > 0 else 1
            alpha_fn = lambda mag: _weight_to_alpha(mag, shared_vmax, min_alpha=min_alpha)

        mask = cat_vals == -1
        rgba[mask, 0:3] = _hex_to_rgb(c_grey)
        rgba[mask, 3] = 1.0

        if mask_green.any():
            rgba[mask_green, 0:3] = _hex_to_rgb(c_green)
            rgba[mask_green, 3] = alpha_fn(green_mag)

        if mask_blue.any():
            rgba[mask_blue, 0:3] = _hex_to_rgb(c_blue)
            rgba[mask_blue, 3] = alpha_fn(blue_mag)

        if mask_yellow.any():
            rgba[mask_yellow, 0:3] = _hex_to_rgb(c_yellow)
            rgba[mask_yellow, 3] = alpha_fn(yellow_mag)

        cax = ax.imshow(rgba)

        marker_edgecolor = "white"
        rows_idx, cols_idx = np.where(mask_green)
        changes = diff_vals[mask_green]
        nonzero = changes != 0
        if nonzero.any():
            up = changes[nonzero] > 0
            down = ~up
            r_nz = rows_idx[nonzero]
            c_nz = cols_idx[nonzero]
            if up.any():
                ax.scatter(c_nz[up], r_nz[up], marker='^', s=marker_size,
                           c=c_yellow, edgecolors=marker_edgecolor, linewidths=0.8, zorder=5)
            if down.any():
                ax.scatter(c_nz[down], r_nz[down], marker='v', s=marker_size,
                           c=c_blue, edgecolors=marker_edgecolor, linewidths=0.8, zorder=5)

    c0, c1, c2 = c_green, c_blue, c_yellow
    ax.axvline(len(common_pollis) - 0.5, color=c0, linestyle='--')
    ax.axvline(len(common_pollis) - 0.4, color=c1, linestyle='--')
    ax.axvline(len(common_pollis) + len(pollis_y1) - 0.5, color=c1, linestyle='--')
    ax.axvline(len(common_pollis) + len(pollis_y1) - 0.4, color=c2, linestyle='--')
    ax.axvline(len(common_pollis) + len(pollis_y1) + len(pollis_y2) - 0.5, color=c2, linestyle='--')

    ax.axhline(len(common_plants) - 0.5, color=c0, linestyle='--')
    ax.axhline(len(common_plants) - 0.4, color=c1, linestyle='--')
    ax.axhline(len(common_plants) + len(plants_y1) - 0.5, color=c1, linestyle='--')
    ax.axhline(len(common_plants) + len(plants_y1) - 0.4, color=c2, linestyle='--')
    ax.axhline(len(common_plants) + len(plants_y1) + len(plants_y2) - 0.5, color=c2, linestyle='--')

    if names:
        ax.set_xticks(range(len(ordered_pollis)))
        ax.set_xticklabels(ordered_pollis, rotation=90, fontsize=fontsize)
        ax.set_yticks(range(len(ordered_plants)))
        ax.set_yticklabels(ordered_plants, fontsize=fontsize)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    return ax, (ordered_plants, ordered_pollis)

def plot_networks_separately_matched_order(web_t1, web_t2, axs=None,
                                            cmap1="Blues", cmap2="Oranges",
                                            fmt="{:.0f}", fontsize=8,
                                            show_labels=True, show_colorbar=False,
                                            shared_vmax=True):
    """
    Plotea web_t1 y web_t2 en dos ejes separados, usando el orden de _joint_order_binary
    -- la MISMA fuente que usa plot_betadiversity_Wmatrix_to_ax.
    """
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    ordered_plants, ordered_pollis = _joint_order_binary(web_t1, web_t2)

    rows_t1 = [p for p in ordered_plants if p in web_t1.index]
    cols_t1 = [a for a in ordered_pollis if a in web_t1.columns]
    rows_t2 = [p for p in ordered_plants if p in web_t2.index]
    cols_t2 = [a for a in ordered_pollis if a in web_t2.columns]

    vmax = None
    if shared_vmax:
        m1 = web_t1.values.max() if web_t1.size else 0
        m2 = web_t2.values.max() if web_t2.size else 0
        vmax = max(m1, m2)

    plot_network_matrix_to_ax(web_t1, ax=axs[0], order=(rows_t1, cols_t1), cmap=cmap1,
                               fmt=fmt, fontsize=fontsize, show_labels=show_labels,
                               show_colorbar=show_colorbar, vmax=vmax)
    plot_network_matrix_to_ax(web_t2, ax=axs[1], order=(rows_t2, cols_t2), cmap=cmap2,
                               fmt=fmt, fontsize=fontsize, show_labels=show_labels,
                               show_colorbar=show_colorbar, vmax=vmax)

    return axs, (ordered_plants, ordered_pollis)


def corr_func_annotate(x, y, axi=None, method="pearson", color="black", fontsize=13,square=False,pvalue=True, **kws):

    axi = axi or plt.gca()

    nan_policy='propagate'
    try:
        nan_policy=kws["nan_policy"]
    except:
        pass    

    if (method=="pearson"):
        try:
            r, p = stats.pearsonr(x, y)
        except:
            r,p=(np.nan,np.nan)    
    elif(method =="spearman"):
        try:
            r, p = stats.spearmanr(x, y,nan_policy=nan_policy)
        except:
            r,p=(np.nan,np.nan)  
    print(r)
    print(p)
    x = 0.2
    y = 0.2
    try:
        x=kws["xy"][0]
        y=kws["xy"][1]
    except:
        pass

    try:
         color=kws["color"]
    except:
        pass

    if (square==True):
        axi.annotate(r'$r^{2} = {%.2f}(%.2f)$' % (pow(r,2),p), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)

    else:
        if (pvalue==True):
            if (method=="spearman"):
                axi.annotate(r'$\rho = {%.2f}(%.2f)$'% (r,p), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
            else:
                axi.annotate(r'$r = {%.2f}(%.2f)$'% (r,p), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
        else:
            if (method=="spearman"):
                axi.annotate(r'$\rho = {%.2f}$'% (r), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
            else:
                axi.annotate(r'$r = {%.2f}$'% (r), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
            
def betadiversity_decomposed(web_t1,web_t2):
    #The previous version only measured rewiring changes between species present in both years, 
    # because those that were present in only one year appeared as NaNs. 
    # To fix this, we modified the code to obtain both things.
    
    web=(web_t1+2*web_t2)#.dropna(how="all").dropna(how="all", axis=1)
    #This introduces NaNs in those with only one year, so we have to put them back in
    core_web=web.copy() #copiar el core
    
    pollis1= set(web_t1.columns)
    plants1= set(web_t1.index)

    pollis2= set(web_t2.columns)
    plants2= set(web_t2.index)

    common_pollis = pollis1.intersection(pollis2)
    common_plants = plants1.intersection(plants2)
    
    pollis_y1=pollis1-pollis2
    plants_y1=plants1-plants2
    pollis_y2=pollis2-pollis1
    plants_y2=plants2-plants1

    for p in plants_y1: #copy Y1 interactions
        for a in pollis_y1:
            web.loc[p,a]=web_t1.loc[p,a]
        for a in common_pollis:
            web.loc[p,a]=web_t1.loc[p,a]
    for a in pollis_y1:
        for p in common_plants:
            web.loc[p,a]=web_t1.loc[p,a]

    for p in plants_y2:#cppy Y2 interactions
        for a in pollis_y2:
            web.loc[p,a]=2*web_t2.loc[p,a]
        for a in common_pollis:
            web.loc[p,a]=2*web_t2.loc[p,a]
    for a in pollis_y2:
        for p in common_plants:
            web.loc[p,a]=2*web_t2.loc[p,a]
    for a in pollis_y2:
        for p in plants_y1:
            web.loc[p,a]=-1
    for a in pollis_y1:
        for p in plants_y2:
            web.loc[p,a]=-1
    
    web.replace(0, np.nan, inplace=True)
    #Now we have a matrix where a cell is 3 if the interaction appears in both years, 1 if it appears only in the first, and 2 if it appears only in the second.
    #Calculate beta diversity of interactions
    p=[]
    dif=web #web_t2-web_t1, with this we only got changes in permanent core
    ind=0
    for x in [3,1,2]: # a,b,c
        #print(x)
        boolean_df=dif==x
        count=boolean_df.sum().sum()
        p.append(count)
        ind+=1
    #print(p)
    beta_wn=((p[0]+p[1]+p[2])/((2*p[0]+p[1]+p[2])/2))-1

    #calculate betadiversity of interactions due to rewiring
    p_rw=[]
    dif_rw=core_web
    ind=0
    for x in [3, 1,2]:
        #print(x)
        boolean_df=dif_rw==x
        count=boolean_df.sum().sum()
        p_rw.append(count)
        ind+=1
    #print(p)
    beta_os=((p_rw[0]+p_rw[1]+p_rw[2])/((2*p_rw[0]+p_rw[1]+p_rw[2])/2))-1


    return beta_wn,beta_os

def matrices_to_long_df(overlap_matrices):
    """
    Convert a dict of overlap matrices to a long-form DataFrame.
    
    Parameters:
    -----------
    overlap_matrices : dict
        Keys are (site_ID, year) tuples.
        Values are DataFrames with plant species as index, 
        pollinator species as columns, and binary overlap values.
    
    Returns:
    --------
    long_df : pandas.DataFrame
        Columns: ['Site_ID', 'Year', 'Plant_gen_sp', 'Pollinator_gen_sp']
        Each row corresponds to a cell==1 in the matrices.
    """
    records = []
    for (site, year), mat in overlap_matrices.items():
        # Stack the DataFrame to long form, drop zeros
        stacked = mat.stack().reset_index()
        stacked.columns = ['Plant_gen_sp', 'Pollinator_gen_sp', 'overlap']
        
        # Filter for overlap == 1
        stacked = stacked[stacked['overlap'] == 1]
        
        # Add site and year
        stacked['Site_ID'] = site
        stacked['Year']    = year
        
        # Keep only the desired columns
        stacked = stacked[['Site_ID', 'Year', 'Plant_gen_sp', 'Pollinator_gen_sp']]
        
        records.append(stacked)
    
    # Concatenate all site-year DataFrames
    long_df = pd.concat(records, ignore_index=True)
    return long_df

def get_metanetworks(Ilist):
    # Step 1: Create the dictionary to hold matrices
    metanetworks = {}

    # Step 2: Loop over each site
    for site, group in Ilist.reset_index().groupby("Site_ID"):
        print(site)
        # Optionally drop duplicates (interactions might be repeated across years)
        interaction_counts = (
            group[["Plant_gen_sp", "Pollinator_gen_sp"]]
            .drop_duplicates()
            .value_counts()
            .reset_index(name="value")
        )

        # Create the incidence matrix using pivot_table
        matrix = interaction_counts.pivot_table(
            index="Plant_gen_sp",
            columns="Pollinator_gen_sp",
            values="value",
            fill_value=0
        )

        # Store in the dictionary
        metanetworks[site] = matrix

    return metanetworks


    # Pseudo R²
    r2_pseudo = ssr / sst
    return r2_pseudo

def concatenate_values(row):
    return f"{row['Estimate']} ({row['p-value']})"




