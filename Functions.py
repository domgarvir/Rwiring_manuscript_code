import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator

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

#plotting bipartite networks
def get_betadiversity_network_by_year(web_t1,web_t2):
    
    web=(web_t1+2*web_t2)#.dropna(how="all").dropna(how="all", axis=1)
    #his introduces NaNs in those with only one year, so we have to put them back in
    #core_web=web.dropna(how="all").dropna(how="all", axis=1)
    core_web=web
    
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

    B1 = nx.Graph()
    B2 = nx.Graph()
    beta_node_dict={}
    for a in common_pollis:
        B1.add_node(a, bipartite='Pollinator', category="common")
        B2.add_node(a, bipartite='Pollinator', category="common")
        beta_node_dict[a]="common"
    for a in pollis_y1:
        B1.add_node(a, bipartite='Pollinator', category="y1" )
        beta_node_dict[a]="y1" 
    for a in pollis_y2:
        B2.add_node(a, bipartite='Pollinator', category="y2" )
        beta_node_dict[a]="y2" 
    for p in common_plants:
        B1.add_node(p, bipartite='Plant', category="common" )
        B2.add_node(p, bipartite='Plant', category="common" )
        beta_node_dict[p]="common"
    for p in plants_y1:
        B1.add_node(p, bipartite='Plant', category="y1" )
        beta_node_dict[p]="y1" 
    for p in plants_y2:
        B2.add_node(p, bipartite='Plant', category="y2" )
        beta_node_dict[p]="y2"     

    for p in plants_y1:
        for a in pollis_y1:
            core_web.loc[p,a]=web_t1.loc[p,a]
        for a in common_pollis:
            core_web.loc[p,a]=web_t1.loc[p,a]
    for a in pollis_y1:
        for p in common_plants:
            core_web.loc[p,a]=web_t1.loc[p,a]

    for p in plants_y2:
        for a in pollis_y2:
            core_web.loc[p,a]=2*web_t2.loc[p,a]
        for a in common_pollis:
            core_web.loc[p,a]=2*web_t2.loc[p,a]
    for a in pollis_y2:
        for p in common_plants:
            core_web.loc[p,a]=2*web_t2.loc[p,a]
    
    
    core_web.replace(0, np.nan, inplace=True)
    Ilist=core_web.stack().to_frame()
    Ilist.columns=["weigth"]
    Ilist1=Ilist[Ilist!=2].dropna()
    Ilist2=Ilist[Ilist>1].dropna()
    edges1=[]
    for (plant, pollinator), weigth in Ilist1['weigth'].items():
        edges1.append((plant, pollinator, weigth))
    B1.add_weighted_edges_from(edges1)

    edges2=[]
    for (plant, pollinator), weigth in Ilist2['weigth'].items():
        edges2.append((plant, pollinator, weigth))
    B2.add_weighted_edges_from(edges2)
      
    return B1, B2, beta_node_dict

def get_betadiversity_Wnetwork_by_year(web_t1, web_t2):
    """
    Construye dos grafos bipartitos (B1 para t1, B2 para t2) a partir de dos
    matrices de interaccion, clasificando:

      - cada NODO (planta o polinizador) segun si es 'common' (presente en
        ambos anios), 'y1' (solo en t1) o 'y2' (solo en t2) -- igual que antes.
      - cada ARISTA segun si esa interaccion concreta (par planta-polinizador)
        es 'common' (persistente, tambien presente en el otro año) o exclusiva
        de este año ('y1' en B1, 'y2' en B2) -- esto es lo nuevo, integrado
        aqui para no tener que llamar aparte a annotate_edge_categories.

    El peso ('weight') de cada arista es siempre el peso REAL de ESE año
    (no se toca, no se promedia, no se combina) -- solo se añade la
    'category' como informacion adicional para poder colorear los links
    (verde=common, color propio=y1/y2) sin depender del valor del peso.

    Devuelve
    -------
    B1, B2 : nx.Graph
        Grafos bipartitos de cada red. Cada nodo tiene 'bipartite' y
        'category' (common/y1/y2). Cada arista tiene 'weight' (peso real de
        ese año) y 'category' (common/y1 en B1; common/y2 en B2).
    beta_node_dict : dict
        {nombre_especie: categoria} para todas las especies de ambas redes.
    """
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

    # ------------------------------------------------------------------
    # 1) Clasificacion de nodos (igual que siempre, sin cambios)
    # ------------------------------------------------------------------
    B1 = nx.Graph()
    B2 = nx.Graph()
    beta_node_dict = {}

    for a in common_pollis:
        B1.add_node(a, bipartite='Pollinator', category="common")
        B2.add_node(a, bipartite='Pollinator', category="common")
        beta_node_dict[a] = "common"
    for a in pollis_y1:
        B1.add_node(a, bipartite='Pollinator', category="y1")
        beta_node_dict[a] = "y1"
    for a in pollis_y2:
        B2.add_node(a, bipartite='Pollinator', category="y2")
        beta_node_dict[a] = "y2"

    for p in common_plants:
        B1.add_node(p, bipartite='Plant', category="common")
        B2.add_node(p, bipartite='Plant', category="common")
        beta_node_dict[p] = "common"
    for p in plants_y1:
        B1.add_node(p, bipartite='Plant', category="y1")
        beta_node_dict[p] = "y1"
    for p in plants_y2:
        B2.add_node(p, bipartite='Plant', category="y2")
        beta_node_dict[p] = "y2"

    # ------------------------------------------------------------------
    # 2) Edges -- construidos directamente desde web_t1/web_t2 (sin decodificar
    #    core_web=web_t1+2*web_t2, para no reintroducir el bug de siempre).
    # ------------------------------------------------------------------
    edges1 = {(p, a): w for (p, a), w in web_t1.stack().items() if w != 0}
    edges2 = {(p, a): w for (p, a), w in web_t2.stack().items() if w != 0}

    B1.add_weighted_edges_from((p, a, w) for (p, a), w in edges1.items())
    B2.add_weighted_edges_from((p, a, w) for (p, a), w in edges2.items())

    # ------------------------------------------------------------------
    # 3) NUEVO -- category por arista: 'common' si el mismo par (planta,
    #    polinizador) tambien interactua en el otro año, 'y1'/'y2' si no.
    #    Comparar como frozenset porque el orden (p,a) vs (a,p) no importa.
    # ------------------------------------------------------------------
    pairs1 = {frozenset(k) for k in edges1}
    pairs2 = {frozenset(k) for k in edges2}

    for u, v in B1.edges():
        B1.edges[u, v]['category'] = "common" if frozenset((u, v)) in pairs2 else "y1"
    for u, v in B2.edges():
        B2.edges[u, v]['category'] = "common" if frozenset((u, v)) in pairs1 else "y2"

    return B1, B2, beta_node_dict

#helper for plotting bipartite networks
def calculate_positions(nodes, node_widths, buffer):
    positions = {}
    current_x = 0
    for node in nodes:
        width = node_widths[node]
        positions[node] = (current_x + width / 2, 0)  # Position at the center of the width
        current_x += width + buffer  # Move to the next position with buffer space
    return positions, current_x

#plot bipartite networks
def plot_bipartite_network_to_ax(B,draw_labels=False,node_cathegories=None,height = 2,w=2,edgecolor="None",h=0.4,buffer=4,lw=1.2,ax=None,colored=True):
    beta_link_dict={1.0:'y1',2.0:'y2',3.0:'common'}
    if colored:
        node_color_dict={'common':'#009455','y1':'#73d2de','y2':"#ffd23f"}
    else:
        node_color_dict={'common':'#535353','y1':'#949494','y2':"#949494"}      

    plants = [n for n, d in B.nodes(data=True) if d['bipartite'] == 'Plant']
    pollinators = [n for n, d in B.nodes(data=True) if d['bipartite'] == 'Pollinator']

    # Calculate degree of nodes
    degree_dict = dict(B.degree(weight='weight'))
    # Define the custom order of categories
    category_order = [ "y2", "y1","common"]
    sorted_plants = sorted(plants, key=lambda x: (degree_dict[x], category_order.index(B.nodes[x]['category'])), reverse=True)
    sorted_pollinators = sorted(pollinators, key=lambda x: (degree_dict[x], category_order.index(B.nodes[x]['category'])), reverse=True)

    # Normalize the degree for node width scaling
    max_degree = max(degree_dict.values())
    w=w
    node_widths = {node: (degree / max_degree) * w for node, degree in degree_dict.items()}  # Adjust scale as needed

    #sorted edges
    edge_categories = set([d['weight'] for u, v, d in B.edges(data=True)])
    # Sort edge categories so higher categories are plotted last
    sorted_edge_categories = sorted(edge_categories, reverse=False)

    # Create a layout for the nodes with buffer space
    pos = {}
    height = height
    buffer = w/buffer
  
    # Position plants at the bottom, sorted by degree
    plant_positions, max_x_plant = calculate_positions(sorted_plants, node_widths, buffer)
    for node, (x, y) in plant_positions.items():
        pos[node] = (x, 0)

    # Position pollinators at the top, sorted by degree
    pollinator_positions, max_x_pollinator = calculate_positions(sorted_pollinators, node_widths, buffer)
    for node, (x, y) in pollinator_positions.items():
        pos[node] = (x, height)

    # Determine the maximum x value for setting plot limits
    max_x = max(max_x_plant, max_x_pollinator)


    for category in sorted_edge_categories:
        edges = [(u, v, d['weight']) for u, v, d in B.edges(data=True) if d['weight'] == category]
        for edge in edges:
            x1, y1 = pos[edge[0]]
            x2, y2 = pos[edge[1]]
            cosa=edge[2]
            if (node_cathegories):
                category = beta_link_dict[cosa]
                link_color = node_color_dict[category]
                if (cosa>2):
                    lwi=lw*1.8
                else:
                    lwi=lw
                ax.plot([x1, x2], [y1, y2], '-', lw=lwi, color=link_color,zorder=1)
            else:
                ax.plot([x1, x2], [y1, y2], '-', lw=lwi, zorder=1)   


    edgecolor=edgecolor
    h=h
    # Draw nodes as rectangles after edges
    for node in sorted_plants:
        x, y = pos[node]
        width = node_widths[node]
        if (node_cathegories):
            category = node_cathegories[node]
            color = node_color_dict[category]
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor=color, zorder=2)
        else:
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor='grey', zorder=2)
        ax.add_patch(rect)

    for node in sorted_pollinators:
        x, y = pos[node]
        width = node_widths[node]
        if (node_cathegories):
            category = node_cathegories[node]
            color = node_color_dict[category]
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor=color, zorder=2)
        else:
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor='black', zorder=2)
        ax.add_patch(rect)

    if (draw_labels):
    # Draw labels for better readability
        for node, (x, y) in pos.items():
            if node in sorted_plants:
                ax.text(x, y - 0.25, node, ha='center', va='top', fontsize=10, rotation=90, zorder=3)  # Below the rectangles for plants
            else:
                ax.text(x, y + 0.25, node, ha='center', va='bottom', fontsize=10, rotation=90, zorder=3)  # Above the rectangles for pollinators

    # Set limits and display the graph
    ax.set_xlim(-1, max_x + buffer)  # Adjust x limits based on cumulative width and buffer
    ax.set_ylim(-0.5, height+0.5)
    ax.set_aspect('equal')
    #plt.title('Bipartite Network of Plants and Pollinators')
    ax.axis('off')
    #plt.tight_layout()
    #filename="../Output/Network.pdf"
    #plt.savefig(filename)
    #plt.show()

    return ax

#plot bipartite weighted networks
def plot_bipartite_Wnetwork_to_ax(B, draw_labels=False, node_cathegories=None,
                                   height=2, w=2, edgecolor="None", h=0.4, buffer=4,
                                   lw=1.2, ax=None, colored=True,
                                   edge_color="#0072B2", edge_color_common="#009455",
                                   min_alpha=0.15, vmax=None, lw_range=None):
    """
    Dibuja UN panel = UN grafo (B3 o B4, uno cada vez -- paneles separados por año).

    Requiere que las aristas de B ya tengan el atributo 'category' puesto por
    annotate_edge_categories(B3, B4) ANTES de llamar a esta funcion:

        B3, B4, beta_dict2 = get_betadiversity_Wnetwork_by_year(web_t1, web_t2)
        B3, B4 = annotate_edge_categories(B3, B4)
        plot_bipartite_Wnetwork_to_ax(B3, ..., edge_color="#0072B2")
        plot_bipartite_Wnetwork_to_ax(B4, ..., edge_color="#E69F00")

    Color de las aristas: 'common' (persistente, tambien en el otro año) ->
    edge_color_common (verde); 'y1'/'y2' (exclusiva de este año) -> edge_color
    (el que pases en la llamada). La intensidad (alpha) y el grosor (si pasas
    lw_range) se calculan a partir del peso REAL de cada arista en ESTE grafo
    -- ya no se usa el peso para decidir el color, solo la categoria.
    """
    if colored:
        node_color_dict = {'common': '#009455', 'y1': '#0072B2', 'y2': '#E69F00'}
    else:
        node_color_dict = {'common': '#535353', 'y1': '#949494', 'y2': '#949494'}

    edge_color_dict = {'common': edge_color_common, 'y1': edge_color, 'y2': edge_color}

    plants = [n for n, d in B.nodes(data=True) if d['bipartite'] == 'Plant']
    pollinators = [n for n, d in B.nodes(data=True) if d['bipartite'] == 'Pollinator']

    degree_dict = dict(B.degree(weight='weight'))
    category_order = ["y2", "y1", "common"]
    sorted_plants = sorted(plants, key=lambda x: (degree_dict[x], category_order.index(B.nodes[x]['category'])), reverse=True)
    sorted_pollinators = sorted(pollinators, key=lambda x: (degree_dict[x], category_order.index(B.nodes[x]['category'])), reverse=True)

    max_degree = max(degree_dict.values())
    node_widths = {node: (degree / max_degree) * w for node, degree in degree_dict.items()}

    pos = {}
    buffer = w / buffer

    plant_positions, max_x_plant = calculate_positions(sorted_plants, node_widths, buffer)
    for node, (x, y) in plant_positions.items():
        pos[node] = (x, 0)

    pollinator_positions, max_x_pollinator = calculate_positions(sorted_pollinators, node_widths, buffer)
    for node, (x, y) in pollinator_positions.items():
        pos[node] = (x, height)

    max_x = max(max_x_plant, max_x_pollinator)

    # ------------------------------------------------------------------
    # ARISTAS -- bucle por categoria (verde primero, luego el color propio
    # encima), intensidad/grosor segun el PESO REAL, no segun la categoria.
    # ------------------------------------------------------------------
    weights = [d['weight'] for u, v, d in B.edges(data=True)]
    vmax_eff = vmax if vmax is not None else (max(weights) if weights else 1)

    sorted_edge_categories = ['common', 'y1', 'y2']  # verde debajo, color propio encima

    for category in sorted_edge_categories:
        edges = [(u, v, d['weight']) for u, v, d in B.edges(data=True) if d['category'] == category]
        if not edges:
            continue
        link_color = edge_color_dict[category]
        # dentro de cada categoria, dibuja antes las de menor peso
        for u, v, wgt in sorted(edges, key=lambda e: e[2]):
            x1, y1_ = pos[u]
            x2, y2_ = pos[v]
            frac = np.clip(wgt / vmax_eff, 0, 1) if vmax_eff > 0 else 1
            alpha = min_alpha + frac * (1 - min_alpha)
            lwi = lw if lw_range is None else lw_range[0] + frac * (lw_range[1] - lw_range[0])
            ax.plot([x1, x2], [y1_, y2_], '-', lw=lwi, color=link_color, alpha=alpha, zorder=1)

    # ------------------------------------------------------------------
    # NODOS: sin cambios
    # ------------------------------------------------------------------
    for node in sorted_plants:
        x, y = pos[node]
        width = node_widths[node]
        if node_cathegories:
            category = node_cathegories[node]
            color = node_color_dict[category]
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor=color, zorder=2)
        else:
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor='grey', zorder=2)
        ax.add_patch(rect)

    for node in sorted_pollinators:
        x, y = pos[node]
        width = node_widths[node]
        if node_cathegories:
            category = node_cathegories[node]
            color = node_color_dict[category]
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor=color, zorder=2)
        else:
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor='black', zorder=2)
        ax.add_patch(rect)

    if draw_labels:
        for node, (x, y) in pos.items():
            if node in sorted_plants:
                ax.text(x, y - 0.25, node, ha='center', va='top', fontsize=10, rotation=90, zorder=3)
            else:
                ax.text(x, y + 0.25, node, ha='center', va='bottom', fontsize=10, rotation=90, zorder=3)

    ax.set_xlim(-1, max_x + buffer)
    ax.set_ylim(-0.5, height + 0.5)
    ax.set_aspect('equal')
    ax.axis('off')

    return ax
#-------------------

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




# PLOTS  ##############################3
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator
label_dict={'CN_CI':'Persistent: constant', 'CN_DI':'Persistent: decreasing strength', 'CN_II': 'Persistent: increasing strength', 'CN_LI':'Lost: rewiring', 
'CN_GI': 'Gained: rewiring', 'L1N_LI':'Lost: one transient sp.', 'LBN_LI': 'Lost: both transient sp.', 
'G1N_GI':'Gained: one transient species', 'GBN_GI':'Gained: both transient sp.'}
label_dict2={'CN_CI':'Persistent:\n constant', 'CN_DI':'Persistent:\n decreasing strength', 'CN_II': 'Persistent:\n increasing strength', 'CN_LI':'Lost:\n rewiring', 
'CN_GI': 'Gained:\n rewiring', 'L1N_LI':'Lost:\n one transient sp.', 'LBN_LI': 'Lost:\n both transient sp.', 
'G1N_GI':'Gained:\n one transient species', 'GBN_GI':'Gained:\n both transient sp.'}
year_label_dict={2016:'2015-16', 2017:'2016-17', 2018:'2017-18', 2019:'2018-19', 2020:'2019-20', 2021:'2020-21', 2022:'2021-22'}

# Plots showing each year ######################################333
def plot_diverging_bars(df_site, value_col='freq_sum', ax=None,
                         loss_metrics=('CN_LI', 'L1N_LI', 'LBN_LI'),
                         gain_metrics=('CN_GI', 'G1N_GI', 'GBN_GI'),
                        loss_colors=('#C6DFEE', '#72B1D4', '#0072B2'),
                        gain_colors=('#F9E9C6', '#F1CA72', '#E69F00'),
                         #
                         year_col='year', metric_col='metric'):
    """
    c_blue = '#0072B2'
    c_yellow = '#E69F00'
    c_green = '#009455'
    Grafico de barras apiladas DIVERGENTE: las metricas de perdida (loss_metrics)
    se apilan hacia ABAJO (valores negativos), las de ganancia (gain_metrics) se
    apilan hacia ARRIBA (valores positivos) -- una barra por año.

    (Sin cambios respecto a la version anterior -- se mantiene para poder
    seguir usandola sola si no quieres el panel superior.)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    df_site = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()

    years = sorted(df_site[year_col].unique())
    x = np.arange(len(years))

    top = np.zeros(len(years))
    for metric, color in zip(gain_metrics, gain_colors):
        vals = np.array([
            df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
        #ax.bar(x, vals, bottom=top, color=color, label=metric, width=0.6)
        ax.bar(x, vals, bottom=top, color=color, label=label_dict[metric], width=0.6)
        top += vals

    bottom = np.zeros(len(years))
    for metric, color in zip(loss_metrics, loss_colors):
        vals = np.array([
            -df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
        #ax.bar(x, vals, bottom=bottom, color=color, label=metric, width=0.6)
        ax.bar(x, vals, bottom=bottom, color=color, label=label_dict[metric], width=0.6)
        bottom += vals



    ax.axhline(0, color='black', lw=1)
    ax.set_xticks(x)
    #ax.set_xticklabels(years)
    ax.set_xticklabels([year_label_dict.get(y, y) for y in years])
    ax.set_ylabel(value_col)
    ax.set_xlabel(year_col)
    #ax.legend(bbox_to_anchor=(1.02, 1), loc='center left', frameon=False)
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)

    return ax


def plot_constant_bars(df_site, value_col='freq_sum', ax=None,
                        constant_metrics=('CN_DI', 'CN_CI', 'CN_II'),
constant_facecolor='#009455',
constant_edgecolors=('#C6DFEE', None, '#F9E9C6'),
                        constant_hatches=('///', None, '///'),
                        hatch_linewidth=2.0,
                        year_col='year', metric_col='metric'):
    """
    Panel SUPERIOR: barras apiladas (NO divergentes, todo hacia arriba) con las
    categorias del bloque de especies permanentes que NO implican perdida/ganancia
    de nodo -- CN_DI (disminuye), CN_CI (constante), CN_II (aumenta).

    Todas usan el MISMO color de fondo gris (constant_facecolor), para que no se
    confundan con los colores solidos rojo/azul del panel divergente. Para poder
    distinguir DI de CI de II sin añadir colores nuevos, se usa un patron de
    rayas (hatch) en el borde:
      - CN_DI: gris de fondo + rayas ROJAS (mismo rojo que las perdidas de abajo)
      - CN_CI: gris liso, SIN rayas (no cambia -- ni perdida ni ganancia)
      - CN_II: gris de fondo + rayas AZULES (mismo azul que las ganancias de abajo)

    Las rayas se dibujan SIN el marco/borde exterior de la barra (linewidth=0
    en la barra), y su grosor se controla aparte via hatch_linewidth -- son dos
    cosas independientes en matplotlib (el borde del rectangulo usa 'linewidth',
    las rayas del hatch usan rcParams['hatch.linewidth']).

    Mismo eje x (años) que plot_diverging_bars, pensado para compartir figura.

    Parameters
    ----------
    constant_facecolor : str
        Color de fondo UNICO para las tres categorias (gris, para no competir
        con el rojo/azul del panel de abajo).
    constant_edgecolors : tuple de (str o None)
        Color de las rayas de cada categoria, en el mismo orden que
        constant_metrics. None = sin rayas (gris liso).
    constant_hatches : tuple de (str o None)
        Patron de rayas de cada categoria (p.ej. '///', 'xxx'), mismo orden
        que constant_metrics. None = sin rayas.
    hatch_linewidth : float
        Grosor de las rayas del hatch (NO del borde de la barra, que se deja
        siempre en 0). Sube este valor para rayas mas gruesas/visibles.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 2.5))

    df_site = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()

    years = sorted(df_site[year_col].unique())
    x = np.arange(len(years))

    bottom = np.zeros(len(years))
    # IMPORTANTE: hatch.linewidth se lee en el momento de RENDERIZAR la figura
    # (plt.savefig/plt.show), no en el momento de crear las barras con ax.bar().
    # Por eso NO usamos un "with rc_context(...)" (se cerraria antes de que el
    # usuario llame a savefig fuera de esta funcion, y el valor se perderia sin
    # tener efecto). En su lugar lo fijamos de forma persistente en rcParams.
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = hatch_linewidth

    for metric, edgecolor, hatch in zip(constant_metrics, constant_edgecolors, constant_hatches):
        vals = np.array([
            df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
        # linewidth=0 siempre -- sin marco/borde exterior; las rayas del hatch
        # se dibujan igualmente (usan hatch.linewidth, no el linewidth del patch)
        ax.bar(x, vals, bottom=bottom, facecolor=constant_facecolor,
               edgecolor=edgecolor if edgecolor else constant_facecolor,
               #hatch=hatch, linewidth=0, label=metric, width=0.6)
               hatch=hatch, linewidth=0, label=label_dict[metric], width=0.6)
        bottom += vals

    ax.set_xticks(x)
    ax.set_ylabel(value_col)
    #ax.legend(bbox_to_anchor=(1.02, 1), loc='center left', frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    order = [1, 2, 0]  # CN_CI, CN_II, CN_DI (indices sobre constant_metrics=('CN_DI','CN_CI','CN_II'))
    ax.legend([handles[i] for i in order], [labels[i] for i in order],
          bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)
    #ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)

    return ax


def plot_diverging_bars_with_constant(df_site, value_col='n_interactions', figsize=(16, 15),
                                       height_ratios=(1, 2), match_scale=True,
                                       ylim_top=None,
                                       loss_metrics=('CN_LI', 'L1N_LI', 'LBN_LI'),
                                       gain_metrics=('CN_GI', 'G1N_GI', 'GBN_GI'),
                                       loss_colors=('#0072B2', '#72B1D4','#C6DFEE' ),
                                       gain_colors=('#E69F00', '#F1CA72','#F9E9C6' ),
                                       constant_metrics=('CN_DI', 'CN_CI', 'CN_II'),
                                       constant_facecolor='#009455',
                                       constant_edgecolors=('#C6DBEF', None, '#FDD0A2'),
                                       constant_hatches=('///', None, '///'),
                                       hatch_linewidth=2.0,
                                       year_col='year', metric_col='metric',
                                       rotate_years=90, year_label_dict=None,
                                       ax_top=None, ax_bottom=None,fontsize=20):
    """
    Figura de UN sitio con DOS paneles, eje x (año) compartido.
 
    NUEVO: acepta ax_top/ax_bottom externos -- si se pasan AMBOS, dibuja
    sobre esos ejes en vez de crear una figura nueva. Se corrige tambien
    el bug de year_label_dict usado sin ser parametro.
 
    Devuelve
    -------
    Si ax_top/ax_bottom son None: fig, (ax_top, ax_bottom)
    Si se pasaron: (ax_top, ax_bottom)
    """
    if year_label_dict is None:
        year_label_dict = {}
 
    created_own_fig = ax_top is None and ax_bottom is None
    if created_own_fig:
        fig, (ax_top, ax_bottom) = plt.subplots(
            2, 1, figsize=figsize, sharex=True,
            gridspec_kw={'height_ratios': height_ratios}
        )
 
    plot_constant_bars(df_site, value_col=value_col, ax=ax_top,
                        constant_metrics=constant_metrics,
                        constant_facecolor=constant_facecolor,
                        constant_edgecolors=constant_edgecolors,
                        constant_hatches=constant_hatches,
                        hatch_linewidth=hatch_linewidth,
                        year_col=year_col, metric_col=metric_col)
 
    plot_diverging_bars(df_site, value_col=value_col, ax=ax_bottom,
                         loss_metrics=loss_metrics, gain_metrics=gain_metrics,
                         loss_colors=loss_colors, gain_colors=gain_colors,
                         year_col=year_col, metric_col=metric_col)
 
    if ylim_top is not None:
        ax_top.set_ylim(*ylim_top)
    elif match_scale:
        ymax_bottom_positive = ax_bottom.get_ylim()[1]
        ax_top.set_ylim(0, ymax_bottom_positive)
 
    ax_top.yaxis.set_major_locator(MultipleLocator(20))
    ax_bottom.yaxis.set_major_locator(MultipleLocator(20))
 
    df_flat = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site
    years = sorted(df_flat[year_col].unique())
    ax_bottom.set_xticklabels([year_label_dict.get(y, y) for y in years], rotation=rotate_years,fontsize=fontsize)
    ax_top.set_xlabel('')
 
    ax_top.set_ylabel('Permanent\n interactions')
    ax_bottom.set_ylabel('Transient\n interactions')
 
    if created_own_fig:
        plt.tight_layout()
        return fig, (ax_top, ax_bottom)
    else:
        return ax_top, ax_bottom
 

def plot_species_turnover(df_site, ax=None, figsize=(8, 5),
                           in_metric='A2', out_metric='A1', permanent_metric='PA',
                           in_color='#E69F00', out_color='#0072B2',
                           total_color='gray', permanent_color='#009455',
                           year_col='year', metric_col='metric', value_col='n_nodes',
                           bar_width=0.6, marker_size=6, rotate_years=90,
                           year_label_dict=None, fontsize=20, show_legend=True):
    """
    Turnover de especies (p.ej. polinizadores) por año, para UN sitio.
 
    Eje x: años. Para cada año:
      - Barra hacia ARRIBA: numero de especies que ENTRAN ese año (in_metric,
        por defecto 'A2' -- exclusivas del año 2 de esa transicion).
      - Barra hacia ABAJO (negativa): numero de especies que SALEN (out_metric,
        por defecto 'A1' -- exclusivas del año 1).
      - Punto GRIS conectado por linea: numero TOTAL de especies ese año
        (permanentes + entran + salen -- PA+A1+A2).
      - Punto de otro color conectado por linea: numero de especies
        PERMANENTES ese año (permanent_metric, por defecto 'PA').
 
    Parameters
    ----------
    df_site : pd.DataFrame
        Ya filtrado a UN Site_ID (p.ej. df.loc["sitio_1"]), con columnas o
        niveles year_col/metric_col/value_col -- salida de rewiring_freq_py
        (la segunda tabla, node_counts) apilada por año.
    in_metric, out_metric, permanent_metric : str
        Nombres de fila en metric_col a usar para entran/salen/permanentes.
        Por defecto para polinizadores: 'A2'/'A1'/'PA'. Para plantas usa
        'P2'/'P1'/'PP'.
    year_label_dict : dict o None
        {año: texto legible} para el eje x.
 
    Devuelve
    -------
    ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if year_label_dict is None:
        year_label_dict = {}
 
    df_site = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()
 
    years = sorted(df_site[year_col].unique())
    x = np.arange(len(years))
 
    def get_vals(metric):
        return np.array([
            df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
 
    in_vals = get_vals(in_metric)
    out_vals = get_vals(out_metric)
    permanent_vals = get_vals(permanent_metric)
    total_vals = permanent_vals + in_vals + out_vals
 
    ax.bar(x, in_vals, color=in_color, width=bar_width, label='Gained', zorder=2)
    ax.bar(x, -out_vals, color=out_color, width=bar_width, label='Lost', zorder=2)
 
    ax.plot(x, total_vals, '-o', color=total_color, markersize=marker_size,
            label='Total', zorder=3)
    ax.plot(x, permanent_vals, '-o', color=permanent_color, markersize=marker_size,
            label='Permanent', zorder=3)
 
    ax.axhline(0, color='black', lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([year_label_dict.get(y, y) for y in years], rotation=rotate_years, fontsize=fontsize)

    ax.set_ylabel(value_col, fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    if show_legend:
        ax.legend(frameon=False, fontsize=fontsize)
 
    return ax

# final plot, combination of the itneractions and species
def plot_species_and_turnover_summary(
        df_site_species, df_site_interactions,
        figsize=(20, 12), height_ratios=(1, 2),
        plant_kwargs=None, pollinator_kwargs=None, diverging_kwargs=None,rotate_years=90,fontsize=24):
    """
    Figura compuesta de UN sitio, 2x2:
      - columna izquierda: plot_species_turnover -- arriba PLANTAS
        (P2/P1/PP), abajo POLINIZADORES (A2/A1/PA).
      - columna derecha: plot_diverging_bars_with_constant -- arriba
        constante/aumenta/disminuye, abajo perdida/ganancia.
 
    Las DOS columnas comparten el MISMO height_ratios: el panel de abajo
    (polinizadores / divergente) ocupa mas alto que el de arriba (plantas /
    constante), en ambas columnas a la vez.
 
    Parameters
    ----------
    df_site_plants, df_site_pollinators : pd.DataFrame
        Tablas de turnover de especies, ya filtradas a UN Site_ID.
    df_site_interactions : pd.DataFrame
        Tabla de cambios de interacciones (Changes_df), ya filtrada a UN
        Site_ID.
    height_ratios : (float, float)
        Alto relativo arriba/abajo, aplicado a AMBAS columnas.
    plant_kwargs, pollinator_kwargs, diverging_kwargs : dict o None
        Parametros extra para cada funcion.
 
    Devuelve
    -------
    fig, (ax_plants, ax_poll, ax_top, ax_bottom)
    """
    if plant_kwargs is None:
        plant_kwargs = {}
    if pollinator_kwargs is None:
        pollinator_kwargs = {}
    if diverging_kwargs is None:
        diverging_kwargs = {}

    years=df_site_interactions.index.get_level_values('year').unique()
 
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=height_ratios,hspace=0.15, wspace=0.35)
 
    ax_plants = fig.add_subplot(gs[0, 0])
    ax_plants.tick_params(axis='x', labelbottom=False)
    ax_poll = fig.add_subplot(gs[1, 0], sharex=ax_plants)
    ax_top = fig.add_subplot(gs[0, 1])
    ax_top.tick_params(axis='x', labelbottom=False)
    ax_bottom = fig.add_subplot(gs[1, 1], sharex=ax_top)
 
    plot_species_turnover(df_site_species, ax=ax_plants,in_metric='P2', out_metric='P1', permanent_metric='PP')
    plot_species_turnover(df_site_species, ax=ax_poll,in_metric='A2', out_metric='A1', permanent_metric='PA')
 
    plot_diverging_bars_with_constant(df_site_interactions, ax_top=ax_top, ax_bottom=ax_bottom,height_ratios=height_ratios, **diverging_kwargs)
 
    ax_plants.set_title('Species', fontsize=fontsize+2)
    ax_top.set_title('Interactions', fontsize=fontsize+2)

    ax_bottom.set_xticklabels([year_label_dict.get(y, y) for y in years], rotation=rotate_years)
    ax_poll.set_xticklabels([year_label_dict.get(y, y) for y in years], rotation=rotate_years)

    ax_plants.set_ylabel('Plants',fontsize=fontsize-2)
    ax_poll.set_ylabel('Pollinators',fontsize=fontsize-2)
    ax_poll.set_xlabel('year',fontsize=fontsize)

    #legend
    if ax_plants.get_legend():
        ax_plants.get_legend().remove()
    if ax_poll.get_legend():
        ax_poll.get_legend().remove()
    handles, labels = ax_plants.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(-0.02, 0.5), frameon=False,fontsize=fontsize)

    fig.subplots_adjust(left=0.2)  # sube este valor si aun se corta

    ax_plants.text(0.02, 0.95, '\u273F', transform=ax_plants.transAxes,fontsize=fontsize, va='top', ha='left')
    ax_poll.text(0.02, 0.95, '\u2653', transform=ax_poll.transAxes, fontsize=fontsize, va='top', ha='left')
    
    plt.tight_layout()
 
    return fig, (ax_plants, ax_poll, ax_top, ax_bottom)
 
# final plot interaction changes trough time for all sites
def plot_diverging_bars_with_constant_all_sites(
        df, value_col='freq_sum', ncol=4, nrow=3, figsize=None,
        height_ratios_per_site=(1, 2), sharey='row', sharex=True,
        match_scale=True, ylim_top=None, ytick_step=20,
        loss_metrics=('CN_LI', 'L1N_LI', 'LBN_LI'),
        gain_metrics=('CN_GI', 'G1N_GI', 'GBN_GI'),
        loss_colors=('#0072B2', '#72B1D4','#C6DFEE' ),
        gain_colors=('#E69F00', '#F1CA72','#F9E9C6' ),
        constant_metrics=('CN_DI', 'CN_CI', 'CN_II'),
        constant_facecolor='#009455',
        constant_edgecolors=('#C6DBEF', None, '#FDD0A2'),
        constant_hatches=('///', None, '///'),
        hatch_linewidth=2.0,
        year_col='year', metric_col='metric', site_level='Site_ID',
        rotate_years=90, year_label_dict=None, label_dict=None,
        hatch_legend_order=(1, 2, 0), legend_ncol=3):
    """
    Version 'todos los sitios' de plot_diverging_bars_with_constant: para CADA
    Site_ID, dibuja los dos paneles (constante arriba + divergente abajo),
    organizados en una grilla nrow x ncol (por defecto 3x4, para 12 sitios).

    Replica los mismos parametros/colores/comportamiento que la version de
    un solo sitio (MultipleLocator para los ticks del eje y, year_label_dict
    para las etiquetas de año, orden del hatch en la leyenda) para que ambas
    figuras sean visualmente coherentes.

    Parameters
    ----------
    year_label_dict : dict o None
        {año: texto legible}. Si None, se usan los años tal cual.
    label_dict : dict o None
        {codigo_metrica: texto legible}, usado en la leyenda combinada.
    hatch_legend_order : tuple de int
        Orden de las 3 entradas del bloque "constante" en la leyenda, como
        indices sobre constant_metrics=('CN_DI','CN_CI','CN_II'). Por defecto
        (1, 2, 0) = CN_CI, CN_II, CN_DI (constante, aumenta, disminuye).
    sharey : {'row', 'all', None}
        Igual que en plot_diverging_bars_all_sites.
    ytick_step : float
        Paso de MultipleLocator para el eje y de TODOS los paneles.

    Devuelve
    -------
    fig, site_axes : dict {site: (ax_top, ax_bottom)}
    """
    if year_label_dict is None:
        year_label_dict = {}
    if label_dict is None:
        label_dict = {}

    df_flat = df.reset_index() if df.index.nlevels > 1 or df.index.name is not None else df.copy()
    sites = sorted(df_flat[site_level].unique())

    if figsize is None:
        figsize = (4.2 * ncol, (height_ratios_per_site[0] + height_ratios_per_site[1]) * 1.5 * nrow)

    fontsize=12

    fig = plt.figure(figsize=figsize)
    height_ratios = list(height_ratios_per_site) * nrow
    gs = gridspec.GridSpec(nrow * 2, ncol, figure=fig, height_ratios=height_ratios,
                            hspace=0.15, wspace=0.15)

    site_axes = {}
    ref_x_by_col = {}
    ref_ytop_by_row = {}
    ref_ybottom_by_row = {}
    ref_ytop_global = [None]
    ref_ybottom_global = [None]

    for i, site in enumerate(sites):
        row_group = i // ncol
        col = i % ncol

        sharex_ax = ref_x_by_col.get(col) if sharex else None
        if sharey == 'row':
            sharey_top_ax = ref_ytop_by_row.get(row_group)
            sharey_bottom_ax = ref_ybottom_by_row.get(row_group)
        elif sharey == 'all':
            sharey_top_ax = ref_ytop_global[0]
            sharey_bottom_ax = ref_ybottom_global[0]
        else:
            sharey_top_ax = None
            sharey_bottom_ax = None

        ax_top = fig.add_subplot(gs[row_group * 2, col], sharex=sharex_ax, sharey=sharey_top_ax)
        ax_bottom = fig.add_subplot(gs[row_group * 2 + 1, col], sharex=sharex_ax, sharey=sharey_bottom_ax)

        if sharex and col not in ref_x_by_col:
            ref_x_by_col[col] = ax_bottom
        if sharey == 'row':
            ref_ytop_by_row.setdefault(row_group, ax_top)
            ref_ybottom_by_row.setdefault(row_group, ax_bottom)
        elif sharey == 'all':
            if ref_ytop_global[0] is None:
                ref_ytop_global[0] = ax_top
            if ref_ybottom_global[0] is None:
                ref_ybottom_global[0] = ax_bottom

        df_site = df_flat[df_flat[site_level] == site]

        _plot_constant_on_ax(df_site, value_col, ax_top, constant_metrics,
                              constant_facecolor, constant_edgecolors,
                              constant_hatches, hatch_linewidth, year_col, metric_col)

        _plot_diverging_on_ax(df_site, value_col, ax_bottom, loss_metrics, gain_metrics,
                               loss_colors, gain_colors, year_col, metric_col)

        ax_top.set_title(str(site), fontsize=fontsize + 1)
        ax_top.set_xlabel('')
        ax_top.tick_params(axis='both', labelsize=fontsize)
        ax_bottom.tick_params(axis='both', labelsize=fontsize)

        years = sorted(df_site[year_col].unique())
        ax_bottom.set_xticks(np.arange(len(years)))
        ax_bottom.set_xticklabels([year_label_dict.get(y, y) for y in years],
                                   rotation=rotate_years, fontsize=fontsize)

        site_axes[site] = (ax_top, ax_bottom)

    # ylim del panel superior, DESPUES de que sharey sincronizo los inferiores
    for site, (ax_top, ax_bottom) in site_axes.items():
        if ylim_top is not None:
            ax_top.set_ylim(*ylim_top)
        elif match_scale:
            ax_top.set_ylim(0, ax_bottom.get_ylim()[1])
        ax_top.yaxis.set_major_locator(MultipleLocator(ytick_step))
        ax_bottom.yaxis.set_major_locator(MultipleLocator(ytick_step))

    # oculta labels interiores (deja solo bordes: año abajo, valores a la izq)
    for site, (ax_top, ax_bottom) in site_axes.items():
        ax_top.label_outer()
        ax_bottom.label_outer()

    # ylabels solo en la columna izquierda (los paneles ya comparten via label_outer,
    # pero el texto del label hay que ponerlo a mano en la primera columna de cada fila)
    for row_group in range(nrow):
        first_site_idx = row_group * ncol
        if first_site_idx < len(sites):
            ax_top, ax_bottom = site_axes[sites[first_site_idx]]
            ax_top.set_ylabel('Permanent\ninteractions', fontsize=fontsize-1)
            ax_bottom.set_ylabel('Transient\ninteractions', fontsize=fontsize-1)

    # leyenda UNICA combinada, debajo de toda la figura, con el orden del hatch aplicado
    handles = []
    ordered_constant = [(constant_metrics[j], constant_edgecolors[j], constant_hatches[j]) for j in hatch_legend_order]
    for metric, edgecolor, hatch in ordered_constant:
        handles.append(mpatches.Patch(
            facecolor=constant_facecolor,
            edgecolor=edgecolor if edgecolor else constant_facecolor,
            hatch=hatch, linewidth=0, label=label_dict.get(metric, metric)))
    for metric, color in zip(loss_metrics, loss_colors):
        handles.append(mpatches.Patch(color=color, label=label_dict.get(metric, metric)))
    for metric, color in zip(gain_metrics, gain_colors):
        handles.append(mpatches.Patch(color=color, label=label_dict.get(metric, metric)))

    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.0),
               ncol=legend_ncol, frameon=False,fontsize=fontsize)

    fig.supxlabel(year_col)

    return fig, site_axes

# final plot species changes trouhg tikme for all sites
def plot_species_turnover_all_sites(
        Nodes_df, ncol=4, nrow=3, figsize=None,
        height_ratios_per_site=(1, 2), site_level='Site_ID',
        year_col='year', metric_col='metric', value_col='n_nodes',
        plant_kwargs=None, pollinator_kwargs=None,
        fontsize=8, rotate_labels=90,
        inner_hspace=0.08, outer_hspace=0.5, wspace=0.3):
    """
    Version 'todos los sitios': para CADA Site_ID, dibuja los DOS paneles
    de plot_species_turnover (plantas arriba, polinizadores abajo).
 
    GridSpec ANIDADO: uno EXTERIOR controla el espacio ENTRE sitios
    (outer_hspace, mayor), y dentro de cada sitio, uno INTERIOR (2x1)
    controla el espacio entre plantas/polinizadores de ESE sitio
    (inner_hspace, menor) -- los dos espaciados se ajustan por separado.
 
    Parameters
    ----------
    height_ratios_per_site : (float, float)
        Alto relativo plantas/polinizadores DENTRO de cada sitio. Por
        defecto (1, 2).
    inner_hspace : float
        Espacio vertical entre plantas y polinizadores del MISMO sitio.
    outer_hspace : float
        Espacio vertical entre una fila de sitios y la siguiente.
    wspace : float
        Espacio horizontal entre columnas de sitios.
 
    Devuelve
    -------
    fig, site_axes : dict {site: (ax_plants, ax_poll)}
    """
    if plant_kwargs is None:
        plant_kwargs = {}
    if pollinator_kwargs is None:
        pollinator_kwargs = {}
 
    df_flat = Nodes_df.reset_index() if Nodes_df.index.nlevels > 1 or Nodes_df.index.name is not None else Nodes_df.copy()
    sites = sorted(df_flat[site_level].unique())
 
    if figsize is None:
        figsize = (3.6 * ncol, (sum(height_ratios_per_site)) * 1.6 * nrow)
 
    fig = plt.figure(figsize=figsize)
    outer_gs = gridspec.GridSpec(nrow, ncol, figure=fig,
                                  hspace=outer_hspace, wspace=wspace)
 
    site_axes = {}
    ref_x_by_col = {}
 
    for i, site in enumerate(sites):
        row_group = i // ncol
        col = i % ncol
 
        inner_gs = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_gs[row_group, col],
            height_ratios=height_ratios_per_site, hspace=inner_hspace)
 
        sharex_ax = ref_x_by_col.get(col)
 
        ax_plants = fig.add_subplot(inner_gs[0], sharex=sharex_ax)
        ax_poll = fig.add_subplot(inner_gs[1], sharex=ax_plants)
 
        if col not in ref_x_by_col:
            ref_x_by_col[col] = ax_poll
 
        df_site = df_flat[df_flat[site_level] == site]
 
        plot_species_turnover(df_site, ax=ax_plants,
                               in_metric='P2', out_metric='P1', permanent_metric='PP',
                               year_col=year_col, metric_col=metric_col, value_col=value_col,
                               fontsize=fontsize, show_legend=False, **plant_kwargs)
        plot_species_turnover(df_site, ax=ax_poll,
                               in_metric='A2', out_metric='A1', permanent_metric='PA',
                               year_col=year_col, metric_col=metric_col, value_col=value_col,
                               fontsize=fontsize, show_legend=False, **pollinator_kwargs)
 
        ax_plants.set_title(str(site), fontsize=fontsize + 2)
        ax_plants.tick_params(axis='x', labelbottom=False)
 
        site_axes[site] = (ax_plants, ax_poll)
 
        years = sorted(df_site[year_col].unique())
        ax_poll.set_xticks(np.arange(len(years)))
        ax_poll.set_xticklabels(years, rotation=rotate_labels, fontsize=fontsize)

        ax_poll.set_ylabel("Polliantors")
        ax_plants.set_ylabel("Plants")
 
    first_ax = site_axes[sites[0]][1]
    handles, labels = first_ax.get_legend_handles_labels()
    if not handles:
        first_ax = site_axes[sites[0]][0]
        handles, labels = first_ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0),
               ncol=4, frameon=False, fontsize=fontsize + 1)
 
    return fig, site_axes
 


# -----------------------------------------------------------------------
# helpers internos -- dibujan sobre un ax dado, sin crear figura/leyenda propia
# -----------------------------------------------------------------------

def _plot_constant_on_ax(df_site, value_col, ax, constant_metrics, constant_facecolor,
                          constant_edgecolors, constant_hatches, hatch_linewidth,
                          year_col, metric_col):
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = hatch_linewidth

    years = sorted(df_site[year_col].unique())
    x = np.arange(len(years))
    bottom = np.zeros(len(years))
    for metric, edgecolor, hatch in zip(constant_metrics, constant_edgecolors, constant_hatches):
        vals = np.array([
            df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
        ax.bar(x, vals, bottom=bottom, facecolor=constant_facecolor,
               edgecolor=edgecolor if edgecolor else constant_facecolor,
               hatch=hatch, linewidth=0, width=0.6)
        bottom += vals
    ax.set_xticks(x)


def _plot_diverging_on_ax(df_site, value_col, ax, loss_metrics, gain_metrics,
                           loss_colors, gain_colors, year_col, metric_col):
    years = sorted(df_site[year_col].unique())
    x = np.arange(len(years))

    bottom = np.zeros(len(years))
    for metric, color in zip(loss_metrics, loss_colors):
        vals = np.array([
            -df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
        ax.bar(x, vals, bottom=bottom, color=color, width=0.6)
        bottom += vals

    top = np.zeros(len(years))
    for metric, color in zip(gain_metrics, gain_colors):
        vals = np.array([
            df_site.loc[(df_site[year_col] == yr) & (df_site[metric_col] == metric), value_col].sum()
            for yr in years
        ])
        ax.bar(x, vals, bottom=top, color=color, width=0.6)
        top += vals

    ax.axhline(0, color='black', lw=1)
    ax.set_xticks(x)

# Plots showing average trough time, for one site
def plot_species_summary_single(
        df_site, value_col='n_nodes', ax=None, figsize=(4.5, 7),
        year_col='year', metric_col='metric',
        plant_permanent_metric='PP', poll_permanent_metric='PA',
        plant_in_metric='P2', plant_out_metric='P1',
        poll_in_metric='A2', poll_out_metric='A1',
        plant_color='#009455', poll_color='#009455',
        plant_in_color='#E69F00', plant_out_color='#0072B2',
        poll_in_color='#E69F00', poll_out_color='#0072B2',
        within_gap=0.65, between_gap=1.5,
        bar_width=0.5, updown_offset=0.08, capsize=4, fontsize=14,
        show_points=True, point_color='black', point_alpha=0.5,
        point_size=14, jitter=0.06,
        plant_icon='\u273F', animal_icon='\u2653',
        icon_fontsize=30, icon_y=-0.03):
    """
    Version de UN SOLO PANEL de plot_species_summary: las 4 categorias
    (Permanent Plants, Transient Plants, Permanent Animals, Transient Animals)
    en el MISMO eje x, en vez de dos paneles separados.

    Espaciado: Permanent-Transient DENTRO del mismo grupo (plantas o animales)
    quedan mas cerca entre si (within_gap); el hueco ENTRE el grupo de plantas
    y el de animales es mayor (between_gap) -- para que se lean como dos
    "bloques" claramente separados.

    Cada categoria "Permanent" es una barra simple (media +- std a traves de
    los años). Cada categoria "Transient" es una barra divergente (entra
    hacia arriba, sale hacia abajo, ligeramente desplazadas entre si via
    updown_offset), igual que en el panel inferior de plot_species_summary.

    Parameters
    ----------
    df_site : pd.DataFrame
        Ya filtrado a UN Site_ID, con year/metric/value_col.
    within_gap : float
        Distancia en x entre la barra "Permanent" y la "Transient" del MISMO
        grupo (plantas o animales).
    between_gap : float
        Distancia en x entre el grupo de plantas y el grupo de animales
        (deberia ser mayor que within_gap).
    show_points : bool
        Superpone los valores individuales (un punto por año) con jitter.
    plant_icon, animal_icon : str
        Caracter Unicode a mostrar debajo de cada bloque de etiquetas.
        Por defecto plant_icon='\u273F', animal_icon='\u2653'
    icon_fontsize : float
        Tamaño de los iconos.
    icon_y : float
        Posicion vertical de los iconos, en fraccion del eje (coordenadas
        de ejes, no de datos) -- mas negativo = mas abajo, para dejar
        hueco a las xticklabels que ya estan debajo del eje.

    Devuelve
    -------
    ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df_flat = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()
    rng = np.random.default_rng(0)

    def vals(metric):
        return df_flat.loc[df_flat[metric_col] == metric, value_col]

    def scatter_points(x0, values, sign=1):
        if not show_points or len(values) == 0:
            return
        jitter_x = x0 + rng.uniform(-jitter, jitter, size=len(values))
        ax.scatter(jitter_x, sign * values.values, facecolors='none', edgecolors=point_color,
                   alpha=point_alpha, s=point_size, linewidths=0.8, zorder=4)

    # ------------------------------------------------------------------
    # posiciones en x: PP -- (within_gap) -- transient_plants -- (between_gap) -- PA -- (within_gap) -- transient_animals
    # ------------------------------------------------------------------
    x_pp = 0
    x_pt = x_pp + within_gap
    x_pa = x_pt + between_gap
    x_at = x_pa + within_gap

    # --- Permanent Plants ---
    pp_vals = vals(plant_permanent_metric)
    ax.bar(x_pp, pp_vals.mean(), yerr=pp_vals.std(), color=plant_color, width=bar_width,
           capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
    scatter_points(x_pp, pp_vals)

    # --- Transient Plants (P2 arriba, P1 abajo) ---
    p2_vals = vals(plant_in_metric)
    p1_vals = vals(plant_out_metric)
    x_p2 = x_pt + updown_offset
    x_p1 = x_pt - updown_offset
    ax.bar(x_p2, p2_vals.mean(), yerr=p2_vals.std(), color=plant_in_color, width=bar_width,
           capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
    ax.bar(x_p1, -p1_vals.mean(), yerr=p1_vals.std(), color=plant_out_color, width=bar_width,
           capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
    scatter_points(x_p2, p2_vals, sign=1)
    scatter_points(x_p1, p1_vals, sign=-1)

    # --- Permanent Animals ---
    pa_vals = vals(poll_permanent_metric)
    ax.bar(x_pa, pa_vals.mean(), yerr=pa_vals.std(), color=poll_color, width=bar_width,
           capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
    scatter_points(x_pa, pa_vals)

    # --- Transient Animals (A2 arriba, A1 abajo) ---
    a2_vals = vals(poll_in_metric)
    a1_vals = vals(poll_out_metric)
    x_a2 = x_at + updown_offset
    x_a1 = x_at - updown_offset
    ax.bar(x_a2, a2_vals.mean(), yerr=a2_vals.std(), color=poll_in_color, width=bar_width,
           capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
    ax.bar(x_a1, -a1_vals.mean(), yerr=a1_vals.std(), color=poll_out_color, width=bar_width,
           capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
    scatter_points(x_a2, a2_vals, sign=1)
    scatter_points(x_a1, a1_vals, sign=-1)

    ax.axhline(0, color='black', lw=1)
    ax.set_xticks([x_pp, x_pt, x_pa, x_at])
    ax.set_xticklabels(['P', 'T',
                         'P', 'T'], fontsize=fontsize)
    ax.set_ylabel("Species richness", fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)

    # iconos debajo de cada bloque (planta bajo el par Permanent/Transient
    # Plants, animal bajo el par Permanent/Transient Animals), centrados.
    # transform=get_xaxis_transform() -- x en coordenadas de DATOS, y en
    # fraccion del eje (0=borde inferior del axes), asi el icono se coloca
    # debajo de las xticklabels independientemente del rango de datos.
    x_plant_center = (x_pp + x_pt) / 2
    x_animal_center = (x_pa + x_at) / 2
    ax.text(x_plant_center, icon_y, plant_icon, fontsize=icon_fontsize,
            ha='center', va='top', transform=ax.get_xaxis_transform())
    ax.text(x_animal_center, icon_y, animal_icon, fontsize=icon_fontsize,
            ha='center', va='top', transform=ax.get_xaxis_transform())


    return ax

def plot_diverging_bars_mean_error(df, value_col='freq_sum', ax=None, figsize=(5, 5),
                                    loss_metrics=('CN_LI', 'L1N_LI', 'LBN_LI'),
                                    gain_metrics=('CN_GI', 'G1N_GI', 'GBN_GI'),
                                    category_labels=('Permanent-\nPermanent', 'Permanent-\nTransient', 'Transient-\nTransient'),
                                    loss_colors=('#C6DFEE', '#72B1D4', '#0072B2'),
                                    gain_colors=('#F9E9C6', '#F1CA72', '#E69F00'),
                                    year_col='year', metric_col='metric', site_level='Site_ID',
                                    bar_width=0.5, updown_offset=0.15, capsize=4, fontsize=10,
                                    show_points=True, point_color='black',
                                    point_alpha=0.35, point_size=12, jitter=0.1):
    """
    Grafico de barras con barra de error. Eje x: las 3 categorias de nodo
    (Permanent-Permanent / Permanent-Transient / Transient-Transient), NO los
    años. Para cada categoria hay DOS barras (perdida hacia abajo, ganancia
    hacia arriba), SIN apilar entre si -- y ligeramente DESPLAZADAS entre si
    (perdida a la izquierda, ganancia a la derecha) via updown_offset, para
    que no queden perfectamente alineadas.

    Altura de cada barra = media de value_col a traves de TODAS las filas de
    df para esa metrica (si pasas el df completo, sitios x años; si pasas
    df de un solo sitio, la media es solo a traves de los años de ese sitio).
    Barra de error = desviacion estandar de esa misma coleccion.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame con Site_ID/year/metric como MultiIndex o columnas. Puede
        ser el df completo (todos los sitios) o ya filtrado a un sitio.
    loss_metrics, gain_metrics : tuple de 3 str
        En el orden (CN, 1-nuevo, ambas-nuevas).
    category_labels : tuple de 3 str
        Etiquetas del eje x, mismo orden que loss_metrics/gain_metrics.
    loss_colors, gain_colors : tuple de 3 str
        Un color por categoria (mismo orden que loss_metrics/gain_metrics).
    updown_offset : float
        Desplazamiento horizontal entre la barra de perdida (se mueve a la
        IZQUIERDA) y la de ganancia (se mueve a la DERECHA) dentro de cada
        categoria. 0 = sin desplazamiento (alineadas).
    show_points : bool
        Si True (por defecto), superpone los valores individuales con jitter.

    Devuelve
    -------
    ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df_flat = df.reset_index() if df.index.nlevels > 1 or df.index.name is not None else df.copy()

    x = np.arange(len(loss_metrics))

    def mean_std_vals(metric):
        vals = df_flat.loc[df_flat[metric_col] == metric, value_col]
        return vals.mean(), vals.std(), vals.values

    loss_stats = [mean_std_vals(m) for m in loss_metrics]
    gain_stats = [mean_std_vals(m) for m in gain_metrics]

    rng = np.random.default_rng(0)  # semilla fija -- jitter reproducible entre llamadas

    for i in range(len(x)):
        loss_mean, loss_std, loss_vals = loss_stats[i]
        gain_mean, gain_std, gain_vals = gain_stats[i]

        x_loss = x[i] - updown_offset
        x_gain = x[i] + updown_offset

        ax.bar(x_loss, -loss_mean, yerr=loss_std, color=loss_colors[i],
               width=bar_width, capsize=capsize,
               label='Lost' if i == 0 else None,
               error_kw={'elinewidth': 1.2}, zorder=2)
        ax.bar(x_gain, gain_mean, yerr=gain_std, color=gain_colors[i],
               width=bar_width, capsize=capsize,
               label='Gained' if i == 0 else None,
               error_kw={'elinewidth': 1.2}, zorder=2)

        if show_points:#ax_top.scatter(jitter_x, m_vals, facecolors='none', edgecolors=point_color, alpha=point_alpha, s=point_size, linewidths=0.8, zorder=3)
            jitter_x_loss = x_loss + rng.uniform(-jitter, jitter, size=len(loss_vals))
            jitter_x_gain = x_gain + rng.uniform(-jitter, jitter, size=len(gain_vals))
            ax.scatter(jitter_x_loss, -loss_vals, facecolors='none', edgecolors=point_color, alpha=point_alpha, s=point_size, linewidths=0.8, zorder=3)
            ax.scatter(jitter_x_gain, gain_vals, facecolors='none', edgecolors=point_color, alpha=point_alpha, s=point_size, linewidths=0.8, zorder=3)

    ax.axhline(0, color='black', lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(category_labels, fontsize=fontsize)
    ax.set_ylabel(value_col, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.legend(frameon=False, fontsize=fontsize)

    return ax

def plot_interactions_summary(
        df_site, value_col='freq_sum', figsize=(5, 7), height_ratios=(1, 1.4),
        year_col='year', metric_col='metric',
        constant_metrics=('CN_CI', 'CN_DI', 'CN_II'),
        constant_facecolor='#009455',
        constant_edgecolors=(None, '#C6DFEE', '#F9E9C6'),
        constant_hatches=(None, '\\\\\\', '///'),
        hatch_linewidth=3.0, top_bar_width=0.4,
        loss_metrics=('CN_LI', 'L1N_LI', 'LBN_LI'),
        gain_metrics=('CN_GI', 'G1N_GI', 'GBN_GI'),
        category_labels=('Permanent-\nPermanent\n(Rewiring)', 'Permanent-\nTransient', 'Transient-\nTransient'),
        loss_colors=('#0072B2', '#72B1D4', '#C6DFEE'),
        gain_colors=('#E69F00', '#F1CA72', '#F9E9C6'),
        bottom_bar_width=0.4, updown_offset=0.15, capsize=4, fontsize=14,
        show_points=True, point_color='black', point_alpha=0.8,
        point_size=12, jitter=0.1, ax_top=None, ax_bottom=None):
    """
    Figura de UN sitio (df_site) con DOS paneles, para las INTERACCIONES
    (a diferencia de plot_species_summary_single, que es para especies/nodos):
 
      - panel superior: 3 barras SIN apilar -- CN_CI (verde liso), CN_DI
        (verde + rayas azules), CN_II (verde + rayas amarillas) -- media +-
        std a traves de los AÑOS de este sitio (df_site ya viene de un
        unico Site_ID).
      - panel inferior: plot_diverging_bars_mean_error aplicada a df_site
        -- como solo se le pasa un sitio, la media +- std tambien es solo
        a traves de los años de ESE sitio. Perdida (azul) desplazada a la
        izquierda, ganancia (amarillo) desplazada a la derecha, dentro de
        cada una de las 3 categorias de nodo (Permanent-Permanent /
        Permanent-Transient / Transient-Transient).
 
    Parameters
    ----------
    df_site : pd.DataFrame
        Ya filtrado a UN Site_ID, con year/metric/value_col para todos los
        años de ese sitio.
    constant_metrics : tuple de 3 str
        (constante, disminuye, aumenta) -- CN_CI/CN_DI/CN_II por defecto.
    constant_edgecolors, constant_hatches : tuple de 3 (str o None)
        Estilo de cada barra del panel superior, mismo orden que
        constant_metrics. None = sin rayas (solo el color solido).
    ax_top, ax_bottom : matplotlib Axes o None
        Si se pasan AMBOS, se dibuja sobre esos ejes en vez de crear una
        figura nueva -- util para componer esta figura dentro de otra mas
        grande (p.ej. junto a plot_species_summary_single). Si se dejan en
        None (por defecto), se crea una figura nueva de 2 paneles como
        siempre.
 
    Devuelve
    -------
    Si ax_top/ax_bottom son None: fig, (ax_top, ax_bottom)
    Si se pasaron ax_top/ax_bottom: (ax_top, ax_bottom) -- sin fig, porque
    la figura ya existia de antes (la creo quien llama).
    """
    df_flat = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()
 
    def mean_std_vals(metric):
        vals = df_flat.loc[df_flat[metric_col] == metric, value_col]
        return vals.mean(), vals.std(), vals.values
 
    created_own_fig = ax_top is None and ax_bottom is None
    if created_own_fig:
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=figsize,
                                                 gridspec_kw={'height_ratios': height_ratios})
 
    # ------------------------------------------------------------------
    # panel superior: CN_CI / CN_DI / CN_II, SIN apilar, media +- std sobre años
    # ------------------------------------------------------------------
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = hatch_linewidth
 
    rng = np.random.default_rng(0)
    x_top = np.arange(len(constant_metrics))
 
    for i, metric in enumerate(constant_metrics):
        m_mean, m_std, m_vals = mean_std_vals(metric)
        edgecolor = constant_edgecolors[i]
        hatch = constant_hatches[i]
        ax_top.bar(x_top[i], m_mean, yerr=m_std, facecolor=constant_facecolor,
                   edgecolor=edgecolor if edgecolor else constant_facecolor,
                   hatch=hatch, linewidth=0, width=top_bar_width,
                   capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
        if show_points:
            jitter_x = x_top[i] + rng.uniform(-jitter, jitter, size=len(m_vals))
            ax_top.scatter(jitter_x, m_vals, facecolors='none', edgecolors=point_color, alpha=point_alpha, s=point_size, linewidths=0.8, zorder=3)
                           #ax.scatter(jitter_x, sign * values.values, facecolors='none', edgecolors=point_color, alpha=point_alpha, s=point_size, linewidths=0.8, zorder=4)
 
    ax_top.set_xticks(x_top)
    ax_top.set_xticklabels(['Constant','Decreasing', 'Increasing'], fontsize=fontsize)
    ax_top.set_ylabel(value_col, fontsize=fontsize)
    ax_top.tick_params(axis='y', labelsize=fontsize)
    ax_top.set_title('Permanent-Permanent: constant / decreasing / increasing', fontsize=fontsize + 1)
 
    # ------------------------------------------------------------------
    # panel inferior: plot_diverging_bars_mean_error, reutilizada tal cual
    # ------------------------------------------------------------------
    plot_diverging_bars_mean_error(
        df_flat, value_col=value_col, ax=ax_bottom,
        loss_metrics=loss_metrics, gain_metrics=gain_metrics,
        category_labels=category_labels, loss_colors=loss_colors, gain_colors=gain_colors,
        year_col=year_col, metric_col=metric_col,
        bar_width=bottom_bar_width, updown_offset=updown_offset, capsize=capsize,
        fontsize=fontsize, show_points=show_points, point_color=point_color,
        point_alpha=point_alpha, point_size=point_size, jitter=jitter)
 
    # en vez de la leyenda automatica (arriba a la izq), colocamos "Lost
    # interactions" abajo a la derecha (junto a las barras de perdida) y
    # "Gained interactions" arriba a la derecha (junto a las de ganancia)
    legend = ax_bottom.get_legend()
    if legend is not None:
        legend.remove()
    ax_bottom.text(0.98, 0.02, 'Lost', transform=ax_bottom.transAxes,
                   ha='right', va='bottom', fontsize=fontsize)
    ax_bottom.text(0.98, 0.98, 'Gained', transform=ax_bottom.transAxes,
                   ha='right', va='top', fontsize=fontsize)
    
    
    if created_own_fig:
        plt.tight_layout()
        return fig, (ax_top, ax_bottom)
    else:
        return ax_top, ax_bottom


def plot_interaction_frequency_bars(
        df_site, freq_col='freq_sum', n_col='n_interactions', ax=None, figsize=(6, 7),
        metrics=('CN_CI', 'CN_DI', 'CN_II', 'CN_GI', 'G1N_GI', 'GBN_GI', 'CN_LI', 'L1N_LI', 'LBN_LI'),
        metric_col='metric', year_col='year',
        colors=('#009455', '#009455', '#009455',
                '#E69F00', '#F1CA72', '#F9E9C6',
                '#0072B2', '#72B1D4', '#C6DFEE' ),
        edgecolors=(None, '#C6DFEE', '#F9E9C6', None, None, None, None, None, None),
        hatches=(None, '\\\\\\', '///', None, None, None, None, None, None),
        hatch_linewidth=3.0,
        label_dict=label_dict2, bar_height=0.6, capsize=4, fontsize=14,
        show_points=True, point_color='black', point_alpha=0.35,
        point_size=14, jitter=0.08):
    """
    Grafico de barras HORIZONTALES: frecuencia media por interaccion
    (freq_sum/n_interactions) de cada una de las 9 categorias, para UN sitio.
 
    Orden de arriba a abajo (por defecto): CN_CI, CN_DI, CN_II -- bloque
    permanente; CN_GI, G1N_GI, GBN_GI -- ganancias; CN_LI, L1N_LI, LBN_LI --
    perdidas.
 
    CN_DI y CN_II llevan el MISMO estilo de rayas (hatch) que en
    plot_interactions_summary: verde de fondo + rayas azules (CN_DI) o
    amarillas (CN_II), en vez de color solido.
 
    Eje y: las 9 categorias (metrics), con su etiqueta a la IZQUIERDA (layout
    estandar de barh). Eje x: frecuencia media por interaccion. Cada fila del
    df_site (un año, una categoria) aporta un valor freq_col/n_col para los
    puntos/barra de error; la barra en si es la razon AGREGADA (suma de
    freq_col entre suma de n_col, sumando todos los años).
 
    Parameters
    ----------
    df_site : pd.DataFrame
        Ya filtrado a UN Site_ID -- Changes_df.loc[site].
    freq_col, n_col : str
        Columnas a dividir (freq_col / n_col).
    metrics : tuple de 9 str
        Orden de las categorias de arriba a abajo en el eje y.
    colors : tuple de 9 str
        Color de fondo (facecolor) de cada categoria, mismo orden que metrics.
    edgecolors, hatches : tuple de 9 (str o None)
        Estilo de rayas por categoria, mismo orden que metrics. None = sin
        rayas (barra solida, usa solo el color de colors).
    hatch_linewidth : float
        Grosor de las rayas del hatch.
    label_dict : dict o None
        {codigo_metrica: texto legible} para las etiquetas del eje y.
 
    Devuelve
    -------
    ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if label_dict is None:
        label_dict = {}
 
    df_flat = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()
    rng = np.random.default_rng(0)
 
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = hatch_linewidth
 
    y = np.arange(len(metrics))
 
    means, stds, all_vals = [], [], []
    for metric in metrics:
        sub = df_flat.loc[df_flat[metric_col] == metric]
        vals = (sub[freq_col] / sub[n_col]).replace([np.inf, -np.inf], np.nan).dropna()
        all_vals.append(vals.values)
        stds.append(vals.std())
 
        # razon AGREGADA (suma de freq_col entre suma de n_col, sumando
        # todos los años) -- no es la media de las razones por año
        total_freq = sub[freq_col].sum()
        total_n = sub[n_col].sum()
        means.append(total_freq / total_n if total_n != 0 else np.nan)
 
    for i in range(len(metrics)):
        edgecolor = edgecolors[i]
        ax.barh(y[i], means[i], xerr=stds[i], facecolor=colors[i],
                edgecolor=edgecolor if edgecolor else colors[i],
                hatch=hatches[i], linewidth=0, height=bar_height,
                capsize=capsize, error_kw={'elinewidth': 1.2}, zorder=2)
 
    if show_points:
        for i, vals in enumerate(all_vals):
            jitter_y = y[i] + rng.uniform(-jitter, jitter, size=len(vals))
            ax.scatter(vals, jitter_y, facecolors='none', edgecolors=point_color,
                       alpha=point_alpha, s=point_size, linewidths=0.8, zorder=3)
 
    ax.set_yticks(y)
    ax.set_yticklabels([label_dict.get(m, m) for m in metrics], fontsize=fontsize)
    ax.yaxis.tick_right()
    ax.invert_yaxis()  # primera categoria de la lista arriba, no abajo
    ax.set_xlim(0, None)
    ax.invert_xaxis()  # primera categoria de la lista arriba, no abajo
    ax.set_xlabel('Average interaction frequency', fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.axvline(0, color='black', lw=1)
 
    return ax
 
def plot_species_and_interactions_summary(
        node_changes_df_site, changes_df_site,
        figsize=(15, 7), width_ratios=(1, 1.4, 1), height_ratios=(1, 1.4),
        species_kwargs=None, interactions_kwargs=None, frequency_kwargs=None):
    """
    Figura compuesta de UN sitio: plot_species_summary_single a la IZQUIERDA
    (un unico panel), plot_interactions_summary en el CENTRO (dos paneles
    apilados), y plot_interaction_frequency_bars a la DERECHA (un unico
    panel, barras horizontales) -- reutilizando las tres funciones tal cual,
    sin duplicar logica.
 
    Parameters
    ----------
    node_changes_df_site : pd.DataFrame
        Ya filtrado a UN Site_ID -- tabla de turnover de ESPECIES/NODOS
        (PP/PA/P1/P2/A1/A2 en metric_col), la misma que usas para
        plot_species_summary_single (tu "node_changes_df").
    changes_df_site : pd.DataFrame
        Ya filtrado a UN Site_ID -- tabla de cambios de INTERACCIONES
        (CN_CI/CN_DI/CN_II/CN_LI/CN_GI/L1N_LI/LBN_LI/G1N_GI/GBN_GI en
        metric_col), la misma que usas para plot_interactions_summary
        (tu "Changes_df"). Se reutiliza TAMBIEN para el panel derecho
        (plot_interaction_frequency_bars) -- misma tabla, otra vista.
    width_ratios : (float, float, float)
        Ancho relativo columna izquierda (species) / centro (interactions) /
        derecha (frequency bars).
    height_ratios : (float, float)
        Alto relativo panel superior/inferior de la columna CENTRAL (se pasa
        tal cual a plot_interactions_summary). Las columnas izquierda y
        derecha ocupan el alto TOTAL (ambas filas) cada una.
    species_kwargs : dict o None
        Parametros extra para plot_species_summary_single (ademas de
        node_changes_df_site y ax) -- p.ej. {'value_col': 'n_nodes', ...}.
    interactions_kwargs : dict o None
        Parametros extra para plot_interactions_summary (ademas de
        changes_df_site, ax_top, ax_bottom) -- p.ej. {'value_col': 'freq_sum', ...}.
    frequency_kwargs : dict o None
        Parametros extra para plot_interaction_frequency_bars (ademas de
        changes_df_site y ax) -- p.ej. {'freq_col': 'freq_sum', 'n_col': 'n_interactions', ...}.
 
    Devuelve
    -------
    fig, (ax_species, ax_int_top, ax_int_bottom, ax_freq)
    """
    if species_kwargs is None:
        species_kwargs = {}
    if interactions_kwargs is None:
        interactions_kwargs = {}
    if frequency_kwargs is None:
        frequency_kwargs = {}
 
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 3, figure=fig, width_ratios=width_ratios,
                            height_ratios=height_ratios, wspace=0.3, hspace=0.1)
 
    ax_species = fig.add_subplot(gs[:, 0])     # columna izquierda, ambas filas
    ax_int_top = fig.add_subplot(gs[0, 1])     # columna central, fila superior
    ax_int_bottom = fig.add_subplot(gs[1, 1])  # columna central, fila inferior
    ax_freq = fig.add_subplot(gs[:, 2])        # columna derecha, ambas filas
 
    plot_species_summary_single(node_changes_df_site, ax=ax_species, **species_kwargs)
    plot_interactions_summary(changes_df_site, ax_top=ax_int_top, ax_bottom=ax_int_bottom,
                               **interactions_kwargs)
    plot_interaction_frequency_bars(changes_df_site, ax=ax_freq, **frequency_kwargs)
 
    ax_species.set_title("Species", fontsize=18)
    ax_int_top.set_title("Interaction counts", fontsize=18)
    ax_int_top.set_ylabel("Permanent interactions")
    ax_int_bottom.set_ylabel("Transient interactions")
    ax_freq.set_title("Interaction frequency", fontsize=18)
 
    ymin, ymax = ax_int_bottom.get_ylim()
    m = max(abs(ymin), abs(ymax))
    ax_int_bottom.set_ylim(-m, m)
 
    ymin, ymax = ax_species.get_ylim()
    m = max(abs(ymin), abs(ymax))
    ax_species.set_ylim(-m, m)

 
    plt.tight_layout()
 
    return fig, (ax_species, ax_int_top, ax_int_bottom, ax_freq)
 

#statisticfal test 

from scipy.stats import kruskal
import scikit_posthocs as sp

def test_interaction_frequency_differences(
        df_site, freq_col='freq_interactions', n_col='n_interactions',
        metrics=('CN_CI', 'CN_DI', 'CN_II', 'CN_LI', 'CN_GI', 'L1N_LI', 'LBN_LI', 'G1N_GI', 'GBN_GI'),
        metric_col='metric', alpha=0.05):
    """
    Test estadistico para saber si la frecuencia media por interaccion
    (freq_col/n_col) difiere entre las 9 categorias.
 
    1. Kruskal-Wallis (no parametrico, equivalente a ANOVA de un factor):
       compara las 9 categorias A LA VEZ -- responde "¿hay AL MENOS una
       categoria que difiere de las demas?", sin decir cual.
       Se usa Kruskal-Wallis en vez de ANOVA porque con pocos años por sitio
       (tipicamente ~4 puntos por categoria) no se puede asumir normalidad
       de forma fiable -- Kruskal-Wallis no la requiere.
 
    2. Test de Dunn (post-hoc, con correccion de Bonferroni): SOLO si el
       Kruskal-Wallis sale significativo, compara TODOS los pares de
       categorias entre si, para saber CUALES son las que difieren.
 
    Parameters
    ----------
    df_site : pd.DataFrame
        Ya filtrado a UN Site_ID (o varios, si quieres agregar mas puntos),
        con metric_col/freq_col/n_col.
    metrics : tuple de str
        Categorias a comparar.
    alpha : float
        Umbral de significacion (0.05 por defecto).
 
    Devuelve
    -------
    kw_result : dict
        {'H': estadistico, 'p_value': float, 'significant': bool}
    posthoc_df : pd.DataFrame o None
        Matriz de p-valores (Dunn, Bonferroni) para cada par de categorias --
        None si el Kruskal-Wallis no fue significativo (no tiene sentido
        mirar los pares si no hay diferencia global).
    """
    df_flat = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()
 
    # una lista de valores (razon freq/n por año) por cada categoria
    groups = []
    labels = []
    for metric in metrics:
        sub = df_flat.loc[df_flat[metric_col] == metric]
        vals = (sub[freq_col] / sub[n_col]).replace([np.inf, -np.inf], np.nan).dropna()
        if len(vals) > 0:
            groups.append(vals.values)
            labels.append(metric)
 
    # ------------------------------------------------------------------
    # 1) Kruskal-Wallis -- diferencia global entre las categorias
    # ------------------------------------------------------------------
    H, p_value = kruskal(*groups)
    significant = p_value < alpha
 
    kw_result = {'H': H, 'p_value': p_value, 'significant': significant}
 
    print(f"Kruskal-Wallis: H={H:.3f}, p={p_value:.4f} "
          f"({'SIGNIFICATIVO' if significant else 'no significativo'} a alpha={alpha})")
 
    # ------------------------------------------------------------------
    # 2) Post-hoc de Dunn -- SOLO si el paso 1 fue significativo
    # ------------------------------------------------------------------
    posthoc_df = None
    if significant:
        # formato "largo" que pide scikit-posthocs: una fila por valor,
        # con su categoria correspondiente
        long_data = []
        for label, vals in zip(labels, groups):
            for v in vals:
                long_data.append({'metric': label, 'value': v})
        long_df = pd.DataFrame(long_data)
 
        posthoc_df = sp.posthoc_dunn(long_df, val_col='value', group_col='metric',
                                      p_adjust='holm')#'bonferroni')
        posthoc_df = posthoc_df.loc[labels, labels]  # reordena segun 'metrics'
 
        print("\nPost-hoc de Dunn (p-valores, correccion Holm/Bonferroni):")
        print(posthoc_df.round(4))
    else:
        print("Como el Kruskal-Wallis no fue significativo, no se hace el post-hoc "
              "(no habria diferencias que localizar).")
 
    return kw_result, posthoc_df
 
# TEMPORALITY 

def _compute_site_year_metrics(
        nodes_df_lag, changes_df_lag, site_level, metric_col, year_col,
        value_col_nodes, value_col_changes,
        permanent_int_categories, change_categories,
        poll_permanent, poll_transient, plant_permanent, plant_transient):
    """
    Helper interno: calcula, por sitio y año, las variables (proporcion de
    polinizadores permanentes, de plantas permanentes, numero total de
    cambios de interaccion, y proporcion de interacciones permanentes) para
    UN lag. Devuelve 4 dicts {site: [valores por año]}.
 
    "Permanente" (interaccion) = la interaccion SIGUE EXISTIENDO en ambos
    tiempos, cambie o no de intensidad -- CN_CI + CN_DI + CN_II.
    "Cambio" = la interaccion aparece o desaparece del todo -- el resto de
    categorias (CN_LI, CN_GI, L1N_LI, LBN_LI, G1N_GI, GBN_GI).
    """
    nodes_flat = nodes_df_lag.reset_index() if nodes_df_lag.index.nlevels > 1 or nodes_df_lag.index.name is not None else nodes_df_lag.copy()
    changes_flat = changes_df_lag.reset_index() if changes_df_lag.index.nlevels > 1 or changes_df_lag.index.name is not None else changes_df_lag.copy()
 
    sites = sorted(nodes_flat[site_level].unique())
 
    def get_node_val(df_site_year, metric):
        v = df_site_year.loc[df_site_year[metric_col] == metric, value_col_nodes]
        return v.sum() if len(v) else 0
 
    poll_ratios = {s: [] for s in sites}
    plant_ratios = {s: [] for s in sites}
    n_changes = {s: [] for s in sites}
    prop_permanent_int = {s: [] for s in sites}
 
    for site in sites:
        sub_nodes = nodes_flat[nodes_flat[site_level] == site]
        sub_changes = changes_flat[changes_flat[site_level] == site]
        years = sorted(sub_nodes[year_col].unique())
 
        for yr in years:
            yr_nodes = sub_nodes[sub_nodes[year_col] == yr]
 
            pa = get_node_val(yr_nodes, poll_permanent)
            a_trans = sum(get_node_val(yr_nodes, m) for m in poll_transient)
            poll_total = pa + a_trans
            poll_ratios[site].append(pa / poll_total if poll_total > 0 else np.nan)
 
            pp = get_node_val(yr_nodes, plant_permanent)
            p_trans = sum(get_node_val(yr_nodes, m) for m in plant_transient)
            plant_total = pp + p_trans
            plant_ratios[site].append(pp / plant_total if plant_total > 0 else np.nan)
 
            yr_changes = sub_changes[sub_changes[year_col] == yr]
            n_permanent_int = yr_changes.loc[yr_changes[metric_col].isin(permanent_int_categories), value_col_changes].sum()
            n_change = yr_changes.loc[yr_changes[metric_col].isin(change_categories), value_col_changes].sum()
            n_all_int = n_permanent_int + n_change
 
            n_changes[site].append(n_change)
            prop_permanent_int[site].append(n_permanent_int / n_all_int if n_all_int > 0 else np.nan)
 
    return poll_ratios, plant_ratios, n_changes, prop_permanent_int
 
#efor one site
from matplotlib.ticker import PercentFormatter
def plot_lag_summary_temporal(
        All_Nodes, All_Changes, site, lags=(1, 2, 3, 4), ax=None, figsize=(20, 4.5),
        site_level='Site_ID', metric_col='metric', year_col='year',
        value_col_nodes='n_nodes',
        value_col_changes_1='n_interactions', value_col_changes_2='freq_interactions',
        permanent_int_categories=('CN_CI', 'CN_DI', 'CN_II'),
        change_categories=('CN_LI', 'CN_GI', 'L1N_LI', 'LBN_LI', 'G1N_GI', 'GBN_GI'),
        poll_permanent='PA', poll_transient=('A1', 'A2'),
        plant_permanent='PP', plant_transient=('P1', 'P2'),
        cmap=None, cmap_range=(1.0, 0.4),
        marker='.', markersize=14, capsize=0, elinewidth=2.5,
        fontsize=20, sharey=True, rotate_labels=0,
        titles=('Permanent floral visitors', 'Permanent flowering plants','Permanent interactions', 'Permanent frequency')):
    """
    Version TEMPORAL de plot_persistence_by_distance -- caso particular de
    plot_lag_summary_combined_horizontal, filtrado SIEMPRE a UN sitio: el
    eje x son directamente los LAGS (t+1, t+2, t+3, t+4) de ESE sitio --
    un punto por lag, sin desplazamiento (ya no compiten varias series en
    la misma posicion x).
 
    4 paneles horizontales, mismo estilo que plot_persistence_by_distance.
 
    Parameters
    ----------
    All_Nodes, All_Changes : dict {lag: pd.DataFrame}
        Tus diccionarios completos.
    site : str
        El Site_ID a representar -- SIEMPRE se filtra a este sitio.
    lags : tuple de int
        Que lags mostrar, en el eje x.
    ax : lista de 4 matplotlib Axes o None
        Si se pasa, dibuja sobre esos 4 ejes. Si no, crea figura nueva.
 
    Devuelve
    -------
    fig (o None si se pasaron ejes), axes (lista de 4), summary_df
    """
    if cmap is None:
        cmap = plt.cm.Greens
    colors = cmap(np.linspace(cmap_range[0], cmap_range[1], len(lags)))
 
    All_Nodes_site = {}
    All_Changes_site = {}
    for lag in lags:
        nodes_flat = All_Nodes[lag]
        nodes_flat = nodes_flat.reset_index() if nodes_flat.index.nlevels > 1 or nodes_flat.index.name is not None else nodes_flat.copy()
        changes_flat = All_Changes[lag]
        changes_flat = changes_flat.reset_index() if changes_flat.index.nlevels > 1 or changes_flat.index.name is not None else changes_flat.copy()
 
        All_Nodes_site[lag] = nodes_flat[nodes_flat[site_level] == site]
        All_Changes_site[lag] = changes_flat[changes_flat[site_level] == site]
 
    plant_means, plant_stds = [], []
    poll_means, poll_stds = [], []
    prop1_means, prop1_stds = [], []
    prop2_means, prop2_stds = [], []
 
    for lag in lags:
        poll_r, plant_r, _, prop_1 = _compute_site_year_metrics(
            All_Nodes_site[lag], All_Changes_site[lag], site_level, metric_col, year_col,
            value_col_nodes, value_col_changes_1,
            permanent_int_categories, change_categories,
            poll_permanent, poll_transient, plant_permanent, plant_transient)
 
        _, _, _, prop_2 = _compute_site_year_metrics(
            All_Nodes_site[lag], All_Changes_site[lag], site_level, metric_col, year_col,
            value_col_nodes, value_col_changes_2,
            permanent_int_categories, change_categories,
            poll_permanent, poll_transient, plant_permanent, plant_transient)
 
        poll_vals = poll_r.get(site, [])
        plant_vals = plant_r.get(site, [])
        prop1_vals = prop_1.get(site, [])
        prop2_vals = prop_2.get(site, [])
 
        poll_means.append(np.nanmean(poll_vals) if len(poll_vals) else np.nan)
        poll_stds.append(np.nanstd(poll_vals) if len(poll_vals) else np.nan)
        plant_means.append(np.nanmean(plant_vals) if len(plant_vals) else np.nan)
        plant_stds.append(np.nanstd(plant_vals) if len(plant_vals) else np.nan)
        prop1_means.append(np.nanmean(prop1_vals) if len(prop1_vals) else np.nan)
        prop1_stds.append(np.nanstd(prop1_vals) if len(prop1_vals) else np.nan)
        prop2_means.append(np.nanmean(prop2_vals) if len(prop2_vals) else np.nan)
        prop2_stds.append(np.nanstd(prop2_vals) if len(prop2_vals) else np.nan)
 
    lag_labels = [f't+{lag}' for lag in lags]
    summary_df = pd.DataFrame({
        'lag': lag_labels,
        'plant_mean': plant_means, 'plant_std': plant_stds,
        'pollinator_mean': poll_means, 'pollinator_std': poll_stds,
        'interaction_1_mean': prop1_means, 'interaction_1_std': prop1_stds,
        'interaction_2_mean': prop2_means, 'interaction_2_std': prop2_stds,
    }).set_index('lag')
 
    created_own_fig = ax is None
    if created_own_fig:
        fig, axes = plt.subplots(1, 4, figsize=figsize, sharey=sharey)
    else:
        fig = None
        axes = ax
 
    all_means = [ poll_means, plant_means,prop1_means, prop2_means]
    all_stds = [ poll_stds, plant_stds,prop1_stds, prop2_stds]
 
    x = np.arange(len(lags))
 
    for panel_ax, means, stds, title in zip(axes, all_means, all_stds, titles):
        for i in range(len(lags)):
            panel_ax.errorbar(x[i], means[i], yerr=stds[i], fmt=marker,
                               color=colors[i], markersize=markersize, capsize=capsize,
                               elinewidth=elinewidth)
        panel_ax.set_title(title, fontsize=fontsize)
        panel_ax.set_xticks(x)
        panel_ax.set_xticklabels(lag_labels, rotation=rotate_labels, fontsize=fontsize)
        panel_ax.tick_params(axis='y', labelsize=fontsize)
        panel_ax.set_xlim(-2,5)
        panel_ax.yaxis.set_major_locator(MultipleLocator(0.2))
        panel_ax.set_xticklabels(['L1', 'L2', 'L3', 'L4'])
        panel_ax.plot([0, 1], [1.12, 1.12], transform=panel_ax.transAxes,color='black', lw=1, clip_on=False)
        panel_ax.plot([0, 0], [1.0, 1.12], transform=panel_ax.transAxes,color='black', lw=1, clip_on=False)
        panel_ax.plot([1, 1], [1.0, 1.12], transform=panel_ax.transAxes,color='black', lw=1, clip_on=False)
        panel_ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
 
    if sharey and created_own_fig:
        for panel_ax in axes[1:]:
            panel_ax.tick_params(axis='y', labelleft=False)
    
    for i, panel_ax in enumerate(axes):
        #panel_ax.spines['top'].set_visible(False)
        panel_ax.spines['right'].set_visible(False)
        if i > 0:
            panel_ax.spines['left'].set_visible(False)
            panel_ax.tick_params(axis='y', left=False)  # tambien quita las marcas de tick, no solo la linea
 
    handles = [mpatches.Patch(color=colors[i], label=f"T{i+1}: {lag}") for i, lag in enumerate(lag_labels)]
    fig.legend(handles=handles, title='Time lag', loc='lower center', bbox_to_anchor=(0.5, -0.2),ncol=len(lag_labels), frameon=False, fontsize=fontsize-2, title_fontsize=fontsize)

    if created_own_fig:
        plt.tight_layout()
        return fig, axes, summary_df
    else:
        return axes, summary_df
 

#el horizontal lo usamos norlamente solo con un sitio 
def plot_lag_summary_combined_horizontal(
        All_Nodes, All_Changes, lags=(1, 2, 3, 4), figsize=(22, 4.5),
        site_level='Site_ID', metric_col='metric', year_col='year',
        value_col_nodes='n_nodes',
        value_col_changes_1='n_interactions', value_col_changes_2='freq_interactions',
        permanent_int_categories=('CN_CI', 'CN_DI', 'CN_II'),
        change_categories=('CN_LI', 'CN_GI', 'L1N_LI', 'LBN_LI', 'G1N_GI', 'GBN_GI'),
        poll_permanent='PA', poll_transient=('A1', 'A2'),
        plant_permanent='PP', plant_transient=('P1', 'P2'),
        offset=0.12, cmap=None, cmap_range=(1.0, 0.4),
        fontsize=10, capsize=0, marker='.', markersize=8,
        rotate_labels=80, label_dict=None, legend_loc='upper right',
        site=None, sharey=True):
    """
    Version HORIZONTAL con 4 PANELES en fila, compartiendo el eje y:
      1. Prop. pollinators permanent
      2. Prop. plants permanent
      3. Prop. interactions permanent -- usando value_col_changes_1
         (por defecto 'n_interactions')
      4. Prop. interactions permanent -- usando value_col_changes_2
         (por defecto 'freq_interactions')
 
    Los paneles 3 y 4 usan la MISMA logica, solo cambia la columna de
    Changes_df -- asi comparas directamente contar vs pesar por frecuencia.
 
    Devuelve
    -------
    fig, (ax1, ax2, ax3, ax4), extra_data
    """
    if cmap is None:
        cmap = plt.cm.Greens
    if label_dict is None:
        label_dict = {
            'poll': 'Prop. pollinators\npermanent',
            'plant': 'Prop. plants\npermanent',
            'prop1': f'Prop. interactions\npermanent ({value_col_changes_1})',
            'prop2': f'Prop. interactions\npermanent ({value_col_changes_2})',
        }
 
    colors = cmap(np.linspace(cmap_range[0], cmap_range[1], len(lags)))
 
    if site is not None:
        All_Nodes_use = {}
        All_Changes_use = {}
        for lag in lags:
            nodes_flat = All_Nodes[lag]
            nodes_flat = nodes_flat.reset_index() if nodes_flat.index.nlevels > 1 or nodes_flat.index.name is not None else nodes_flat.copy()
            changes_flat = All_Changes[lag]
            changes_flat = changes_flat.reset_index() if changes_flat.index.nlevels > 1 or changes_flat.index.name is not None else changes_flat.copy()
 
            All_Nodes_use[lag] = nodes_flat[nodes_flat[site_level] == site]
            All_Changes_use[lag] = changes_flat[changes_flat[site_level] == site]
    else:
        All_Nodes_use = All_Nodes
        All_Changes_use = All_Changes
 
    nodes_flat0 = All_Nodes_use[lags[0]]
    nodes_flat0 = nodes_flat0.reset_index() if nodes_flat0.index.nlevels > 1 or nodes_flat0.index.name is not None else nodes_flat0
    sites = sorted(nodes_flat0[site_level].unique())
    x = np.arange(len(sites))
 
    n = len(lags)
    offsets = [(i - (n - 1) / 2) * offset for i in range(n)]
 
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=figsize, sharey=sharey)
 
    extra_data = {}
    for i, lag in enumerate(lags):
        poll_ratios, plant_ratios, _, prop_1 = _compute_site_year_metrics(
            All_Nodes_use[lag], All_Changes_use[lag], site_level, metric_col, year_col,
            value_col_nodes, value_col_changes_1,
            permanent_int_categories, change_categories,
            poll_permanent, poll_transient, plant_permanent, plant_transient)
        #plant_ratios = {s: (1 - np.array(vals)).tolist() for s, vals in plant_ratios.items()}
        _, _, _, prop_2 = _compute_site_year_metrics(
            All_Nodes_use[lag], All_Changes_use[lag], site_level, metric_col, year_col,
            value_col_nodes, value_col_changes_2,
            permanent_int_categories, change_categories,
            poll_permanent, poll_transient, plant_permanent, plant_transient)
        
        extra_data[lag] = {'prop_1': prop_1, 'prop_2': prop_2}
 
        x_off = x + offsets[i]
        color = colors[i]
 
        for ax, data_dict in zip((ax1, ax2, ax3, ax4), (poll_ratios, plant_ratios, prop_1, prop_2)):
            means = [np.nanmean(data_dict.get(s, [np.nan])) for s in sites]
            stds = [np.nanstd(data_dict.get(s, [np.nan])) for s in sites]
            ax.errorbar(x_off, means, yerr=stds, fmt=marker, capsize=capsize,
                        color=color, markersize=markersize, elinewidth=1.2)
 
    ax1.set_title(label_dict['poll'], fontsize=fontsize)
    ax2.set_title(label_dict['plant'], fontsize=fontsize)
    ax3.set_title(label_dict['prop1'], fontsize=fontsize)
    ax4.set_title(label_dict['prop2'], fontsize=fontsize)
 
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xticks(x)
        ax.set_xticklabels(sites, rotation=rotate_labels, ha='right', fontsize=fontsize)
        ax.set_xlabel(site_level, fontsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
        ax.set_xlim(-0.5, len(sites) - 0.5)
 
    if sharey:
        for ax in (ax2, ax3, ax4):
            ax.tick_params(axis='y', labelleft=False)
 
    handles = [mpatches.Patch(color=colors[i], label=f't+{lags[i]}') for i in range(n)]
    ax4.legend(handles=handles, title='Time lag', loc=legend_loc, frameon=False, fontsize=fontsize)
 
    plt.tight_layout()
 
    return fig, (ax1, ax2, ax3, ax4), extra_data 

# All sites 
def plot_lag_summary_combined_vertical(
        All_Nodes, All_Changes, lags=(1, 2, 3, 4), figsize=(10, 15),
        site_level='Site_ID', metric_col='metric', year_col='year',
        value_col_nodes='n_nodes',
        value_col_changes_1='n_interactions', value_col_changes_2='freq_interactions',
        permanent_int_categories=('CN_CI', 'CN_DI', 'CN_II'),
        change_categories=('CN_LI', 'CN_GI', 'L1N_LI', 'LBN_LI', 'G1N_GI', 'GBN_GI'),
        poll_permanent='PA', poll_transient=('A1', 'A2'),
        plant_permanent='PP', plant_transient=('P1', 'P2'),
        offset=0.12, cmap=None, cmap_range=(1.0, 0.4),
        fontsize=10, capsize=0, marker='.', markersize=8,
        rotate_labels=80, label_dict=None, legend_loc='upper right',
        site=None, sharey=True):
    """
    Version HORIZONTAL con 4 PANELES en fila, compartiendo el eje y:
      1. Prop. pollinators permanent
      2. Prop. plants permanent
      3. Prop. interactions permanent -- usando value_col_changes_1
         (por defecto 'n_interactions')
      4. Prop. interactions permanent -- usando value_col_changes_2
         (por defecto 'freq_interactions')
 
    Los paneles 3 y 4 usan la MISMA logica, solo cambia la columna de
    Changes_df -- asi comparas directamente contar vs pesar por frecuencia.
 
    Devuelve
    -------
    fig, (ax1, ax2, ax3, ax4), extra_data
    """
    if cmap is None:
        cmap = plt.cm.Greens
    if label_dict is None:
        label_dict = {
            'poll': 'Prop. pollinators\npermanent',
            'plant': 'Prop. plants\npermanent',
            'prop1': f'Prop. interactions\npermanent ({value_col_changes_1})',
            'prop2': f'Prop. interactions\npermanent ({value_col_changes_2})',
        }
 
    colors = cmap(np.linspace(cmap_range[0], cmap_range[1], len(lags)))
 
    if site is not None:
        All_Nodes_use = {}
        All_Changes_use = {}
        for lag in lags:
            nodes_flat = All_Nodes[lag]
            nodes_flat = nodes_flat.reset_index() if nodes_flat.index.nlevels > 1 or nodes_flat.index.name is not None else nodes_flat.copy()
            changes_flat = All_Changes[lag]
            changes_flat = changes_flat.reset_index() if changes_flat.index.nlevels > 1 or changes_flat.index.name is not None else changes_flat.copy()
 
            All_Nodes_use[lag] = nodes_flat[nodes_flat[site_level] == site]
            All_Changes_use[lag] = changes_flat[changes_flat[site_level] == site]
    else:
        All_Nodes_use = All_Nodes
        All_Changes_use = All_Changes
 
    nodes_flat0 = All_Nodes_use[lags[0]]
    nodes_flat0 = nodes_flat0.reset_index() if nodes_flat0.index.nlevels > 1 or nodes_flat0.index.name is not None else nodes_flat0
    sites = sorted(nodes_flat0[site_level].unique())
    x = np.arange(len(sites))
 
    n = len(lags)
    offsets = [(i - (n - 1) / 2) * offset for i in range(n)]
 
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=figsize, sharey=sharey, sharex=True)
 
    extra_data = {}
    for i, lag in enumerate(lags):
        poll_ratios, plant_ratios, _, prop_1 = _compute_site_year_metrics(
            All_Nodes_use[lag], All_Changes_use[lag], site_level, metric_col, year_col,
            value_col_nodes, value_col_changes_1,
            permanent_int_categories, change_categories,
            poll_permanent, poll_transient, plant_permanent, plant_transient)
        #plant_ratios = {s: (1 - np.array(vals)).tolist() for s, vals in plant_ratios.items()}
        _, _, _, prop_2 = _compute_site_year_metrics(
            All_Nodes_use[lag], All_Changes_use[lag], site_level, metric_col, year_col,
            value_col_nodes, value_col_changes_2,
            permanent_int_categories, change_categories,
            poll_permanent, poll_transient, plant_permanent, plant_transient)
        
        extra_data[lag] = {'prop_1': prop_1, 'prop_2': prop_2}
 
        x_off = x + offsets[i]
        color = colors[i]
 
        for ax, data_dict in zip((ax1, ax2, ax3, ax4), (poll_ratios, plant_ratios, prop_1, prop_2)):
            means = [np.nanmean(data_dict.get(s, [np.nan])) for s in sites]
            stds = [np.nanstd(data_dict.get(s, [np.nan])) for s in sites]
            ax.errorbar(x_off, means, yerr=stds, fmt=marker, capsize=capsize,
                        color=color, markersize=markersize, elinewidth=1.2)
 
    ax1.set_ylabel(label_dict['poll'], fontsize=fontsize)
    ax2.set_ylabel(label_dict['plant'], fontsize=fontsize)
    ax3.set_ylabel(label_dict['prop1'], fontsize=fontsize)
    ax4.set_ylabel(label_dict['prop2'], fontsize=fontsize)
 
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xticks(x)
        ax.set_xticklabels(sites, rotation=rotate_labels, ha='right', fontsize=fontsize)
        #ax.set_xlabel(site_level, fontsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
        
 
    if sharey:
        for ax in (ax2, ax3, ax4):
            ax.tick_params(axis='y', labelleft=False)
 
    handles = [mpatches.Patch(color=colors[i], label=f't+{lags[i]}') for i in range(n)]
    labels = [h.get_label() for h in handles]
    #ax4.legend(handles=handles, title='Time lag', loc=legend_loc, frameon=False, fontsize=fontsize)
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5),ncol=1, frameon=False, fontsize=fontsize + 1, title="Time lag")
 
    plt.tight_layout()
 
    return fig, (ax1, ax2, ax3, ax4), extra_data
 

# Spatial decay 
from matplotlib.ticker import PercentFormatter
def plot_persistence_by_distance(
        Nodes_df, Changes_df, ax=None, figsize=(20, 4.5),
        bin_level='bin', cmap=None, cmap_range=(1.0, 0.4),
        marker='.', markersize=10, capsize=4, elinewidth=1.2,
        fontsize=20, sharey=True, rotate_labels=0,
        titles=('Permanent floral visitors', 'Permanent flowering plants',
                'Permanent interactions', 'Permanent frequency')):
    """
    Analisis de persistencia EN FUNCION DE LA DISTANCIA (bins de distancia
    entre comparaciones) -- NO es el analisis por lag temporal (analisis
    distinto), aunque el ESTILO visual es el mismo: 4 paneles en fila, eje
    y compartido, punto+barra de error por bin, coloreados con degradado.
 
    El calculo (persistencia = permanentes / total) es EXACTAMENTE el que
    ya tenias -- no se modifica ninguna formula.
 
    Parameters
    ----------
    Nodes_df : pd.DataFrame
        Con columna 'n_nodes' y MultiIndex que incluye 'metric' y bin_level.
    Changes_df : pd.DataFrame
        Con columnas 'n_interactions'/'freq_interactions' y MultiIndex que
        incluye 'metric' y bin_level.
    bin_level : str
        Nivel del indice a usar para agrupar (por defecto 'bin').
    cmap, cmap_range : colormap y rango
        Degradado de color a traves de los bins.
    ax : lista de 4 matplotlib Axes o None
        Si se pasa, dibuja sobre esos 4 ejes. Si no, crea figura nueva.
 
    Devuelve
    -------
    fig (o None si se pasaron ejes), axes (lista de 4), summary_df
    """
    nodes = Nodes_df["n_nodes"].unstack("metric").astype(float)
 
    plant_persistence = (
        nodes["PP"] /
        (nodes["PP"] + nodes["P1"] + nodes["P2"])
    )
 
    pollinator_persistence = (
        nodes["PA"] /
        (nodes["PA"] + nodes["A1"] + nodes["A2"])
    )
 
    changes = Changes_df[
        ["n_interactions", "freq_interactions"]
    ].unstack("metric")
 
    interaction_persistence_n = (
        changes["n_interactions"][["CN_CI", "CN_II", "CN_DI"]].sum(axis=1)
        /
        changes["n_interactions"][[
            "CN_CI", "CN_II", "CN_DI",
            "CN_GI", "CN_LI",
            "L1N_LI", "G1N_GI",
            "LBN_LI", "GBN_GI"
        ]].sum(axis=1)
    )
 
    interaction_persistence_freq = (
        changes["freq_interactions"][["CN_CI", "CN_II", "CN_DI"]].sum(axis=1)
        /
        changes["freq_interactions"][[
            "CN_CI", "CN_II", "CN_DI",
            "CN_GI", "CN_LI",
            "L1N_LI", "G1N_GI",
            "LBN_LI", "GBN_GI"
        ]].sum(axis=1)
    )
 
    plant_mean = plant_persistence.groupby(level=bin_level).mean()
    plant_std = plant_persistence.groupby(level=bin_level).std()
 
    pollinator_mean = pollinator_persistence.groupby(level=bin_level).mean()
    pollinator_std = pollinator_persistence.groupby(level=bin_level).std()
 
    interaction_n_mean = interaction_persistence_n.groupby(level=bin_level).mean()
    interaction_n_std = interaction_persistence_n.groupby(level=bin_level).std()
 
    interaction_freq_mean = interaction_persistence_freq.groupby(level=bin_level).mean()
    interaction_freq_std = interaction_persistence_freq.groupby(level=bin_level).std()
 
    bins = list(plant_mean.index)
 
    summary_df = pd.DataFrame({
        'plant_mean': plant_mean, 'plant_std': plant_std,
        'pollinator_mean': pollinator_mean, 'pollinator_std': pollinator_std,
        'interaction_n_mean': interaction_n_mean, 'interaction_n_std': interaction_n_std,
        'interaction_freq_mean': interaction_freq_mean, 'interaction_freq_std': interaction_freq_std,
    })
 
    if cmap is None:
        cmap = plt.cm.Greens
    colors = cmap(np.linspace(cmap_range[0], cmap_range[1], len(bins)))
 
    created_own_fig = ax is None
    if created_own_fig:
        fig, axes = plt.subplots(1, 4, figsize=figsize, sharey=sharey)
    else:
        fig = None
        axes = ax
 
    all_means = [plant_mean, pollinator_mean, interaction_n_mean, interaction_freq_mean]
    all_stds = [plant_std, pollinator_std, interaction_n_std, interaction_freq_std]
 
    x = np.arange(len(bins))
 
    for panel_ax, means, stds, title in zip(axes, all_means, all_stds, titles):
        for i, b in enumerate(bins):
            panel_ax.errorbar(x[i], means.loc[b], yerr=stds.loc[b], fmt=marker,
                               color=colors[i], markersize=markersize, capsize=capsize,
                               elinewidth=elinewidth)
        panel_ax.set_title(title, fontsize=fontsize)
        panel_ax.set_xticks(x)
        panel_ax.set_xticklabels(bins, rotation=rotate_labels, fontsize=fontsize)
        panel_ax.tick_params(axis='y', labelsize=fontsize)
        panel_ax.set_xlim(-2,5)
        panel_ax.yaxis.set_major_locator(MultipleLocator(0.2))
        panel_ax.set_xticklabels(['D1', 'D2', 'D3', 'D4'])
        panel_ax.plot([0, 1], [1.1, 1.1], transform=panel_ax.transAxes,color='black', lw=1, clip_on=False)
        panel_ax.plot([0, 0], [1.0, 1.1], transform=panel_ax.transAxes,color='black', lw=1, clip_on=False)
        panel_ax.plot([1, 1], [1.0, 1.1], transform=panel_ax.transAxes,color='black', lw=1, clip_on=False)
        panel_ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
 
    if sharey and created_own_fig:
        for panel_ax in axes[1:]:
            panel_ax.tick_params(axis='y', labelleft=False)
            
    for i, panel_ax in enumerate(axes):
        #panel_ax.spines['top'].set_visible(False)
        panel_ax.spines['right'].set_visible(False)
        if i > 0:
            panel_ax.spines['left'].set_visible(False)
            panel_ax.tick_params(axis='y', left=False)  # tambien quita las marcas de tick, no solo la linea
    
    handles = [mpatches.Patch(color=colors[i], label=f"D{i+1}: {b} km") for i, b in enumerate(bins)]
    fig.legend(handles=handles, title='Spatial distance', loc='lower center',bbox_to_anchor=(0.5, -0.2), ncol=len(bins), frameon=False, fontsize=fontsize-2,title_fontsize=fontsize)

    if created_own_fig:
        plt.tight_layout()
        return fig, axes, summary_df
    else:
        return axes, summary_df


#otros

def test_interaction_frequency_differences2(
        df_site, freq_col='freq_interactions', n_col='n_interactions',
        metrics=('CN_CI', 'CN_DI', 'CN_II', 'CN_LI', 'CN_GI', 'L1N_LI', 'LBN_LI', 'G1N_GI', 'GBN_GI'),
        metric_col='metric', alpha=0.05, p_adjust='holm',
        label_dict=None, generate_latex=True, latex_caption='Post-hoc Dunn test (p-values)',
        latex_label='tab:dunn_posthoc'):
    """
    Test estadistico para saber si la frecuencia media por interaccion
    (freq_col/n_col) difiere entre las categorias.
 
    1. Kruskal-Wallis: compara todas las categorias A LA VEZ -- responde
       "¿hay AL MENOS una categoria que difiere de las demas?".
    2. Post-hoc de Dunn (con correccion por comparaciones multiples): SOLO
       si el Kruskal-Wallis sale significativo, compara TODOS los pares.
 
    NUEVO: si generate_latex=True, ademas devuelve el string de una tabla
    LaTeX (formato booktabs) con la matriz de p-valores de Dunn, lista para
    pegar en tu documento.
 
    Parameters
    ----------
    df_site : pd.DataFrame
        Ya filtrado a UN Site_ID (o varios, si quieres agregar mas puntos).
    metrics : tuple de str
        Categorias a comparar.
    alpha : float
        Umbral de significacion.
    p_adjust : str
        Metodo de correccion para el post-hoc de Dunn -- 'holm' (por
        defecto) o 'bonferroni'.
    label_dict : dict o None
        {codigo_metrica: texto legible} -- si se pasa, se usa para las
        filas/columnas de la tabla LaTeX en vez de los codigos crudos.
    generate_latex : bool
        Si True (por defecto), genera tambien la tabla LaTeX.
    latex_caption, latex_label : str
        Caption y label (\\label{}) de la tabla LaTeX.
 
    Devuelve
    -------
    kw_result : dict
        {'H': estadistico, 'p_value': float, 'significant': bool}
    posthoc_df : pd.DataFrame o None
        Matriz de p-valores (Dunn) -- None si el Kruskal-Wallis no fue
        significativo.
    latex_table : str o None
        La tabla en LaTeX (string) -- None si generate_latex=False o si el
        Kruskal-Wallis no fue significativo (no hay tabla que generar).
    """
    df_flat = df_site.reset_index() if df_site.index.nlevels > 1 or df_site.index.name is not None else df_site.copy()
    if label_dict is None:
        label_dict = {}
 
    groups = []
    labels = []
    for metric in metrics:
        sub = df_flat.loc[df_flat[metric_col] == metric]
        vals = (sub[freq_col] / sub[n_col]).replace([np.inf, -np.inf], np.nan).dropna()
        if len(vals) > 0:
            groups.append(vals.values)
            labels.append(metric)
 
    # ------------------------------------------------------------------
    # 1) Kruskal-Wallis
    # ------------------------------------------------------------------
    H, p_value = kruskal(*groups)
    significant = p_value < alpha
    kw_result = {'H': H, 'p_value': p_value, 'significant': significant}
 
    print(f"Kruskal-Wallis: H={H:.3f}, p={p_value:.4f} "
          f"({'SIGNIFICATIVO' if significant else 'no significativo'} a alpha={alpha})")
 
    # ------------------------------------------------------------------
    # 2) Post-hoc de Dunn -- SOLO si el paso 1 fue significativo
    # ------------------------------------------------------------------
    posthoc_df = None
    latex_table = None
 
    if significant:
        long_data = []
        for label, vals in zip(labels, groups):
            for v in vals:
                long_data.append({'metric': label, 'value': v})
        long_df = pd.DataFrame(long_data)
 
        posthoc_df = sp.posthoc_dunn(long_df, val_col='value', group_col='metric', p_adjust=p_adjust)
        posthoc_df = posthoc_df.loc[labels, labels]
 
        print(f"\nPost-hoc de Dunn (p-valores, correccion {p_adjust}):")
        print(posthoc_df.round(4))
 
        if generate_latex:
            latex_table = _posthoc_to_latex(posthoc_df, label_dict, latex_caption, latex_label, alpha)
    else:
        print("Como el Kruskal-Wallis no fue significativo, no se hace el post-hoc.")
 
    return kw_result, posthoc_df, latex_table
 
 
def _posthoc_to_latex(posthoc_df, label_dict, caption, label, alpha):
    """
    Convierte la matriz de p-valores de Dunn en una tabla LaTeX (booktabs),
    con los p-valores significativos marcados en negrita y con asteriscos.
    """
    def format_cell(p, is_diagonal):
        if is_diagonal:
            return '--'
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < alpha else ''
        formatted = f"{p:.3f}{stars}"
        if p < alpha:
            formatted = f"\\textbf{{{formatted}}}"
        return formatted
 
    display_labels = [label_dict.get(m, m).replace('_', '\\_') for m in posthoc_df.index]
 
    lines = []
    lines.append("\\begin{table}[ht]")
    lines.append("\\centering")
    lines.append("\\begin{tabular}{l" + "c" * len(posthoc_df.columns) + "}")
    lines.append("\\toprule")
    lines.append(" & " + " & ".join(display_labels) + " \\\\")
    lines.append("\\midrule")
 
    for i, row_metric in enumerate(posthoc_df.index):
        row_cells = []
        for j, col_metric in enumerate(posthoc_df.columns):
            p = posthoc_df.loc[row_metric, col_metric]
            row_cells.append(format_cell(p, i == j))
        lines.append(display_labels[i] + " & " + " & ".join(row_cells) + " \\\\")
 
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append(f"\\caption{{{caption}. Bold values indicate significant differences "
                  f"(* p$<${alpha}, ** p$<$0.01, *** p$<$0.001).}}")
    lines.append(f"\\label{{{label}}}")
    lines.append("\\end{table}")
 
    latex_str = "\n".join(lines)
    print("\n" + "=" * 60)
    print("Tabla LaTeX:")
    print("=" * 60)
    print(latex_str)
 
    return latex_str
 


