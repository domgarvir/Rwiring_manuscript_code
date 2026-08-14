from Functions import *

######### REWIRING FUNCTIONS
##Helper to quantify changes
def distance(x_t1, x_t2, method="euclidean"):
    # Combine x_t1 and x_t2 into a DataFrame
    x = pd.concat([pd.DataFrame([x_t1]), pd.DataFrame([x_t2])])
    
    # Calculate the pairwise distance
    y = pdist(x, method)
    
    # Return the first distance value
    return y[0]

#quantifying changes in interactions from year to year
def rewiring_py_new(A, B, quantitative=False, method="euclidean"):
    A_bin = (A != 0).astype(int)
    B_bin = (B != 0).astype(int)

    if not quantitative:
        A = A_bin
        B = B_bin

    common_col = A.columns.intersection(B.columns)
    common_row = A.index.intersection(B.index)

    # --- NEW: species exclusive to each network (needed for LBN / GBN) ---
    rows_only_A = A.index.difference(common_row)     # plants lost from t1 -> t2
    cols_only_A = A.columns.difference(common_col)   # pollinators lost
    rows_only_B = B.index.difference(common_row)     # plants gained
    cols_only_B = B.columns.difference(common_col)   # pollinators gained

    CN1 = A.loc[common_row, common_col]
    CN2 = B.loc[common_row, common_col]

    CN1_bin = A_bin.loc[common_row, common_col]
    CN2_bin = B_bin.loc[common_row, common_col]

    index_CN_SI = np.argwhere((CN1.values > 0) & (CN2.values > 0))
    index_CN_SNI = np.argwhere((CN1.values == 0) & (CN1.values == CN2.values))
    index_CN_LI = np.argwhere((CN1.values > 0) & (CN2.values == 0))
    index_CN_GI = np.argwhere((CN1.values == 0) & (CN2.values > 0))

    CN_SI = len(index_CN_SI)
    CN_SNI = len(index_CN_SNI)
    CN_LI = len(index_CN_LI)
    CN_GI = len(index_CN_GI)

    if CN_SI > 0:
        d_q_CN_SI = distance(CN1.values[index_CN_SI[:, 0], index_CN_SI[:, 1]], CN2.values[index_CN_SI[:, 0], index_CN_SI[:, 1]], method=method)
        d_B_CN_SI = distance(CN1_bin.values[index_CN_SI[:, 0], index_CN_SI[:, 1]], CN2_bin.values[index_CN_SI[:, 0], index_CN_SI[:, 1]], method=method)
        W_CN_SI = (1 + d_q_CN_SI) / (1 + d_B_CN_SI)
    else:
        W_CN_SI = 1

    if CN_SNI > 0:
        d_q_CN_SNI = distance(CN1.values[index_CN_SNI[:, 0], index_CN_SNI[:, 1]], CN2.values[index_CN_SNI[:, 0], index_CN_SNI[:, 1]], method=method)
        d_B_CN_SNI = distance(CN1.values[index_CN_SNI[:, 0], index_CN_SNI[:, 1]], CN2.values[index_CN_SNI[:, 0], index_CN_SNI[:, 1]], method=method)
        W_CN_SNI = (1 + d_q_CN_SNI) / (1 + d_B_CN_SNI)
    else:
        W_CN_SNI = 1

    if CN_LI > 0:
        d_q_CN_LI = distance(CN1.values[index_CN_LI[:, 0], index_CN_LI[:, 1]], CN2.values[index_CN_LI[:, 0], index_CN_LI[:, 1]], method=method)
        d_B_CN_LI = distance(CN1_bin.values[index_CN_LI[:, 0], index_CN_LI[:, 1]], CN2_bin.values[index_CN_LI[:, 0], index_CN_LI[:, 1]], method=method)
        W_CN_LI = (1 + d_q_CN_LI) / (1 + d_B_CN_LI)
    else:
        W_CN_LI = 1

    if CN_GI > 0:
        d_q_CN_GI = distance(CN1.values[index_CN_GI[:, 0], index_CN_GI[:, 1]], CN2.values[index_CN_GI[:, 0], index_CN_GI[:, 1]], method=method)
        d_B_CN_GI = distance(CN1_bin.values[index_CN_GI[:, 0], index_CN_GI[:, 1]], CN2_bin.values[index_CN_GI[:, 0], index_CN_GI[:, 1]], method=method)
        W_CN_GI = (1 + d_q_CN_GI) / (1 + d_B_CN_GI)
    else:
        W_CN_GI = 1

    # --- LN: interactions lost where only ONE of the two species leaves ---
    LN = A.copy()
    LN.loc[common_row, common_col] = np.nan          # hide core cells (CN_*)
    LN.loc[rows_only_A, cols_only_A] = np.nan        # NEW: hide "both species lost" cells (-> LBN)

    index_L1N_LI = np.argwhere((LN.values > 0) & ~np.isnan(LN.values))
    index_LN_SNI = np.argwhere((LN.values == 0) & ~np.isnan(LN.values))

    L1N_LI = len(index_L1N_LI)   # lose ONE node -> lost interactions
    LN_SNI = len(index_LN_SNI)

    LN_bin = LN.copy()
    LN_bin[(LN_bin != 0) & ~LN_bin.isna()] = 1

    W_LN_SNI = 1

    if L1N_LI > 0:
        d_q_L1N_LI = distance(LN.values[index_L1N_LI[:, 0], index_L1N_LI[:, 1]], 0 * LN.values[index_L1N_LI[:, 0], index_L1N_LI[:, 1]], method=method)
        d_B_L1N_LI = distance(LN_bin.values[index_L1N_LI[:, 0], index_L1N_LI[:, 1]], 0 * LN_bin.values[index_L1N_LI[:, 0], index_L1N_LI[:, 1]], method=method)
        W_L1N_LI = (1 + d_q_L1N_LI) / (1 + d_B_L1N_LI)
    else:
        W_L1N_LI = 1

    # --- NEW: LBN: interactions lost because BOTH species leave ---
    LBN = A.loc[rows_only_A, cols_only_A]
    LBN_bin = A_bin.loc[rows_only_A, cols_only_A]

    index_LBN_LI = np.argwhere(LBN.values > 0)
    index_LBN_SNI = np.argwhere(LBN.values == 0)

    LBN_LI = len(index_LBN_LI)
    LBN_SNI = len(index_LBN_SNI)

    W_LBN_SNI = 1

    if LBN_LI > 0:
        d_q_LBN_LI = distance(LBN.values[index_LBN_LI[:, 0], index_LBN_LI[:, 1]], 0 * LBN.values[index_LBN_LI[:, 0], index_LBN_LI[:, 1]], method=method)
        d_B_LBN_LI = distance(LBN_bin.values[index_LBN_LI[:, 0], index_LBN_LI[:, 1]], 0 * LBN_bin.values[index_LBN_LI[:, 0], index_LBN_LI[:, 1]], method=method)
        W_LBN_LI = (1 + d_q_LBN_LI) / (1 + d_B_LBN_LI)
    else:
        W_LBN_LI = 1

    # --- GN: interactions gained where only ONE of the two species is new ---
    GN = B.copy()
    GN.loc[common_row, common_col] = np.nan          # hide core cells (CN_*)
    GN.loc[rows_only_B, cols_only_B] = np.nan        # NEW: hide "both species gained" cells (-> GBN)

    index_G1N_GI = np.argwhere((GN.values > 0) & ~np.isnan(GN.values))
    index_GN_SNI = np.argwhere((GN.values == 0) & ~np.isnan(GN.values))

    G1N_GI = len(index_G1N_GI)   # gain ONE node -> gained interactions
    GN_SNI = len(index_GN_SNI)

    GN_bin = GN.copy()
    GN_bin[(GN_bin != 0) & ~GN_bin.isna()] = 1

    if G1N_GI > 0:
        d_q_G1N_GI = distance(0 * GN.values[index_G1N_GI[:, 0], index_G1N_GI[:, 1]], GN.values[index_G1N_GI[:, 0], index_G1N_GI[:, 1]], method=method)
        d_B_G1N_GI = distance(0 * GN_bin.values[index_G1N_GI[:, 0], index_G1N_GI[:, 1]], GN_bin.values[index_G1N_GI[:, 0], index_G1N_GI[:, 1]], method=method)
        W_G1N_GI = (1 + d_q_G1N_GI) / (1 + d_B_G1N_GI)
    else:
        W_G1N_GI = 1

    W_GN_SNI = 1

    # --- NEW: GBN: interactions gained because BOTH species are new ---
    GBN = B.loc[rows_only_B, cols_only_B]
    GBN_bin = B_bin.loc[rows_only_B, cols_only_B]

    index_GBN_GI = np.argwhere(GBN.values > 0)
    index_GBN_SNI = np.argwhere(GBN.values == 0)

    GBN_GI = len(index_GBN_GI)
    GBN_SNI = len(index_GBN_SNI)

    W_GBN_SNI = 1

    if GBN_GI > 0:
        d_q_GBN_GI = distance(0 * GBN.values[index_GBN_GI[:, 0], index_GBN_GI[:, 1]], GBN.values[index_GBN_GI[:, 0], index_GBN_GI[:, 1]], method=method)
        d_B_GBN_GI = distance(0 * GBN_bin.values[index_GBN_GI[:, 0], index_GBN_GI[:, 1]], GBN_bin.values[index_GBN_GI[:, 0], index_GBN_GI[:, 1]], method=method)
        W_GBN_GI = (1 + d_q_GBN_GI) / (1 + d_B_GBN_GI)
    else:
        W_GBN_GI = 1

    # --- atomic categories: every network cell falls into exactly ONE of these ---
    atomic_metrics = ["CN_LI", "L1N_LI", "LBN_LI", "CN_GI", "G1N_GI", "GBN_GI",
                       "CN_SI", "CN_SNI", "LN_SNI", "GN_SNI"]
    atomic_values = [CN_LI, L1N_LI, LBN_LI, CN_GI, G1N_GI, GBN_GI,
                      CN_SI, CN_SNI, LN_SNI, GN_SNI]
    atomic_weights = [W_CN_LI, W_L1N_LI, W_LBN_LI, W_CN_GI, W_G1N_GI, W_GBN_GI,
                       W_CN_SI, W_CN_SNI, W_LN_SNI, W_GN_SNI]
    atomic_values_weighted = [v * w for v, w in zip(atomic_values, atomic_weights)]

    # total_links / total_weighted are computed from atomic categories ONLY,
    # so every cell of A and B is counted exactly once
    total_links = sum(atomic_values)
    total_weighted = sum(atomic_values_weighted)

    # --- derived totals: LN_LI = L1N_LI + LBN_LI ; GN_GI = G1N_GI + GBN_GI ---
    # these are convenience rows, NOT included in total_links/total_weighted
    # (their cells are already counted in the atomic rows above)
    LN_LI = L1N_LI + LBN_LI
    LN_LI_weighted = (L1N_LI * W_L1N_LI) + (LBN_LI * W_LBN_LI)
    W_LN_LI = LN_LI_weighted / LN_LI if LN_LI > 0 else 1

    GN_GI = G1N_GI + GBN_GI
    GN_GI_weighted = (G1N_GI * W_G1N_GI) + (GBN_GI * W_GBN_GI)
    W_GN_GI = GN_GI_weighted / GN_GI if GN_GI > 0 else 1

    derived_metrics = ["LN_LI", "GN_GI"]
    derived_values = [LN_LI, GN_GI]
    derived_weights = [W_LN_LI, W_GN_GI]
    derived_values_weighted = [LN_LI_weighted, GN_GI_weighted]

    metrics = atomic_metrics + derived_metrics
    values = atomic_values + derived_values
    weights = atomic_weights + derived_weights
    values_weighted = atomic_values_weighted + derived_values_weighted
    is_total = [False] * len(atomic_metrics) + [True] * len(derived_metrics)

    result = pd.DataFrame({
        'metric': metrics,
        'value': values,
        'weight': weights,
        'value_weighted': values_weighted,
        'is_derived_total': is_total,
        'total_links': [total_links] * len(metrics),
        'total_weighted': [total_weighted] * len(metrics)
    })

    return result

def rewiring_freq_py(A, B):
    """
    Cuantifica el cambio entre dos redes (A=t1, B=t2) usando la SUMA de frecuencias
    involucradas en cada categoria, en vez de un ratio de distancia cuantitativa/binaria.

    Para cada categoria, "freq_sum" es:
    - Interacciones que se PIERDEN (CN_LI, L1N_LI, LBN_LI): suma de sus frecuencias en A (t1).
      Ej: se pierden 3 interacciones con frecuencias 2, 3, 1 -> freq_sum = 2+3+1 = 6
    - Interacciones que se GANAN (CN_GI, G1N_GI, GBN_GI): suma de sus frecuencias en B (t2).
      Ej: se ganan 4 interacciones con frecuencias 2, 1, 3, 1 -> freq_sum = 2+1+3+1 = 7
    - Interacciones que PERSISTEN (CN_SI): suma de |freq_t2 - freq_t1| de cada una.
      Ej: t1=[1,2,3,1], t2=[1,3,2,1] -> freq_sum = |1-1|+|3-2|+|2-3|+|1-1| = 0+1+1+0 = 2
      Se desglosa ademas por direccion del cambio:
      - CN_II (increased interactions): cuantas de esas persistentes AUMENTARON de peso,
        y la suma de cuanto aumentaron. En el ejemplo: 1 interaccion (2->3), freq_sum=1.
      - CN_DI (decreased interactions): cuantas DISMINUYERON, y la suma de cuanto.
        En el ejemplo: 1 interaccion (3->2), freq_sum=1.
      - CN_CI (constant interactions): cuantas NO cambiaron de frecuencia (t2==t1), y
        cuanta frecuencia acumulan esas interacciones constantes (suma de su frecuencia,
        que es la misma en t1 y t2). En el ejemplo: 2 interacciones (1->1, 1->1),
        freq_sum = 1+1 = 2.
      Las que no cambian (diff=0) no cuentan en CN_II ni CN_DI, pero siguen dentro de CN_SI
      (y ahora si cuentan, exclusivamente, en CN_CI).
    - Ausencias estables (CN_SNI, LN_SNI, GN_SNI): siempre 0 (no hay frecuencia que sumar).

    Se reporta tambien "n_interactions" (el conteo, igual que en rewiring_py) para que
    puedas comparar cuantas interacciones hay en cada categoria frente a cuanta frecuencia
    representan -- dos interacciones pueden ser el mismo "n" pero muy distinto "freq_sum".

    Parameters
    ----------
    A, B : pd.DataFrame
        Matrices de interaccion (filas=plantas, columnas=polinizadores) en t1 y t2.
        Se usan los valores CRUDOS (frecuencias), sin binarizar salvo para clasificar
        que celda cae en que categoria (presencia/ausencia).

    Returns
    -------
    pd.DataFrame con columnas: metric, n_interactions, freq_sum, is_derived_total
    """
    A_bin = (A != 0).astype(int)
    B_bin = (B != 0).astype(int)

    common_col = A.columns.intersection(B.columns)
    common_row = A.index.intersection(B.index)

    rows_only_A = A.index.difference(common_row)
    cols_only_A = A.columns.difference(common_col)
    rows_only_B = B.index.difference(common_row)
    cols_only_B = B.columns.difference(common_col)

    # ------------------------------------------------------------------
    # bloque comun (ambas especies presentes en A y B)
    # ------------------------------------------------------------------
    CN1 = A.loc[common_row, common_col]
    CN2 = B.loc[common_row, common_col]

    mask_CN_SI = (CN1.values > 0) & (CN2.values > 0)
    mask_CN_SNI = (CN1.values == 0) & (CN2.values == 0)
    mask_CN_LI = (CN1.values > 0) & (CN2.values == 0)
    mask_CN_GI = (CN1.values == 0) & (CN2.values > 0)

    n_CN_SI = int(mask_CN_SI.sum())
    n_CN_SNI = int(mask_CN_SNI.sum())
    n_CN_LI = int(mask_CN_LI.sum())
    n_CN_GI = int(mask_CN_GI.sum())

    freq_CN_SI = float(np.abs(CN2.values[mask_CN_SI] - CN1.values[mask_CN_SI]).sum())
    freq_CN_SNI = 0.0
    freq_CN_LI = float(CN1.values[mask_CN_LI].sum())
    freq_CN_GI = float(CN2.values[mask_CN_GI].sum())

    # --- desglose de CN_SI por direccion del cambio: cuantas AUMENTAN vs DISMINUYEN ---
    # (las que no cambian -diff=0- no cuentan en ninguna de las dos; siguen dentro de CN_SI)
    diff_CN_SI = CN2.values[mask_CN_SI] - CN1.values[mask_CN_SI]
    mask_increase = diff_CN_SI > 0
    mask_decrease = diff_CN_SI < 0

    n_CN_II = int(mask_increase.sum())   # constant nodes, INCREASED interaction strength
    n_CN_DI = int(mask_decrease.sum())   # constant nodes, DECREASED interaction strength
    freq_CN_II = float(diff_CN_SI[mask_increase].sum())
    freq_CN_DI = float(np.abs(diff_CN_SI[mask_decrease]).sum())

    # constant nodes, CONSTANT interactions (freq_t2 == freq_t1, no change at all)
    mask_constant = diff_CN_SI == 0
    n_CN_CI = int(mask_constant.sum())
    freq_CN_CI = float(CN1.values[mask_CN_SI][mask_constant].sum())  # = freq acumulada (t1=t2 aqui)

    # ------------------------------------------------------------------
    # LN: se pierde una interaccion donde solo UNA especie desaparece
    # ------------------------------------------------------------------
    LN = A.copy()
    LN.loc[common_row, common_col] = np.nan
    LN.loc[rows_only_A, cols_only_A] = np.nan

    mask_L1N_LI = (LN.values > 0) & ~np.isnan(LN.values)
    mask_LN_SNI = (LN.values == 0) & ~np.isnan(LN.values)

    n_L1N_LI = int(mask_L1N_LI.sum())
    n_LN_SNI = int(mask_LN_SNI.sum())
    freq_L1N_LI = float(LN.values[mask_L1N_LI].sum())
    freq_LN_SNI = 0.0

    # ------------------------------------------------------------------
    # LBN: se pierde una interaccion donde AMBAS especies desaparecen
    # ------------------------------------------------------------------
    LBN = A.loc[rows_only_A, cols_only_A]
    mask_LBN_LI = LBN.values > 0
    mask_LBN_SNI = LBN.values == 0

    n_LBN_LI = int(mask_LBN_LI.sum())
    n_LBN_SNI = int(mask_LBN_SNI.sum())
    freq_LBN_LI = float(LBN.values[mask_LBN_LI].sum())
    freq_LBN_SNI = 0.0

    # ------------------------------------------------------------------
    # GN: se gana una interaccion donde solo UNA especie es nueva
    # ------------------------------------------------------------------
    GN = B.copy()
    GN.loc[common_row, common_col] = np.nan
    GN.loc[rows_only_B, cols_only_B] = np.nan

    mask_G1N_GI = (GN.values > 0) & ~np.isnan(GN.values)
    mask_GN_SNI = (GN.values == 0) & ~np.isnan(GN.values)

    n_G1N_GI = int(mask_G1N_GI.sum())
    n_GN_SNI = int(mask_GN_SNI.sum())
    freq_G1N_GI = float(GN.values[mask_G1N_GI].sum())
    freq_GN_SNI = 0.0

    # ------------------------------------------------------------------
    # GBN: se gana una interaccion donde AMBAS especies son nuevas
    # ------------------------------------------------------------------
    GBN = B.loc[rows_only_B, cols_only_B]
    mask_GBN_GI = GBN.values > 0
    mask_GBN_SNI = GBN.values == 0

    n_GBN_GI = int(mask_GBN_GI.sum())
    n_GBN_SNI = int(mask_GBN_SNI.sum())
    freq_GBN_GI = float(GBN.values[mask_GBN_GI].sum())
    freq_GBN_SNI = 0.0

    # ------------------------------------------------------------------
    # tabla final: categorias atomicas (sin solapamiento) + totales derivados
    # ------------------------------------------------------------------
    atomic_metrics = ["CN_LI", "L1N_LI", "LBN_LI", "CN_GI", "G1N_GI", "GBN_GI",
                       "CN_SI", "CN_SNI", "LN_SNI", "GN_SNI"]
    atomic_n = [n_CN_LI, n_L1N_LI, n_LBN_LI, n_CN_GI, n_G1N_GI, n_GBN_GI,
                n_CN_SI, n_CN_SNI, n_LN_SNI, n_GN_SNI]
    atomic_freq = [freq_CN_LI, freq_L1N_LI, freq_LBN_LI, freq_CN_GI, freq_G1N_GI, freq_GBN_GI,
                   freq_CN_SI, freq_CN_SNI, freq_LN_SNI, freq_GN_SNI]

    # totales derivados (informativos, NO se suman aparte para no duplicar celdas)
    # LN_LI/GN_GI = descomponen por tipo de nodo; CN_II/CN_DI = descomponen CN_SI por direccion
    derived_metrics = ["LN_LI", "GN_GI", "CN_II", "CN_DI"]
    derived_n = [n_L1N_LI + n_LBN_LI, n_G1N_GI + n_GBN_GI, n_CN_II, n_CN_DI]
    derived_freq = [freq_L1N_LI + freq_LBN_LI, freq_G1N_GI + freq_GBN_GI, freq_CN_II, freq_CN_DI]

    # ------------------------------------------------------------------
    # tabla final: SOLO estas 8 categorias, en este orden (por ahora)
    # ------------------------------------------------------------------
    metrics = [ "CN_CI","CN_DI", "CN_II", "CN_LI", "CN_GI", "L1N_LI", "LBN_LI", "G1N_GI", "GBN_GI"]
    n_interactions = [n_CN_CI, n_CN_DI, n_CN_II,  n_CN_LI, n_CN_GI, n_L1N_LI, n_LBN_LI, n_G1N_GI, n_GBN_GI]
    freq_sum = [freq_CN_CI,freq_CN_DI, freq_CN_II,  freq_CN_LI, freq_CN_GI, freq_L1N_LI, freq_LBN_LI, freq_G1N_GI, freq_GBN_GI]

    result = pd.DataFrame({
        'metric': metrics,
        'n_interactions': n_interactions,
        'freq_sum': freq_sum,
    })

    return result