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

def rewiring_freq_py(A, B):
    """
    Quantify the change between two networks (A=t1, B=t2) using the SUM of the
    frequencies involved in each category, instead of a quantitative/binary
    distance ratio.

    For each category, "freq_sum" is:
    - Interactions that are LOST (CN_LI, L1N_LI, LBN_LI): sum of their frequencies in A (t1).
      E.g., 3 interactions with frequencies 2, 3, 1 are lost -> freq_sum = 2+3+1 = 6
    - Interactions that are GAINED (CN_GI, G1N_GI, GBN_GI): sum of their frequencies in B (t2).
      E.g., 4 interactions with frequencies 2, 1, 3, 1 are gained -> freq_sum = 2+1+3+1 = 7
    - Interactions that PERSIST (CN_SI): sum of |freq_t2 - freq_t1| for each one.
      E.g., t1=[1,2,3,1], t2=[1,3,2,1] -> freq_sum = |1-1|+|3-2|+|2-3|+|1-1| = 0+1+1+0 = 2
      This is further broken down by the direction of the change:
      - CN_II (increased interactions): how many of these persistent interactions
        INCREASED in strength, and the sum of the increases. In the example:
        1 interaction (2->3), freq_sum = 1.
      - CN_DI (decreased interactions): how many DECREASED, and the sum of the decreases.
        In the example: 1 interaction (3->2), freq_sum = 1.
      - CN_CI (constant interactions): how many did NOT change in frequency (t2 == t1),
        and the total frequency these constant interactions add up to (the sum of their
        frequency, which is the same in t1 and t2). In the example: 2 interactions
        (1->1, 1->1), freq_sum = 1+1 = 2.
      Interactions with no change (diff = 0) are not counted in CN_II or CN_DI, but
      remain within CN_SI (and are now counted, exclusively, in CN_CI).
    - Stable absences (CN_SNI, LN_SNI, GN_SNI): always 0 (there is no frequency to sum).

    "n_interactions" (the count, as in rewiring_py) is also reported so that you can
    compare how many interactions each category contains against how much frequency
    they represent -- two categories can have the same "n" but very different "freq_sum".

    Parameters
    ----------
    A, B : pd.DataFrame
        Interaction matrices (rows = plants, columns = pollinators) at t1 and t2.
        RAW values (frequencies) are used; they are binarized only to classify which
        category each cell falls into (presence/absence).

    Returns
    -------
    pd.DataFrame with columns: metric, n_interactions, freq_sum, is_derived_total
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
  # common block (both species present in A and B)
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

    # --- breakdown of CN_SI by direction of change: how many INCREASE vs DECREASE ---
    # (those that do not change -diff=0- are counted in neither; they remain within CN_SI)
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
    freq_CN_CI = float(CN1.values[mask_CN_SI][mask_constant].sum())  # = freq cumulative (t1=t2)

    # ------------------------------------------------------------------
    # LN: an interaction is lost where only ONE species disappears
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
    # LBN: an interaction is lost where BOTH species disappear
    # ------------------------------------------------------------------
    LBN = A.loc[rows_only_A, cols_only_A]
    mask_LBN_LI = LBN.values > 0
    mask_LBN_SNI = LBN.values == 0

    n_LBN_LI = int(mask_LBN_LI.sum())
    n_LBN_SNI = int(mask_LBN_SNI.sum())
    freq_LBN_LI = float(LBN.values[mask_LBN_LI].sum())
    freq_LBN_SNI = 0.0

    # ------------------------------------------------------------------
    # GN: an interaction is gained where only ONE species is new
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
    # GBN: an interaction is gained where BOTH species are new
    # ------------------------------------------------------------------
    GBN = B.loc[rows_only_B, cols_only_B]
    mask_GBN_GI = GBN.values > 0
    mask_GBN_SNI = GBN.values == 0

    n_GBN_GI = int(mask_GBN_GI.sum())
    n_GBN_SNI = int(mask_GBN_SNI.sum())
    freq_GBN_GI = float(GBN.values[mask_GBN_GI].sum())
    freq_GBN_SNI = 0.0

    # ------------------------------------------------------------------
    # final table: atomic categories (no overlap) + derived totals
    # ------------------------------------------------------------------
    atomic_metrics = ["CN_LI", "L1N_LI", "LBN_LI", "CN_GI", "G1N_GI", "GBN_GI",
                       "CN_SI", "CN_SNI", "LN_SNI", "GN_SNI"]
    atomic_n = [n_CN_LI, n_L1N_LI, n_LBN_LI, n_CN_GI, n_G1N_GI, n_GBN_GI,
                n_CN_SI, n_CN_SNI, n_LN_SNI, n_GN_SNI]
    atomic_freq = [freq_CN_LI, freq_L1N_LI, freq_LBN_LI, freq_CN_GI, freq_G1N_GI, freq_GBN_GI,
                   freq_CN_SI, freq_CN_SNI, freq_LN_SNI, freq_GN_SNI]

    # derived totals (informative, NOT summed separately to avoid double-counting cells)
    # LN_LI/GN_GI = broken down by node type; CN_II/CN_DI = break down CN_SI by direction
    derived_metrics = ["LN_LI", "GN_GI", "CN_II", "CN_DI"]
    derived_n = [n_L1N_LI + n_LBN_LI, n_G1N_GI + n_GBN_GI, n_CN_II, n_CN_DI]
    derived_freq = [freq_L1N_LI + freq_LBN_LI, freq_G1N_GI + freq_GBN_GI, freq_CN_II, freq_CN_DI]

    # ------------------------------------------------------------------
    # final table: ONLY these 8 categories, in this order (for nowa)
    # ------------------------------------------------------------------
    metrics = [ "CN_CI","CN_DI", "CN_II", "CN_LI", "CN_GI", "L1N_LI", "LBN_LI", "G1N_GI", "GBN_GI"]
    n_interactions = [n_CN_CI, n_CN_DI, n_CN_II,  n_CN_LI, n_CN_GI, n_L1N_LI, n_LBN_LI, n_G1N_GI, n_GBN_GI]
    freq_sum = [freq_CN_CI,freq_CN_DI, freq_CN_II,  freq_CN_LI, freq_CN_GI, freq_L1N_LI, freq_LBN_LI, freq_G1N_GI, freq_GBN_GI]

    result = pd.DataFrame({
        'metric': metrics,
        'n_interactions': n_interactions,
        'freq_sum': freq_sum,
    })

    # ------------------------------------------------------------------
    # Second table: count of permanent/exclusive nodes (species).
    # columns = pollinators, index = plants (convention used throughout the function).
    # ------------------------------------------------------------------
    node_metrics = ["PA", "PP", "A1", "P1", "A2", "P2"]
    node_counts_values = [
        len(common_col),   # PA: permanent pollinators
        len(common_row),   # PP: permanent plants
        len(cols_only_A),  # A1: transient pollinators (only in t1)
        len(rows_only_A),  # P1: transient plants (only in t1)
        len(cols_only_B),  # A2: transient pollinators (only in t2)
        len(rows_only_B),  # P2: transient plants (only in t2)
    ]
 
    node_counts = pd.DataFrame({
        'metric': node_metrics,
        'n_species': node_counts_values,
    })


    return result, node_counts