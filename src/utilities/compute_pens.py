import numpy as np
import scipy
from scipy import stats
import pandas as pd
import sqlite3

# Selection criteria:
# c1 is whether cohort 1 is included, analogous for c2, c3, c4

# lvwt_th and sbp_th divide the space of individuals into four quadrants
# quadrant selects which one is the cases

#         sbp_th
#           |
#      2    |    1
#           |
#   --------+-------- lvwt_th
#           |
#      3    |    4
#           |

### returns population, cases, controls given selection criteria
def get_data(bdc, c1, c2, c3, c4, lvwt_th, sbp_th, quadrant):
    cohort = []
    if c1:
        cohort.append("1")
    if c2:
        cohort.append("2")
    if c3:
        cohort.append("3")
    if c4:
        cohort.append("4")
    population = bdc[bdc["cohort"].isin(cohort)]
    
    if quadrant == 1:
        case_condition = (population["lvwt"] >= lvwt_th) & \
                         (population["systolic.bp"] >= sbp_th)
    elif quadrant == 2:
        case_condition = (population["lvwt"] <= lvwt_th) & \
                         (population["systolic.bp"] >= sbp_th)
    elif quadrant == 3:
        case_condition = (population["lvwt"] <= lvwt_th) & \
                         (population["systolic.bp"] <= sbp_th)
    else:
        case_condition = (population["lvwt"] >= lvwt_th) & \
                         (population["systolic.bp"] <= sbp_th)
    
    case = population[case_condition]
    control = population[~case_condition]
    return population, case, control

### computes penetrance given selection criteria
def compute_penetrance(bdc, c1, c2, c3, c4, lvwt_th, sbp_th, quadrant, mode, pathogenic_ids):
    population, case, control = get_data(bdc, c1, c2, c3, c4, lvwt_th, sbp_th, quadrant)
    
    num_dg = case[pathogenic_ids].any("columns").count()      # how many cases have pathogenic variants
    num_d = case.size                                         # how many cases
    num_g = population[pathogenic_ids].any("columns").count() # how many people have pathogenic variants
    num_all = population.size                                 # how many people total

    gld_rv = scipy.stats.beta(num_dg + 1, num_d - num_dg + 1)      # P(G|D)
    d_rv = scipy.stats.beta(num_d + 1, num_all - num_d + 1)        # P(D)
    g_rv = scipy.stats.beta(num_g + 1, num_all - num_g + 1)        # P(G)

    pen_samples = []
    for i in range(10000):
        gld = gld_rv.rvs()
        d = d_rv.rvs()
        g = g_rv.rvs()
        dlg = gld * d / g
        pen_samples.append(dlg)
    
    if mode == "map":
        hist, edges = np.histogram(pen_samples, bins=1000)
        mode_bin = np.argmax(hist)
        map_pen = edges[mode_bin]
        return map_pen
    elif mode == "5":
        return np.percentile(pen_samples, 0.025)
    elif mode == "95":
        return np.percentile(pen_samples, 0.975)


def compute_pens(bdc_database, variants_database):

    cnx_bdc_test = sqlite3.connect(bdc_database)
    cnx_variants = sqlite3.connect(variants_database)

    bdc = pd.read_sql_query("SELECT * FROM bdc", cnx_bdc_test)
    variants = pd.read_sql_query("SELECT * FROM variants", cnx_variants)
    clinvar = pd.read_sql_query("SELECT * FROM clinvar", cnx_variants)

    # find the rsID's of the pathogenic variants
    pathogenic_ids = []
    for c in bdc.columns:
        if c[0] == "r":
            chromosome = variants[variants["rsID"] == c]["Chromosome"].iloc[0]
            position = variants[variants["rsID"] == c]["Position"].iloc[0]
            clinvar_variants = clinvar[(clinvar["Chromosome"] == chromosome) &                                    (clinvar["Position"] == position)]
            pathogenic = any(clinvar_variants["Pathogenic"])
            if pathogenic:
                pathogenic_ids.append(c)


    # place to store results
    pens = pd.DataFrame(columns=["mode", "penetrance", 
                                 "cohort1", "cohort2", "cohort3", "cohort4",
                                 "LVWT threshold", "SBP threshold",
                                 "Quadrant"])

    choice_set = [(mode, c1, c2, c3, c4, lvwt_th, sbp_th, quad) \
                    for mode in ["map", "5", "95"] \
                    for c1 in [True, False] \
                    for c2 in [True, False] \
                    for c3 in [True, False] \
                    for c4 in [True, False] \
                    for lvwt_th in np.arange(0.8, 1, 0.03) \
                    for sbp_th in np.arange(0.8, 1, 0.03) \
                    for quad in [1]]
                    
    for choice in choice_set:
        mode, c1, c2, c3, c4, lvwt_th, sbp_th, quadrant = choice
        
        # skip if no cohorts are included
        if not (c1 or c2 or c3 or c4):
            continue
            
        pen = compute_penetrance(bdc, c1, c2, c3, c4, lvwt_th, sbp_th, quadrant, mode, pathogenic_ids)
        row = pd.DataFrame([[mode, pen, 
                             c1, c2, c3, c4, 
                             lvwt_th, sbp_th, quadrant]], 
                           columns=["mode", "penetrance", 
                                    "cohort1", "cohort2", "cohort3", "cohort4",
                                    "LVWT threshold", "SBP threshold",
                                    "Quadrant"])
        pens = pens.append(row, ignore_index=True)

    return pens
