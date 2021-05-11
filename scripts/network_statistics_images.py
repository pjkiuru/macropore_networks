# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 18:17:23 2020

@author: pkiuru

Script for calculating statistical tests for metwork metrics and drawing plots

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.stats as stats
import more_models as mm

try:
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    
    import statsmodels.stats.oneway as ao
    import statsmodels.stats.multicomp as mc
    import hypothetical as hp
    import scikit_posthocs as sp
except Exception:
    pass



fontname = 'Arial'
fontsize = 8

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

#rcParams['mathtext.fontset'] = 'cm'

rcParams['font.size'] = fontsize

rcParams['axes.titlesize'] = fontsize
rcParams['axes.labelsize'] = fontsize
rcParams['axes.titleweight'] = 'normal'

rcParams['xtick.labelsize'] = fontsize
rcParams['ytick.labelsize'] = fontsize

datafile = '../IO/METNET network measures.xlsx'
#Read data from file
dataread = 1

# Variance homogeneity tests by scipy.stats: Bartlett, Levene
variance_homogen = 0

# ANOVAS:
# (1) Simple standard ANOVA by scipy.stats
# (2) Standard ANOVA through OLS fit by statsmodels
# and residual normality tests by scipy.stats: (a) Shapiro, (b) Anderson, 
# (c) probability plot
# (3) Welch ANOVA by statsmodels 
# (4) Kruskal-Wallis test by scipy.stats
anova = 0

#Tukey pairwise multicomparison test by statsmodels
tukey = 0

#Games-Howell pairwise multicomparison test by hypothetical
games = 0

#Dunn pairwise multicomparison test by scikit_posthocs
dunn = 0

#Create boxplots of network metrics
boxpl = 1

#Create regression plots of network metrics
linfit = 1

#Calculate correlation coefficients between vertical tortuosity and other metrics
tort_corr = 1


#normality = 0


if dataread:
    
    df=pd.read_excel(datafile, sheet_name='Sheet1',
                     header=0)
    
    p_conn_1 = df['Poros. of conn. cluster'][0:7]
    p_conn_2 = df['Poros. of conn. cluster'][7:14]
    p_conn_3 = df['Poros. of conn. cluster'][14:21]
    
    np_1 = df['Number of pores'][0:7]
    np_2 = df['Number of pores'][7:14]
    np_3 = df['Number of pores'][14:21]
    
    cn_1 = df['Coordination number'][0:7]
    cn_2 = df['Coordination number'][7:14]
    cn_3 = df['Coordination number'][14:21]
    
    cc_1 = df['Clustering coefficient'][0:7]
    cc_2 = df['Clustering coefficient'][7:14]
    cc_3 = df['Clustering coefficient'][14:21]
    
    gt_1 = df['Geometr. tortuosity'][0:7]
    gt_2 = df['Geometr. tortuosity'][7:14].dropna()
    gt_3 = df['Geometr. tortuosity'][14:21]
    
    ccntr_1 = df['Closeness centrality'][0:7]
    ccntr_2 = df['Closeness centrality'][7:14]
    ccntr_3 = df['Closeness centrality'][14:21]
    
    bcntr_1 = df['Betweenness centr.'][0:7]
    bcntr_2 = df['Betweenness centr.'][7:14]
    bcntr_3 = df['Betweenness centr.'][14:21]
    
    tdbc_1 = df['Top-down betw. centr.'][0:7]
    tdbc_2 = df['Top-down betw. centr.'][7:14].dropna()
    tdbc_3 = df['Top-down betw. centr.'][14:21]
    
    pt_1 = df['Total porosity'][0:7]
    pt_2 = df['Total porosity'][7:14]
    pt_3 = df['Total porosity'][14:21]
    
    vfcc_1 = df['Vol. fract. of conn. cluster'][0:7]
    vfcc_2 = df['Vol. fract. of conn. cluster'][7:14]
    vfcc_3 = df['Vol. fract. of conn. cluster'][14:21]
    
    gt_mod = np.asarray(df['Geometr. tortuosity'].dropna())
    tdbc_mod = np.asarray(df['Top-down betw. centr.'].dropna())
    depths_mod = np.array([1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3])

#np_1 = stats.norm.rvs(10, 2, size=1000)
#np_2 = stats.norm.rvs(10, 2, size=1000)
#np_3 = stats.norm.rvs(10, 2, size=1000)
#np_2 = stats.gamma.rvs(0.8, loc=0, scale=1, size=10)

# Homogeneity of variance

if variance_homogen:

    
    print('')
    print('BARTLETT')
    print('')
    
    print('Porosity of connected cluster')
    
    print(stats.bartlett(p_conn_1, p_conn_2, p_conn_3))
    print([np.var(x, ddof=1) for x in [p_conn_1, p_conn_2, p_conn_3]])
    
    print('Number of pores')
    
    print(stats.bartlett(np_1, np_2, np_3))
    print([np.var(x, ddof=1) for x in [np_1, np_2, np_3]])
    
    print('Coordination number')
    
    print(stats.bartlett(cn_1, cn_2, cn_3))
    print([np.var(x, ddof=1) for x in [cn_1, cn_2, cn_3]])
    
    print('Clustering coefficient')
    
    print(stats.bartlett(cc_1, cc_2, cc_3))
    print([np.var(x, ddof=1) for x in [cc_1, cc_2, cc_3]])
    
    print('Geometric tortuosity')
    
    print(stats.bartlett(gt_1, gt_2, gt_3))
    print([np.var(x, ddof=1) for x in [gt_1, gt_2, gt_3]])
    
    print('Closeness centrality')
    
    print(stats.bartlett(ccntr_1, ccntr_2, ccntr_3))
    print([np.var(x, ddof=1) for x in [ccntr_1, ccntr_2, ccntr_3]])
    
    print('Betweenness centrality')
    
    print(stats.bartlett(bcntr_1, bcntr_2, bcntr_3))
    print([np.var(x, ddof=1) for x in [bcntr_1, bcntr_2, bcntr_3]])
    
    print('Top-down betweenness centrality')
    
    print(stats.bartlett(tdbc_1, tdbc_2, tdbc_3))
    print([np.var(x, ddof=1) for x in [tdbc_1, tdbc_2, tdbc_3]])
    
    print('Total porosity')
    
    print(stats.bartlett(pt_1, pt_2, pt_3))
    print([np.var(x, ddof=1) for x in [pt_1, pt_2, pt_3]])
    
    print('Volume fraction of connected cluster')
    
    print(stats.bartlett(vfcc_1, vfcc_2, vfcc_3)),
    print([np.var(x, ddof=1) for x in [vfcc_1, vfcc_2, vfcc_3]])

    
    print('')
    print('LEVENE')
    print('')
    
    print('Porosity of connected cluster')
    
    print(stats.levene(p_conn_1, p_conn_2, p_conn_3))
    print('2 & 3')
    print(stats.levene(p_conn_2, p_conn_3))
    print([np.var(x, ddof=1) for x in [p_conn_1, p_conn_2, p_conn_3]])
    
    print('Number of pores')
    
    print(stats.levene(np_1, np_2, np_3))
    print('2 & 3')
    print(stats.levene(np_2, np_3))
    print([np.var(x, ddof=1) for x in [np_1, np_2, np_3]])
    
    print('Coordination number')
    
    print(stats.levene(cn_1, cn_2, cn_3))
    print('2 & 3')
    print(stats.levene(cn_2, cn_3))
    print([np.var(x, ddof=1) for x in [cn_1, cn_2, cn_3]])
    
    print('Clustering coefficient')
    
    print(stats.levene(cc_1, cc_2, cc_3))
    print('2 & 3')
    print(stats.levene(cc_2, cc_3))
    print([np.var(x, ddof=1) for x in [cc_1, cc_2, cc_3]])
    
    print('Geometric tortuosity')
    
    print(stats.levene(gt_1, gt_2, gt_3))
    print('2 & 3')
    print(stats.levene(gt_2, gt_3))
    print([np.var(x, ddof=1) for x in [gt_1, gt_2, gt_3]])
    
    print('Closeness centrality')
    
    print(stats.levene(ccntr_1, ccntr_2, ccntr_3))
    print('2 & 3')
    print(stats.levene(ccntr_2, ccntr_3))
    print([np.var(x, ddof=1) for x in [ccntr_1, ccntr_2, ccntr_3]])
    
    print('Betweenness centrality')
    
    print(stats.levene(bcntr_1, bcntr_2, bcntr_3))
    print('2 & 3')
    print(stats.levene(bcntr_2, bcntr_3))
    print([np.var(x, ddof=1) for x in [bcntr_1, bcntr_2, bcntr_3]])
    
    print('Top-down betweenness centrality')
    
    print(stats.levene(tdbc_1, tdbc_2, tdbc_3))
    print('2 & 3')
    print(stats.levene(tdbc_2, tdbc_3))
    print([np.var(x, ddof=1) for x in [tdbc_1, tdbc_2, tdbc_3]])
    
    print('Total porosity')
    
    print(stats.levene(pt_1, pt_2, pt_3))
    print('2 & 3')
    print(stats.levene(pt_2, pt_3))
    print([np.var(x, ddof=1) for x in [pt_1, pt_2, pt_3]])
    
    print('Volume fraction of connected cluster')
    
    print(stats.levene(vfcc_1, vfcc_2, vfcc_3))
    print('2 & 3')
    print(stats.levene(vfcc_2, vfcc_3))
    print([np.var(x, ddof=1) for x in [vfcc_1, vfcc_2, vfcc_3]])
# ANOVA


if anova:
    
    print('')
    print('Standard ANOVA, F-ONEWAY')
    print('')
    
    print('Porosity of connected cluster')
    
    print(stats.f_oneway(p_conn_1, p_conn_2, p_conn_3))
    
    print('Number of pores')
    
    print(stats.f_oneway(np_1, np_2, np_3))
    
    print('Coordination number')
    
    print(stats.f_oneway(cn_1, cn_2, cn_3))
    
    print('Clustering coefficient')
    
    print(stats.f_oneway(cc_1, cc_2, cc_3))
    
    print('Geometric tortuosity')
    
    print(stats.f_oneway(gt_1, gt_2, gt_3))
    
    print('Closeness centrality')
    
    print(stats.f_oneway(ccntr_1, ccntr_2, ccntr_3))
    
    print('Betweenness centrality')
    
    print(stats.f_oneway(bcntr_1, bcntr_2, bcntr_3))
    
    print('Top-down betweenness centrality')
    
    print(stats.f_oneway(tdbc_1, tdbc_2, tdbc_3))
    
    print('Total porosity')
    
    print(stats.f_oneway(pt_1, pt_2, pt_3))
    
    print('Volume fraction of connected cluster')
    
    print(stats.f_oneway(vfcc_1, vfcc_2, vfcc_3))
    
    
    
    if True:
        print('')
        print('Standard ANOVA ols fit')
        print('')
        
        for w in range(1,12):
            df.rename(columns={df.columns[w]: df.columns[w].replace(" ", "x")},inplace=True)
            df.rename(columns={df.columns[w]: df.columns[w].replace(".", "q")},inplace=True)
            df.rename(columns={df.columns[w]: df.columns[w].replace("-", "Q")},inplace=True)
        
        st_p_conn = ols('Porosqxofxconnqxcluster ~ C(Depth)', data=df).fit()
        print('-------- Porosity of connected cluster ----------')
        print(sm.stats.anova_lm(st_p_conn, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_p_conn.resid)))
        print(stats.anderson(st_p_conn.resid))
        mm.probability_plot_ols(1, st_p_conn.resid)
        
        st_np = ols('Numberxofxpores ~ C(Depth)', data=df).fit()
        print('-------- Number of pores ----------')
        print(sm.stats.anova_lm(st_np, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_np.resid)))
        print(stats.anderson(st_np.resid))
        mm.probability_plot_ols(2, st_np.resid)
        
        st_cn = ols('Coordinationxnumber ~ C(Depth)', data=df).fit()
        print('-------- Coordination number ----------')
        print(sm.stats.anova_lm(st_cn, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_cn.resid)))
        print(stats.anderson(st_cn.resid))
        mm.probability_plot_ols(3, st_cn.resid)
        
        st_cc = ols('Clusteringxcoefficient ~ C(Depth)', data=df).fit()
        print('-------- Clustering coefficient ----------')
        print(sm.stats.anova_lm(st_cc, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_cc.resid)))
        print(stats.anderson(st_cc.resid))
        mm.probability_plot_ols(4, st_cc.resid)
        
        st_gt = ols('Geometrqxtortuosity ~ C(Depth)', data=df).fit()
        print('-------- Geometric tortuosity ----------')
        print(sm.stats.anova_lm(st_gt, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_gt.resid)))
        print(stats.anderson(st_gt.resid))
        mm.probability_plot_ols(5, st_gt.resid)
        
        st_ccntr = ols('Closenessxcentrality ~ C(Depth)', data=df).fit()
        print('-------- Closeness centrality ----------')
        print(sm.stats.anova_lm(st_ccntr, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_ccntr.resid)))
        print(stats.anderson(st_ccntr.resid))
        mm.probability_plot_ols(6, st_ccntr.resid)
        
        
        st_bcntr = ols('Betweennessxcentrq ~ C(Depth)', data=df).fit()
        print('-------- Betweenness centrality ----------')
        print(sm.stats.anova_lm(st_bcntr, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_bcntr.resid)))
        print(stats.anderson(st_bcntr.resid))
        mm.probability_plot_ols(7, st_bcntr.resid)
        
        st_tdbc = ols('TopQdownxbetwqxcentrq ~ C(Depth)', data=df).fit()
        print('-------- Top-down betweenness centrality ----------')
        print(sm.stats.anova_lm(st_tdbc, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_tdbc.resid)))
        print(stats.anderson(st_tdbc.resid))
        mm.probability_plot_ols(8, st_tdbc.resid)
        
        st_pt = ols('Totalxporosity ~ C(Depth)', data=df).fit()
        print('-------- Total porosity ----------')
        print(sm.stats.anova_lm(st_pt, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_pt.resid)))
        print(stats.anderson(st_pt.resid))
        mm.probability_plot_ols(9, st_pt.resid)
        
        st_vfcc = ols('Volqxfractqxofxconnqxcluster ~ C(Depth)', data=df).fit()
        print('-------- Vol. fract of conn. cluster ----------')
        print(sm.stats.anova_lm(st_vfcc, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_vfcc.resid)))
        print(stats.anderson(st_vfcc.resid))
        mm.probability_plot_ols(10, st_vfcc.resid)
        
        for w in range(1,12):
            df.rename(columns={df.columns[w]: df.columns[w].replace("x", " ")},inplace=True)
            df.rename(columns={df.columns[w]: df.columns[w].replace("q", ".")},inplace=True)
            df.rename(columns={df.columns[w]: df.columns[w].replace("Q", "-")},inplace=True)
    
    #except Exception:
    
    #    pass
    
    
    try:
        print('')
        print('Welch ANOVA')
        print('')
        
        welch = ao.anova_oneway(df['Poros. of conn. cluster'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Porosity of connected cluster ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Number of pores'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Number of pores ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Coordination number'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Coordination number ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Clustering coefficient'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Clustering coefficient ----------')
        print(welch)
        
        welch = ao.anova_oneway(gt_mod, depths_mod,
                           use_var='unequal', welch_correction=True)
        print('-------- Geometric tortuosity ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Closeness centrality'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Closeness centrality ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Betweenness centr.'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Betweenness centrality ----------')
        print(welch)
        
        welch = ao.anova_oneway(tdbc_mod, depths_mod,
                           use_var='unequal', welch_correction=True)
        print('-------- Top-down betweenness centrality ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Total porosity'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Total porosity ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['Vol. fract. of conn. cluster'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- Vol. fract of conn. cluster ----------')
        print(welch)
    
    except Exception:
    
        pass
    
    
    print('')
    print('KRUSKAL')
    print('')
    
    print('Porosity of connected cluster')
    
    print(stats.kruskal(p_conn_1, p_conn_2, p_conn_3))
    
    print('Number of pores')
    
    print(stats.kruskal(np_1, np_2, np_3))
    
    print('Coordination number')
    
    print(stats.kruskal(cn_1, cn_2, cn_3))
    
    print('Clustering coefficient')
    
    print(stats.kruskal(cc_1, cc_2, cc_3))
    
    print('Geometric tortuosity')
    
    print(stats.kruskal(gt_1, gt_2, gt_3))
    
    print('Closeness centrality')
    
    print(stats.kruskal(ccntr_1, ccntr_2, ccntr_3))
    
    print('Betweenness centrality')
    
    print(stats.kruskal(bcntr_1, bcntr_2, bcntr_3))
    
    print('Top-down betweenness centrality')
    
    print(stats.kruskal(tdbc_1, tdbc_2, tdbc_3))
    
    print('Total porosity')
    
    print(stats.kruskal(pt_1, pt_2, pt_3))
    
    print('Volume fraction of connected cluster')
    
    print(stats.kruskal(vfcc_1, vfcc_2, vfcc_3))

if tukey:
    
    print('Tukey')
    print('')
    
    comp = mc.MultiComparison(df['Poros. of conn. cluster'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Poros. of conn. cluster')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Number of pores'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Number of pores')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Coordination number'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Coordination number')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Clustering coefficient'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Clustering coefficient')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(gt_mod, depths_mod)
    post_hoc_res = comp.tukeyhsd()
    print('Geometr. tortuosity')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Closeness centrality'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Closeness centrality')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Betweenness centr.'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Betweenness centr.')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(tdbc_mod, groups=depths_mod)
    post_hoc_res = comp.tukeyhsd()
    print('Top-down betw. centr.')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Total porosity'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Total porosity')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['Vol. fract. of conn. cluster'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('Vol. fract. of conn. cluster')
    print(post_hoc_res.summary())
    

if games:
    
    print('Games')
    print('')
    
    hp_games = hp.posthoc.GamesHowell(df['Poros. of conn. cluster'], group=df['Depth'], alpha=0.05)
    print('Porosity of connected cluster')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Number of pores'], group=df['Depth'], alpha=0.05)
    print('Number of pores')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Coordination number'], group=df['Depth'], alpha=0.05)
    print('Coordination number')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Clustering coefficient'], group=df['Depth'], alpha=0.05)
    print('Clustering coefficient')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(gt_mod, group=depths_mod, alpha=0.05)
    print('Geometric tortuosity')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Closeness centrality'], group=df['Depth'], alpha=0.05)
    print('Closeness centrality')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Betweenness centr.'], group=df['Depth'], alpha=0.05)
    print('Betweenness centrality')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(tdbc_mod, group=depths_mod, alpha=0.05)
    print('Top-down betweenness centrality')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Total porosity'], group=df['Depth'], alpha=0.05)
    print('Total porosity')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['Vol. fract. of conn. cluster'], group=df['Depth'], alpha=0.05)
    print('Volume fraction of connected cluster')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
if dunn:
    
    print('')
    print('Dunn-Bonferroni')
    print('')
    
    print('Coordination number')
    ph_dunn = sp.posthoc_dunn(df, val_col='Coordination number', group_col='Depth', p_adjust = 'holm')    
    print(ph_dunn)
    
    print('Betweenness centrality')
    ph_dunn = sp.posthoc_dunn(df, val_col='Betweenness centr.', group_col='Depth', p_adjust = 'holm')    
    print(ph_dunn)
    
    print('Top-down betweenness centrality')
    ph_dunn = sp.posthoc_dunn(df, val_col='Top-down betw. centr.', group_col='Depth', p_adjust = 'holm')    
    print(ph_dunn)
    
    print('Vol. fract. of conn. cluster')
    ph_dunn = sp.posthoc_dunn(df, val_col='Vol. fract. of conn. cluster', group_col='Depth', p_adjust = 'holm')    
    print(ph_dunn)
    
if boxpl:
    '''
    for i in range(2,3):
    
        fig = plt.figure(num=100+i-1)
        fig.set_size_inches(3,3)
        plt.clf()
        ax = fig.add_subplot(111)
        
        data = [df[df.columns[i]][df['Depth'] == 1],
                df[df.columns[i]][df['Depth'] == 2].dropna(),
                df[df.columns[i]][df['Depth'] == 3]]
        
        ax.boxplot(data, labels= ['0-5 cm', '20-25 cm', '40-45 cm'],
                   showmeans= True)
        
        ax.set_xlabel("Depth")
        ax.set_ylabel(df.columns[i])
        
        plt.tight_layout()
     '''
    bbdim = [0.14, 0.36]
    
    bppos = [
            [0.06, 0.57, bbdim[0],bbdim[1]],
            [0.255, 0.57, bbdim[0],bbdim[1]],
            [0.45, 0.57, bbdim[0],bbdim[1]],
            [0.645, 0.57, bbdim[0],bbdim[1]],
            [0.84, 0.57, bbdim[0],bbdim[1]],
            [0.06, 0.095, bbdim[0],bbdim[1]],
            [0.255, 0.095, bbdim[0],bbdim[1]],
            [0.45, 0.095, bbdim[0],bbdim[1]],
            [0.645, 0.095, bbdim[0],bbdim[1]],
            [0.84, 0.095, bbdim[0],bbdim[1]]]
    
    subtpos = [0.01,1.06]
    
    subt = [('(a)'), ('(b)'), ('(c)'), ('(d)'),
            ('(e)'), ('(f)'), ('(g)'), ('(h)'), ('(i)')]
    
    ylabs = [
            'Network porosity',
            'Number of pores (\u00D710$^3$)',
            'Coordination number',
            'Clustering coefficient',
            'Vertical geom. tortuosity',
            'Closeness centrality',
            'Betweenness centrality (\u00D710$^{-3}$)',
            'Top\u2013bottom betw. centr. (\u00D710$^{-3}$)',
            '',
            '',
            '',
            '',
            'Horizontal geom. tortuosity',
            ]
    
    ylc = [-0.22, -0.16, -0.16, -0.21 ,-0.16,
           -0.18, -0.20, -0.16, -0.00, -0.00, 
           -0.00, -0.00, -0.21]
    
    medianprops = dict(linestyle='-', linewidth=1, color='blue')
    flierprops = dict(marker='o', markersize=4)#, linestyle='none')
    
    fig = plt.figure(num=101)
    fig.set_size_inches(7,3.2)
    plt.clf()
    
    bporder = np.array([2,3,4,5,7,8,9,6,14])
    
    for j in range(0,9):
        
        ax = fig.add_subplot(2,5,j+1)
        
        i = bporder[j]
        
        if i == 3:
            factor = 0.001
        elif i == 8 or i == 9:
            factor = 1000
        else:
            factor = 1
        
        data = [factor*df[df.columns[i]][df['Depth'] == 1],
                factor*df[df.columns[i]][df['Depth'] == 2].dropna(),
                factor*df[df.columns[i]][df['Depth'] == 3]]
        
        ax.boxplot(data, labels= ['0-5', '20-25', '40-45'],
                   showmeans=False, medianprops=medianprops,
                   flierprops=flierprops)
        
        
        
        if i-2 > -2:
            ax.text(0.13, 0.9, "A", transform=ax.transAxes)
            ax.text(0.47, 0.9, "B", transform=ax.transAxes)
            ax.text(0.8, 0.9, "B", transform=ax.transAxes)
        
        #ax.set_title(subt[i-2], horizontalalignment='right')
        
        ax.text(subtpos[0], subtpos[1], subt[j], transform=ax.transAxes)
        
        if i==6:
            ax.set_ylim([0, 1.2 * factor * np.max(df[df.columns[6]])])
        elif i==14:
            ax.set_ylim([0, 1.2 * factor * np.max(df[df.columns[6]])])
        else:
            ax.set_ylim([0, 1.2 * factor * np.max(df[df.columns[i]])])
        
        if i-1 > 4:
            ax.set_xlabel("Depth (cm)")
            ax.xaxis.set_label_coords(0.5, -0.15)
         
        ax.xaxis.set_tick_params(length=0)
        
        ax.set_ylabel(ylabs[i-2], labelpad=0.9)
        ax.yaxis.set_label_coords(ylc[i-2], 0.5)
        
        ax.set_position(bppos[j])
        
        #plt.savefig('netmetrics.pdf')
    
if linfit:
    
    bbdim = [0.41, 0.18]
    
    ds = ['0\u20135 cm', '20\u201325 cm', '40\u201345 cm']
    
    bppos = [
            [0.07, 0.765, bbdim[0],bbdim[1]],
            [0.56, 0.765, bbdim[0],bbdim[1]],
            [0.07, 0.525, bbdim[0],bbdim[1]],
            [0.56, 0.525, bbdim[0],bbdim[1]],
            [0.07, 0.285, bbdim[0],bbdim[1]],
            [0.56, 0.285, bbdim[0],bbdim[1]],
            [0.07, 0.045, bbdim[0],bbdim[1]],
            [0.56, 0.045, bbdim[0],bbdim[1]],
            ]
    
    subt = [('(a)'), ('(b)'), ('(c)'), ('(d)'),
            ('(e)'), ('(f)'), ('(g)'), ('(h)')]
    
    subtpos = [0.03, 0.92]
    
    regrvertpos = [0.91, 0.81, 0.71]
    
    ylabs = [
            'Connected cluster porosity',
            'Number of pores (\u00D710$^3$)',
            'Coordination number',
            'Clustering coefficient',
            'Geometric tortuosity',
            'Closeness centrality',
            'Betweenness centrality (\u00D710$^{-3}$)',
            'Top-bottom betw. centr. (\u00D710$^{-3}$)'
            ]
    
    lineqs = {
            "P_vs_Np_1": "y = 0.041x \u2212 0.068, R$^2$ = 0.75",
            "P_vs_Np_2": "y = 0.018x \u2212 0.002, R$^2$ = 0.91",
            "P_vs_Np_3": "y = 0.021x \u2212 0.016, R$^2$ = 0.94",
            
            "P_vs_CN_1": "y = 8.0x + 3.3, R$^2$ = 0.86",
            "P_vs_CN_2": "y = 7.0x + 2.9, R$^2$ = 0.40",
            "P_vs_CN_3": "y = 4.5x + 2.7, R$^2$ = 0.57",
            
            "CN_vs_CC_1": "y = 0.056x + 0.093, R$^2$ = 0.95",
            "CN_vs_CC_2": "y = 0.071x + 0.012, R$^2$ = 0.43",
            "CN_vs_CC_3": "y = 0.10x \u2212 0.083, R$^2$ = 0.46",
            
            "P_vs_GT_1": "y = \u22120.72x + 1.8, R$^2$ = 0.62",
            "P_vs_GT_2": "y = 23x + 1.1, R$^2$ = 0.96",
            "P_vs_GT_3": "y = \u22121.0x + 2.4, R$^2$ = 0.01",
            
            "GT_vs_GTH_1": "y = 0.55x + 0.62, R$^2$ = 0.91",
            "GT_vs_GTH_2": "y = \u22120.19x + 2.4, R$^2$ = 0.33",
            "GT_vs_GTH_3": "y = \u22120.076x + 2.3, R$^2$ = 0.01",
            
            "CN_vs_Ccntr_1": "y = 1.5x + 28, R$^2$ = 0.37",
            "CN_vs_Ccntr_2": "y = 7.5x + 1.2, R$^2$ = 0.38",
            "CN_vs_Ccntr_3": "y = 9.9x \u2212 3.9, R$^2$ = 0.65",
            
            "CN_vs_Bcntr_1": "y = \u22120.48x + 4.4, R$^2$ = 0.30",
            "CN_vs_Bcntr_2": "y = \u22125.3x + 27, R$^2$ = 0.08",
            "CN_vs_Bcntr_3": "y = \u221213x + 46, R$^2$ = 0.75",
            
            "GT_vs_TDBC_1": "y = 6.6x \u2212 7.5, R$^2$ = 0.17",
            "GT_vs_TDBC_2": "y = \u22123.2x + 20, R$^2$ = 0.50",
            "GT_vs_TDBC_3": "y = 2.1x + 5.4, R$^2$ = 0.02",
            
            }
    
    fig = plt.figure(num=202)
    fig.set_size_inches(6,7)
    plt.clf()
    
    ax1 = fig.add_subplot(4,2,1)
    print('')
    print('P - Np')
    print('')
    dummyx = 0.001*df['Number of pores'][df['Depth'] == 1]
    dummyy = df['Poros. of conn. cluster'][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax1.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0] ,zorder=4)
    ax1.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit', zorder=1)
    ax1.text(0.12, regrvertpos[0], ds[0] + ': ' + lineqs['P_vs_Np_1'],
             transform=ax1.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = 0.001*df['Number of pores'][df['Depth'] == 2]
    dummyy = df['Poros. of conn. cluster'][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax1.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1] ,zorder=5)
    ax1.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit', zorder=2)
    ax1.text(0.12, regrvertpos[1], ds[1] + ': ' + lineqs['P_vs_Np_2'],
             transform=ax1.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = 0.001*df['Number of pores'][df['Depth'] == 3]
    dummyy = df['Poros. of conn. cluster'][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax1.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2] ,zorder=6)
    ax1.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit' ,zorder=3)
    ax1.text(0.12, regrvertpos[2], ds[2] + ': ' + lineqs['P_vs_Np_3'],
             transform=ax1.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    ax1.set_xlim([0, 14])
    ax1.set_ylim([-0.05, 0.88])
    
    ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)
    
    #handles, labels = ax1.get_legend_handles_labels()
    #order = [3,4,5]
    #ax1.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax1.legend(fontsize=fontsize)
    
    ax1.set_xlabel('Number of pores (\u00D710$^3$)')
    ax1.set_ylabel('Network porosity')
    
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    ax1.set_position(bppos[0])   
    
    #----------
    
    ax2 = fig.add_subplot(4,2,2)
    
    xcol = 'Poros. of conn. cluster'
    ycol = 'Coordination number'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax2.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0] ,zorder=4)
    ax2.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit' ,zorder=1)
    ax2.text(0.12, regrvertpos[0], ds[0] + ': ' + lineqs['P_vs_CN_1'],
             transform=ax2.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax2.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1] ,zorder=5)
    ax2.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit' ,zorder=2)
    ax2.text(0.12, regrvertpos[1], ds[1] + ': ' + lineqs['P_vs_CN_2'],
             transform=ax2.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax2.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2] ,zorder=6)
    ax2.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit' ,zorder=3)
    ax2.text(0.12, regrvertpos[2], ds[2] + ': ' + lineqs['P_vs_CN_3'],
             transform=ax2.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    ax2.set_xlim([0, 0.6])
    ax2.set_ylim([1, 11])
    ax2.yaxis.set_label_coords(-0.07, 0.5, transform=ax2.transAxes)
    
    ax2.text(subtpos[0], subtpos[1], subt[1], transform=ax2.transAxes)
    
    #handles, labels = ax2.get_legend_handles_labels()
    #order = [3,4,5]
    #ax2.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax2.legend(fontsize=fontsize)
    
    ax2.set_xlabel('Network porosity')
    ax2.set_ylabel('Coordination number')
    
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    ax2.set_position(bppos[1])    
    
    #-------------------
    
    ax3 = fig.add_subplot(4,2,3)
    
    xcol = 'Coordination number'
    ycol = 'Clustering coefficient'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax3.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0], zorder=4)
    ax3.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit' ,zorder=1)
    ax3.text(0.12, regrvertpos[0], ds[0] + ': ' + lineqs['CN_vs_CC_1'],
             transform=ax3.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax3.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1] ,zorder=5)
    ax3.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit' ,zorder=2)
    ax3.text(0.12, regrvertpos[1], ds[1] + ': ' + lineqs['CN_vs_CC_2'],
             transform=ax3.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax3.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2] ,zorder=6)
    ax3.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit' ,zorder=3)
    ax3.text(0.12, regrvertpos[2], ds[2] + ': ' + lineqs['CN_vs_CC_3'],
             transform=ax3.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    
    ax3.set_xlim([2, 9])
    ax3.set_ylim([0.1, 0.8])
    
    ax3.text(subtpos[0], subtpos[1], subt[2], transform=ax3.transAxes)
    
    #handles, labels = ax3.get_legend_handles_labels()
    #order = [3,4,5]
    #ax3.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax3.legend(fontsize=fontsize)
    
    ax3.set_xlabel('Coordination number')
    ax3.set_ylabel('Clustering coefficient')
    
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    
    ax3.set_position(bppos[2])
    
    #----------
    
    
    ax4 = fig.add_subplot(4,2,4)
    
    xcol = 'Coordination number'
    ycol = 'Closeness centrality'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax4.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0], zorder=4)
    ax4.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit', zorder=1)
    ax4.text(0.12, regrvertpos[0], ds[0] + ': ' + lineqs['CN_vs_Ccntr_1'],
             transform=ax4.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax4.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1], zorder=5)
    ax4.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit', zorder=2)
    ax4.text(0.12, regrvertpos[1], ds[1] + ': ' + lineqs['CN_vs_Ccntr_2'],
             transform=ax4.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax4.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2], zorder=6)
    ax4.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit', zorder=3)
    ax4.text(0.12, regrvertpos[2], ds[2] + ': ' + lineqs['CN_vs_Ccntr_3'],
             transform=ax4.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    
    ax4.set_xlim([2, 9])
    ax4.set_ylim([15.0, 55.0])
    
    ax4.text(subtpos[0], subtpos[1], subt[3], transform=ax4.transAxes)
    
    #handles, labels = ax4.get_legend_handles_labels()
    #order = [3,4,5]
    #ax4.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax4.legend(fontsize=fontsize)
    
    ax4.set_xlabel('Coordination number')
    ax4.set_ylabel('Closeness centrality')
    
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    
    ax4.set_position(bppos[3])    
    
    #----------
    
    ax5 = fig.add_subplot(4,2,5)
    
    xcol = 'Coordination number'
    ycol = 'Betweenness centr.'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = 1000 * df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax5.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax5.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax5.text(0.25, regrvertpos[0], ds[0] + ': ' + lineqs['CN_vs_Bcntr_1'],
             transform=ax5.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = 1000 * df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax5.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax5.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit')
    ax5.text(0.25, regrvertpos[1], ds[1] + ': ' + lineqs['CN_vs_Bcntr_2'],
             transform=ax5.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = 1000 * df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax5.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2])
    ax5.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit')
    ax5.text(0.25, regrvertpos[2], ds[2] + ': ' + lineqs['CN_vs_Bcntr_3'],
             transform=ax5.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    
    ax5.set_xlim([2, 9])
    ax5.set_ylim([-1, 18])
    ax5.yaxis.set_label_coords(-0.06, 0.5, transform=ax5.transAxes)
    
    handles, labels = ax5.get_legend_handles_labels()
    order = [3,4,5]
    ax5.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
               ncol = 3, fontsize=fontsize,  bbox_to_anchor=(0.98, 0.99),
               bbox_transform=plt.gcf().transFigure,
               columnspacing=1, handletextpad=0.25)
    
    #ax5.legend(fontsize=fontsize)
    
    ax5.set_xlabel('Coordination number')
    ax5.set_ylabel('Betweenness centrality (\u00D710$^{-3}$)')
    
    ax5.text(subtpos[0], subtpos[1], subt[4], transform=ax5.transAxes)
    
    ax5.spines["top"].set_visible(False)
    ax5.spines["right"].set_visible(False)
    
    ax5.set_position(bppos[4])    
    
     #--------------------
    
    ax6 = fig.add_subplot(4,2,6)
    
    xcol = 'Poros. of conn. cluster'
    ycol = 'Geometr. tortuosity'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax6.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax6.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax6.text(0.25, regrvertpos[0], ds[0] + ': ' + lineqs['P_vs_GT_1'],
             transform=ax6.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][np.array([0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0]).astype('bool')]
    dummyy = df[ycol][df['Depth'] == 2].dropna()
    lr_1 = stats.linregress(dummyx,dummyy)
    ax6.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax6.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit')
    ax6.text(0.25, regrvertpos[1], ds[1] + ': ' + lineqs['P_vs_GT_2'],
             transform=ax6.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax6.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2])
    ax6.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit')
    ax6.text(0.25, regrvertpos[2], ds[2] + ': ' + lineqs['P_vs_GT_3'],
             transform=ax6.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    
    ax6.set_xlim([0, 0.6])
    ax6.set_ylim([0.9, 4.0])
    ax6.yaxis.set_label_coords(-0.07, 0.5, transform=ax6.transAxes)
    
    ax6.text(subtpos[0], subtpos[1], subt[5], transform=ax6.transAxes)
    
    #handles, labels = ax6.get_legend_handles_labels()
    #order = [3,4,5]
    #ax6.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax6.legend(fontsize=fontsize)
    
    ax6.set_xlabel('Network porosity')
    ax6.set_ylabel('Vertical geometrical tortuosity')
    
    ax6.spines["top"].set_visible(False)
    ax6.spines["right"].set_visible(False)
    
    ax6.set_position(bppos[5])    
    
     #--------------------
    
    ax7 = fig.add_subplot(4,2,7)
    
    xcol = 'Geometr. tortuosity'
    ycol = 'Geomtorthorall'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax7.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax7.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax7.text(0.12, regrvertpos[0], ds[0] + ': ' + lineqs['GT_vs_GTH_1'],
             transform=ax7.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 2].dropna()
    dummyy = df[ycol][np.array([0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0]).astype('bool')]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax7.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax7.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit')
    ax7.text(0.12, regrvertpos[1], ds[1] + ': ' + lineqs['GT_vs_GTH_2'],
             transform=ax7.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax7.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2])
    ax7.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit')
    ax7.text(0.12, regrvertpos[2], ds[2] + ': ' + lineqs['GT_vs_GTH_3'],
             transform=ax7.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    
    ax7.set_xlim([0.9, 4.0])
    ax7.set_ylim([0.9, 4.0])
    ax7.yaxis.set_label_coords(-0.07, 0.5, transform=ax7.transAxes)
    
    ax7.text(subtpos[0], subtpos[1], subt[6], transform=ax7.transAxes)
    
    #handles, labels = ax7.get_legend_handles_labels()
    #order = [3,4,5]
    #ax7.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax7.legend(fontsize=fontsize)
    
    ax7.set_xlabel('Vertical geometrical tortuosity')
    ax7.set_ylabel('Horizontal geometrical tortuosity')
    
    ax7.spines["top"].set_visible(False)
    ax7.spines["right"].set_visible(False)
    
    ax7.set_position(bppos[6])    
    
    #----------
    
    ax8 = fig.add_subplot(4,2,8)
    
    xcol = 'Geometr. tortuosity'
    ycol = 'Top-down betw. centr.'
    print('')
    print(xcol + ' - ' + ycol)
    print('')
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = 1000 * df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax8.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax8.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax8.text(0.12, regrvertpos[0], ds[0] + ': ' + lineqs['GT_vs_TDBC_1'],
             transform=ax8.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 2].dropna()
    dummyy = 1000 * df[ycol][df['Depth'] == 2].dropna()
    lr_1 = stats.linregress(dummyx,dummyy)
    ax8.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax8.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b-', lw=1, label='fit')
    ax8.text(0.12, regrvertpos[1], ds[1] + ': ' + lineqs['GT_vs_TDBC_2'],
             transform=ax8.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = 1000 * df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax8.scatter(dummyx, dummyy, marker='^', c='r', s=10, label=ds[2])
    ax8.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-', lw=1, label='fit')
    ax8.text(0.12, regrvertpos[2], ds[2] + ': ' + lineqs['GT_vs_TDBC_3'],
             transform=ax8.transAxes, fontsize=fontsize)
    print(lr_1)
    print('R2 = ' + str(lr_1.rvalue**2))
    print('')
    
    ax8.set_xlim([1, 4])
    ax8.set_ylim([-1, 26])
    
    ax8.text(subtpos[0], subtpos[1], subt[7], transform=ax8.transAxes)
    
    #handles, labels = ax8.get_legend_handles_labels()
    #order = [3,4,5]
    #ax8.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, bbox_to_anchor=(0.8, 0.1),
    #           bbox_transform=plt.gcf().transFigure)
    
    #ax8.legend(fontsize=fontsize)
    
    ax8.set_xlabel('Vertical geometrical tortuosity')
    ax8.set_ylabel('Top\u2013bottom betw. centr. (\u00D710$^{-3}$)')
    
    ax8.spines["top"].set_visible(False)
    ax8.spines["right"].set_visible(False)
    
    ax8.set_position(bppos[7])    
    
    
    #----------
    
    ax1.xaxis.set_label_coords(0.5, -0.11)
    ax2.xaxis.set_label_coords(0.5, -0.13)
    ax3.xaxis.set_label_coords(0.5, -0.13)
    ax4.xaxis.set_label_coords(0.5, -0.13)
    ax5.xaxis.set_label_coords(0.5, -0.13)
    ax6.xaxis.set_label_coords(0.5, -0.13)
    ax7.xaxis.set_label_coords(0.5, -0.13)
    ax8.xaxis.set_label_coords(0.5, -0.13)
    
    ax1.yaxis.set_label_coords(-0.10, 0.5)
    ax3.yaxis.set_label_coords(-0.10, 0.5)
    ax5.yaxis.set_label_coords(-0.10, 0.5)
    ax7.yaxis.set_label_coords(-0.10, 0.5)
    
    #ax5.yaxis jää pdf:ssä liian lähelle, jos x-arvo on pienempi kuin muilla;
    #kuvaikkunassa on eri asetus!
    
    ax2.yaxis.set_label_coords(-0.08, 0.5)
    ax4.yaxis.set_label_coords(-0.08, 0.5)
    ax6.yaxis.set_label_coords(-0.08, 0.5)
    ax8.yaxis.set_label_coords(-0.08, 0.5)
    
    #plt.savefig('metricscorrs.pdf')
    
    #----------

if tort_corr:
    
    metrs = ['Coordination number', 'Clustering coefficient', 'Closeness centrality' ,'Betweenness centr.'] 
    ds = [1,2,3]
    
    for k, l in enumerate(metrs):
    
        xcol = metrs[k]
        ycol = 'Geometr. tortuosity'
        print('')
        print(xcol + ' - ' + ycol)
        print('')
        for m in ds:
            if m == 2:
                dummyx = df[xcol][np.array([0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0]).astype('bool')]
                dummyy = df[ycol][df['Depth'] == 2].dropna()
            else:
                dummyx = df[xcol][df['Depth'] == m]
                dummyy = df[ycol][df['Depth'] == m]
            lr_1 = stats.linregress(dummyx,dummyy)
            print(m)
            print(lr_1.rvalue)
            print(lr_1.pvalue)
            print('')
        
'''


comp = mc.MultiComparison(df['Closeness centrality'], df['Depth'])
post_hoc_res = comp.tukeyhsd()
print(post_hoc_res.summary())


#Welch ANOVA https://statisticsbyjim.com/anova/welchs-anova-compared-to-classic-one-way-anova/
#https://www.statsmodels.org/devel/generated/statsmodels.stats.oneway.anova_oneway.html#statsmodels.stats.oneway.anova_oneway
aone = ao.anova_oneway(df['Closeness centrality'], df['Depth'],
                       use_var='unequal', welch_correction=True)
print(aone)
#Standard ANOVA
#aone = ao.anova_oneway(df['Closeness centrality'], df['Depth'],
#                       use_var='equal', welch_correction=False)

#https://pypi.org/project/hypothetical/
#https://aaronschlegel.me/games-howell-post-hoc-multiple-comparisons-test-python.html


aov = hp.aov.AnovaOneWay(df['Closeness centrality'], group=df['Depth'])

print(aov.test_summary)


hp_games = hp.posthoc.GamesHowell(df['Closeness centrality'], group=df['Depth'], alpha=0.05)

print(hp_games.test_result)

hp_tukey = hp.posthoc.TukeysTest(np.asarray(df['Closeness centrality']), group=np.asarray(df['Depth']))

print(hp_tukey.test_result)

'''

'''

#NORMALITY TESTS, INACCURATELY PERFORMED:

if normality:

    print('SHAPIRO')
    
    print('1 Porosity of connected cluster')
    
    print(stats.shapiro(p_conn_1))
    print(stats.shapiro(p_conn_2))
    print(stats.shapiro(p_conn_3))
    
    print('2 Number of pores')
    
    print(stats.shapiro(np_1))
    print(stats.shapiro(np_2))
    print(stats.shapiro(np_3))
    
    print('3 Coordination number')
    
    print(stats.shapiro(cn_1))
    print(stats.shapiro(cn_2))
    print(stats.shapiro(cn_3))
    
    print('4 Clustering coefficient')
    
    print(stats.shapiro(cc_1))
    print(stats.shapiro(cc_2))
    print(stats.shapiro(cc_3))
    
    print('5 Geometric tortuosity')
    
    print(stats.shapiro(gt_1))
    print(stats.shapiro(gt_2))
    print(stats.shapiro(gt_3))
    
    print('6 Closeness centrality')
    
    print(stats.shapiro(ccntr_1))
    print(stats.shapiro(ccntr_2))
    print(stats.shapiro(ccntr_3))
    
    print('7 Betweenness centrality')
    
    print(stats.shapiro(bcntr_1))
    print(stats.shapiro(bcntr_2))
    print(stats.shapiro(bcntr_3))
    
    print('8 Top-down betweenness centrality')
    
    print(stats.shapiro(tdbc_1))
    print(stats.shapiro(tdbc_2))
    print(stats.shapiro(tdbc_3))
    
    print('9 Total porosity')
    
    print(stats.shapiro(pt_1))
    print(stats.shapiro(pt_2))
    print(stats.shapiro(pt_3))
    
    print('10 Volume fraction of connected cluster')
    
    print(stats.shapiro(vfcc_1))
    print(stats.shapiro(vfcc_2))
    print(stats.shapiro(vfcc_3))
    
    
    print('ANDERSON')
    
    print('Porosity of connected cluster')
    
    print(stats.anderson(p_conn_1))
    print(stats.anderson(p_conn_2))
    print(stats.anderson(p_conn_3))
    
    print('Number of pores')
    
    print(stats.anderson(np_1))
    print(stats.anderson(np_2))
    print(stats.anderson(np_3))
    
    print('Coordination number')
    
    print(stats.anderson(cn_1))
    print(stats.anderson(cn_2))
    print(stats.anderson(cn_3))
    
    print('Clustering coefficient')
    
    print(stats.anderson(cc_1))
    print(stats.anderson(cc_2))
    print(stats.anderson(cc_3))
    
    print('Geometric tortuosity')
    
    print(stats.anderson(gt_1))
    print(stats.anderson(gt_2))
    print(stats.anderson(gt_3))
    
    print('Closeness centrality')
    
    print(stats.anderson(ccntr_1))
    print(stats.anderson(ccntr_2))
    print(stats.anderson(ccntr_3))
    
    print('Betweenness centrality')
    
    print(stats.anderson(bcntr_1))
    print(stats.anderson(bcntr_2))
    print(stats.anderson(bcntr_3))
    
    print('Top-down betweenness centrality')
    
    print(stats.anderson(tdbc_1))
    print(stats.anderson(tdbc_2))
    print(stats.anderson(tdbc_3))
    
    print('Total porosity')
    
    print(stats.anderson(pt_1))
    print(stats.anderson(pt_2))
    print(stats.anderson(pt_3))
    
    print('Volume fraction of connected cluster')
    
    print(stats.anderson(vfcc_1))
    print(stats.anderson(vfcc_2))
    print(stats.anderson(vfcc_3))
    
    fignum = 1
    fig = plt.figure(num=fignum)
    fig.set_size_inches(15,5)
    plt.clf()
    ax1 = fig.add_subplot(131)
    stats.probplot(np_1, plot= plt, rvalue= True)
    ax1.set_title(str(fignum))
    ax2 = fig.add_subplot(132)
    stats.probplot(np_2, plot= plt, rvalue= True)
    ax2.set_title(str(fignum))
    ax3 = fig.add_subplot(133)
    stats.probplot(np_3, plot= plt, rvalue= True)
    ax3.set_title(str(fignum))
    
    mm.probability_plot(1, p_conn_1, p_conn_2, p_conn_3)
    
    mm.probability_plot(2, np_1, np_2, np_3)
    
    mm.probability_plot(3, cn_1, cn_2, cn_3)
    
    mm.probability_plot(4, cc_1, cc_2, cc_3)
    
    mm.probability_plot(5, gt_1, gt_2, gt_3)
    
    mm.probability_plot(6, ccntr_1, ccntr_2, ccntr_3)
    
    mm.probability_plot(7, bcntr_1, bcntr_2, bcntr_3)
    
    mm.probability_plot(8, tdbc_1, tdbc_2, tdbc_3)
    
    mm.probability_plot(9, pt_1, pt_2, pt_3)
    
    mm.probability_plot(10, vfcc_1, vfcc_2, vfcc_3)

'''







