# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 15:28:29 2020

@author: pkiuru

Script for calculating statistical tests for pore sizes and measured
volumetric water contents (water retention curves)
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

rcParams['axes.titlesize'] = fontsize + 2
rcParams['axes.labelsize'] = fontsize
rcParams['axes.titleweight'] = 'normal'

rcParams['xtick.labelsize'] = fontsize
rcParams['ytick.labelsize'] = fontsize


datafile = '../Data/Network_metrics.xlsx'
dataread = 1

variance_homogen = 1
anova = 1
tukey = 0
games = 0
dunn = 0
boxpl = 0
linfit = 0

# Basic statistics:
# df['vF'].describe()

if dataread:
    
    #Water retention measurements
    #Original: vwc;      square root-transformed: vwc log10
    df=pd.read_excel(datafile, sheet_name='vwc',
                     header=0)
    
    p_conn_1 = df['vF'][0:7]
    p_conn_2 = df['vF'][7:14]
    p_conn_3 = df['vF'][14:21]
    
    np_1 = df['vwc1kPa'][0:7]
    np_2 = df['vwc1kPa'][7:14]
    np_3 = df['vwc1kPa'][14:21]
    
    cn_1 = df['vwc3kPa'][0:7]
    cn_2 = df['vwc3kPa'][7:14]
    cn_3 = df['vwc3kPa'][14:21]
    
    cc_1 = df['vwc6kPa'][0:7]
    cc_2 = df['vwc6kPa'][7:14]
    cc_3 = df['vwc6kPa'][14:21]
    
    gt_1 = df['vwc10kPa'][0:7]
    gt_2 = df['vwc10kPa'][7:14]
    gt_3 = df['vwc10kPa'][14:21]
    
    ccntr_1 = df['pdiammedian'][0:7]
    ccntr_2 = df['pdiammedian'][7:14]
    ccntr_3 = df['pdiammedian'][14:21]
    
    bcntr_1 = df['pdiammean'][0:7]
    bcntr_2 = df['pdiammean'][7:14]
    bcntr_3 = df['pdiammean'][14:21]
    
    tdbc_1 = df['tdiammean'][0:7]
    tdbc_2 = df['tdiammean'][7:14]
    tdbc_3 = df['tdiammean'][14:21]
    
    pt_1 = df['tdiamintmean'][0:7]
    pt_2 = df['tdiamintmean'][7:14]
    pt_3 = df['tdiamintmean'][14:21]
    
    vfcc_1 = df['rho_b'][0:7]
    vfcc_2 = df['rho_b'][7:14]
    vfcc_3 = df['rho_b'][14:21]

#np_1 = stats.norm.rvs(10, 2, size=1000)
#np_2 = stats.norm.rvs(10, 2, size=1000)
#np_3 = stats.norm.rvs(10, 2, size=1000)
#np_2 = stats.gamma.rvs(0.8, loc=0, scale=1, size=10)


# Homogeneity of variance

if variance_homogen:

    
    print('')
    print('BARTLETT')
    print('')
    
    print('vF')
    
    print(stats.bartlett(p_conn_1, p_conn_2, p_conn_3))
    print([np.var(x, ddof=1) for x in [p_conn_1, p_conn_2, p_conn_3]])
    
    print('vwc1kPa')
    
    print(stats.bartlett(np_1, np_2, np_3))
    print([np.var(x, ddof=1) for x in [np_1, np_2, np_3]])
    
    print('vwc3kPa')
    
    print(stats.bartlett(cn_1, cn_2, cn_3))
    print([np.var(x, ddof=1) for x in [cn_1, cn_2, cn_3]])
    
    print('vwc6kPa')
    
    print(stats.bartlett(cc_1, cc_2, cc_3))
    print([np.var(x, ddof=1) for x in [cc_1, cc_2, cc_3]])
    
    print('vwc10kPa')
    
    print(stats.bartlett(gt_1, gt_2, gt_3))
    print([np.var(x, ddof=1) for x in [gt_1, gt_2, gt_3]])
    
    print('pdiammedian')
    
    print(stats.bartlett(ccntr_1, ccntr_2, ccntr_3))
    print([np.var(x, ddof=1) for x in [ccntr_1, ccntr_2, ccntr_3]])
    
    print('pdiammean')
    
    print(stats.bartlett(bcntr_1, bcntr_2, bcntr_3))
    print([np.var(x, ddof=1) for x in [bcntr_1, bcntr_2, bcntr_3]])
    
    print('tdiammean')
    
    print(stats.bartlett(tdbc_1, tdbc_2, tdbc_3))
    print([np.var(x, ddof=1) for x in [tdbc_1, tdbc_2, tdbc_3]])
    
    print('tdiamintmean')
    
    print(stats.bartlett(pt_1, pt_2, pt_3))
    print([np.var(x, ddof=1) for x in [pt_1, pt_2, pt_3]])
    
    print('rho_b')
    
    print(stats.bartlett(vfcc_1, vfcc_2, vfcc_3)),
    print([np.var(x, ddof=1) for x in [vfcc_1, vfcc_2, vfcc_3]])

    
    print('')
    print('LEVENE')
    print('')
    
    print('vF')
    
    print(stats.levene(p_conn_1, p_conn_2, p_conn_3))
    print('2 & 3')
    print(stats.levene(p_conn_2, p_conn_3))
    print([np.var(x, ddof=1) for x in [p_conn_1, p_conn_2, p_conn_3]])
    
    print('vwc1kPa')
    
    print(stats.levene(np_1, np_2, np_3))
    print('2 & 3')
    print(stats.levene(np_2, np_3))
    print([np.var(x, ddof=1) for x in [np_1, np_2, np_3]])
    
    print('vwc3kPa')
    
    print(stats.levene(cn_1, cn_2, cn_3))
    print('2 & 3')
    print(stats.levene(cn_2, cn_3))
    print([np.var(x, ddof=1) for x in [cn_1, cn_2, cn_3]])
    
    print('vwc6kPa')
    
    print(stats.levene(cc_1, cc_2, cc_3))
    print('2 & 3')
    print(stats.levene(cc_2, cc_3))
    print([np.var(x, ddof=1) for x in [cc_1, cc_2, cc_3]])
    
    print('vwc10kPa')
    
    print(stats.levene(gt_1, gt_2, gt_3))
    print('2 & 3')
    print(stats.levene(gt_2, gt_3))
    print([np.var(x, ddof=1) for x in [gt_1, gt_2, gt_3]])
    
    print('pdiammedian')
    
    print(stats.levene(ccntr_1, ccntr_2, ccntr_3))
    print('2 & 3')
    print(stats.levene(ccntr_2, ccntr_3))
    print([np.var(x, ddof=1) for x in [ccntr_1, ccntr_2, ccntr_3]])
    
    print('pdiammean')
    
    print(stats.levene(bcntr_1, bcntr_2, bcntr_3))
    print('2 & 3')
    print(stats.levene(bcntr_2, bcntr_3))
    print([np.var(x, ddof=1) for x in [bcntr_1, bcntr_2, bcntr_3]])
    
    print('tdiammean')
    
    print(stats.levene(tdbc_1, tdbc_2, tdbc_3))
    print('2 & 3')
    print(stats.levene(tdbc_2, tdbc_3))
    print([np.var(x, ddof=1) for x in [tdbc_1, tdbc_2, tdbc_3]])
    
    print('tdiamintmean')
    
    print(stats.levene(pt_1, pt_2, pt_3))
    print('2 & 3')
    print(stats.levene(pt_2, pt_3))
    print([np.var(x, ddof=1) for x in [pt_1, pt_2, pt_3]])
    
    print('rho_b')
    
    print(stats.levene(vfcc_1, vfcc_2, vfcc_3))
    print('2 & 3')
    print(stats.levene(vfcc_2, vfcc_3))
    print([np.var(x, ddof=1) for x in [vfcc_1, vfcc_2, vfcc_3]])
# ANOVA


if anova:
    
    print('')
    print('Standard ANOVA, F-ONEWAY')
    print('')
    
    print('vF')
    
    print(stats.f_oneway(p_conn_1, p_conn_2, p_conn_3))
    
    print('vwc1kPa')
    
    print(stats.f_oneway(np_1, np_2, np_3))
    
    print('vwc3kPa')
    
    print(stats.f_oneway(cn_1, cn_2, cn_3))
    
    print('vwc6kPa')
    
    print(stats.f_oneway(cc_1, cc_2, cc_3))
    
    print('vwc10kPa')
    
    print(stats.f_oneway(gt_1, gt_2, gt_3))
    
    print('pdiaammedian')
    
    print(stats.f_oneway(ccntr_1, ccntr_2, ccntr_3))
    
    print('pdiammean')
    
    print(stats.f_oneway(bcntr_1, bcntr_2, bcntr_3))
    
    print('tdiammean')
    
    print(stats.f_oneway(tdbc_1, tdbc_2, tdbc_3))
    
    print('tdiamintmean')
    
    print(stats.f_oneway(pt_1, pt_2, pt_3))
    
    print('rho_b')
    
    print(stats.f_oneway(vfcc_1, vfcc_2, vfcc_3))
    
    
    
    if True:
        print('')
        print('Standard ANOVA ols fit')
        print('')
        
        
        st_p_conn = ols('vF ~ C(Depth)', data=df).fit()
        print('-------- vF ----------')
        print(sm.stats.anova_lm(st_p_conn, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_p_conn.resid)))
        print(stats.anderson(st_p_conn.resid))
        mm.probability_plot_ols(1, st_p_conn.resid)
        
        st_np = ols('vwc1kPa ~ C(Depth)', data=df).fit()
        print('-------- vwc1kPa ----------')
        print(sm.stats.anova_lm(st_np, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_np.resid)))
        print(stats.anderson(st_np.resid))
        mm.probability_plot_ols(2, st_np.resid)
        
        st_cn = ols('vwc3kPa ~ C(Depth)', data=df).fit()
        print('-------- vwc3kPa ----------')
        print(sm.stats.anova_lm(st_cn, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_cn.resid)))
        print(stats.anderson(st_cn.resid))
        mm.probability_plot_ols(3, st_cn.resid)
        
        st_cc = ols('vwc6kPa ~ C(Depth)', data=df).fit()
        print('-------- vwc6kPa ----------')
        print(sm.stats.anova_lm(st_cc, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_cc.resid)))
        print(stats.anderson(st_cc.resid))
        mm.probability_plot_ols(4, st_cc.resid)
        
        st_gt = ols('vwc10kPa ~ C(Depth)', data=df).fit()
        print('-------- vwc10kPa ----------')
        print(sm.stats.anova_lm(st_gt, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_gt.resid)))
        print(stats.anderson(st_gt.resid))
        mm.probability_plot_ols(5, st_gt.resid)
        
        st_ccntr = ols('pdiammedian ~ C(Depth)', data=df).fit()
        print('-------- pdiammedian ----------')
        print(sm.stats.anova_lm(st_ccntr, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_ccntr.resid)))
        print(stats.anderson(st_ccntr.resid))
        mm.probability_plot_ols(6, st_ccntr.resid)
        
        
        st_bcntr = ols('pdiammean ~ C(Depth)', data=df).fit()
        print('-------- pdiammean ----------')
        print(sm.stats.anova_lm(st_bcntr, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_bcntr.resid)))
        print(stats.anderson(st_bcntr.resid))
        mm.probability_plot_ols(7, st_bcntr.resid)
        
        st_tdbc = ols('tdiammean ~ C(Depth)', data=df).fit()
        print('-------- tdiammean ----------')
        print(sm.stats.anova_lm(st_tdbc, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_tdbc.resid)))
        print(stats.anderson(st_tdbc.resid))
        mm.probability_plot_ols(8, st_tdbc.resid)
        
        st_pt = ols('tdiamintmean ~ C(Depth)', data=df).fit()
        print('-------- tdiamintmean ----------')
        print(sm.stats.anova_lm(st_pt, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_pt.resid)))
        print(stats.anderson(st_pt.resid))
        mm.probability_plot_ols(9, st_pt.resid)
        
        st_vfcc = ols('rho_b ~ C(Depth)', data=df).fit()
        print('-------- rho_b ----------')
        print(sm.stats.anova_lm(st_vfcc, typ=2))
        print('Residual Shapiro ' + str(stats.shapiro(st_vfcc.resid)))
        print(stats.anderson(st_vfcc.resid))
        mm.probability_plot_ols(10, st_vfcc.resid)
        
    
    #except Exception:
    
    #    pass
    
    
    try:
        print('')
        print('Welch ANOVA')
        print('')
        
        welch = ao.anova_oneway(df['vF'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- vF ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['vwc1kPa'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- vwc1kPa ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['vwc3kPa'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- vwc3kPa ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['vwc6kPa'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- vwc6kPa ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['vwc10kPa'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- vwc10kPa ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['pdiammedian'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- pdiammedian ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['pdiammean'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- pdiammean ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['tdiammean'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- tdiammean ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['tdiamintmean'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- tdiamintmean ----------')
        print(welch)
        
        welch = ao.anova_oneway(df['rho_b'], df['Depth'],
                           use_var='unequal', welch_correction=True)
        print('-------- rho_b ----------')
        print(welch)
    
    except Exception:
    
        pass
    
    
    print('')
    print('KRUSKAL')
    print('')
    
    print('vF')
    
    print(stats.kruskal(p_conn_1, p_conn_2, p_conn_3))
    
    print('vwc1kPa')
    
    print(stats.kruskal(np_1, np_2, np_3))
    
    print('vwc3kPa')
    
    print(stats.kruskal(cn_1, cn_2, cn_3))
    
    print('vwc6kPa')
    
    print(stats.kruskal(cc_1, cc_2, cc_3))
    
    print('vwc10kPa')
    
    print(stats.kruskal(gt_1, gt_2, gt_3))
    
    print('pdiammedian')
    
    print(stats.kruskal(ccntr_1, ccntr_2, ccntr_3))
    
    print('pdiammean')
    
    print(stats.kruskal(bcntr_1, bcntr_2, bcntr_3))
    
    print('tdiammean')
    
    print(stats.kruskal(tdbc_1, tdbc_2, tdbc_3))
    
    print('tdiamintmean')
    
    print(stats.kruskal(pt_1, pt_2, pt_3))
    
    print('rho_b')
    
    print(stats.kruskal(vfcc_1, vfcc_2, vfcc_3))

if tukey:
    
    print('Tukey')
    print('')
    
    comp = mc.MultiComparison(df['vF'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('vF')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['vwc1kPa'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('vwc1kPa')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['vwc3kPa'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('vwc3kPa')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['vwc6kPa'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('vwc6kPa')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['vwc10kPa'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('vwc10kPa')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['pdiammedian'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('pdiammedian')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['pdiammean'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('pdiammean')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['tdiammean'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('tdiammean')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['tdiamintmean'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('tdiamintmean')
    print(post_hoc_res.summary())
    
    comp = mc.MultiComparison(df['rho_b'], df['Depth'])
    post_hoc_res = comp.tukeyhsd()
    print('rho_b')
    print(post_hoc_res.summary())
    

if games:
    
    print('Games')
    print('')
    
    hp_games = hp.posthoc.GamesHowell(df['vF'], group=df['Depth'], alpha=0.05)
    print('vF')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['vwc1kPa'], group=df['Depth'], alpha=0.05)
    print('vwc1kPa')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['vwc3kPa'], group=df['Depth'], alpha=0.05)
    print('vwc3kPa')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['vwc6kPa'], group=df['Depth'], alpha=0.05)
    print('vwc6kPa')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['vwc10kPa'], group=df['Depth'], alpha=0.05)
    print('vwc10kPa')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['pdiammedian'], group=df['Depth'], alpha=0.05)
    print('pdiammedian')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['pdiammean'], group=df['Depth'], alpha=0.05)
    print('pdiammean')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['tdiammean'], group=df['Depth'], alpha=0.05)
    print('tdiammean')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['tdiamintmean'], group=df['Depth'], alpha=0.05)
    print('tdiamintmean')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
    hp_games = hp.posthoc.GamesHowell(df['rho_b'], group=df['Depth'], alpha=0.05)
    print('rho_b')
    print(hp_games.test_result)
    print(hp_games.test_result['p_value'])
    
if dunn:
    
    p_adjust =  'holm'# 'bonferroni' 'sidak'
    
    
    
    print('')
    print('Dunn')
    print('')
    print('p-adjust is now ')
    print(p_adjust)
    print('')
    
    print('vF')
    ph_dunn = sp.posthoc_dunn(df, val_col='vF', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('vwc1kPa')
    ph_dunn = sp.posthoc_dunn(df, val_col='vwc1kPa', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('vwc3kPa')
    ph_dunn = sp.posthoc_dunn(df, val_col='vwc3kPa', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('vwc6kPa')
    ph_dunn = sp.posthoc_dunn(df, val_col='vwc6kPa', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('vwc10kPa')
    ph_dunn = sp.posthoc_dunn(df, val_col='vwc10kPa', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('pdiammedian')
    ph_dunn = sp.posthoc_dunn(df, val_col='pdiammedian', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('pdiammean')
    ph_dunn = sp.posthoc_dunn(df, val_col='pdiammean', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('tdiammean')
    ph_dunn = sp.posthoc_dunn(df, val_col='tdiammean', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('tdiamintmean')
    ph_dunn = sp.posthoc_dunn(df, val_col='tdiamintmean', group_col='Depth', p_adjust = p_adjust)    
    print(ph_dunn)
    
    print('rho_b')
    ph_dunn = sp.posthoc_dunn(df, val_col='rho_b', group_col='Depth', p_adjust = p_adjust)    
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
    bbdim = [0.17, 0.37]
    
    bppos = [
            [0.07, 0.56, bbdim[0],bbdim[1]],
            [0.31, 0.56, bbdim[0],bbdim[1]],
            [0.55, 0.56, bbdim[0],bbdim[1]],
            [0.79, 0.56, bbdim[0],bbdim[1]],
            [0.07, 0.08, bbdim[0],bbdim[1]],
            [0.31, 0.08, bbdim[0],bbdim[1]],
            [0.55, 0.08, bbdim[0],bbdim[1]],
            [0.79, 0.08, bbdim[0],bbdim[1]],
            ]
    
    subt = [('a)          '), ('b)          '), ('c)          '), ('d)          '),
            ('e)          '), ('f)          '), ('g)          '), ('h)          ')]
    
    ylabs = [
            'Network porosity',
            'vwc1kPa (\u00D710$^3$)',
            'vwc3kPa',
            'vwc6kPa',
            'vwc10kPa',
            'pdiammedian',
            'pdiammean',
            'tdiammean'
            ]
    
    fig = plt.figure(num=100)
    fig.set_size_inches(14,9)
    plt.clf()
    
    for i in range(2,10):
         
        ax = fig.add_subplot(2,4,i-1)
        
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
                   showmeans=False)
        
        if i-2 > -2:
            ax.text(0.13, 0.9, "A", transform=ax.transAxes)
            ax.text(0.47, 0.9, "B", transform=ax.transAxes)
            ax.text(0.8, 0.9, "B", transform=ax.transAxes)
        
        ax.set_title(subt[i-2], horizontalalignment='right')
        
        #ax.set_ylim([0, 1.2 * factor * np.max(df[df.columns[i]])])
        
        if i-1 > 4:
            ax.set_xlabel("Depth (cm)")
        ax.set_ylabel(ylabs[i-2], labelpad=0.9)
        
        ax.set_position(bppos[i-2])
    
if linfit:
    
    bbdim = [0.38, 0.18]
    
    ds = ['0-5 cm', '20-25 cm', '40-45 cm']
    
    bppos = [
            [0.08, 0.76, bbdim[0],bbdim[1]],
            [0.55, 0.76, bbdim[0],bbdim[1]],
            [0.08, 0.52, bbdim[0],bbdim[1]],
            [0.55, 0.52, bbdim[0],bbdim[1]],
            [0.08, 0.28, bbdim[0],bbdim[1]],
            [0.55, 0.28, bbdim[0],bbdim[1]],
            [0.08, 0.04, bbdim[0],bbdim[1]],
            [0.55, 0.04, bbdim[0],bbdim[1]],
            ]
    
    subt = [('a)          '), ('b)          '), ('c)          '), ('d)          '),
            ('e)          '), ('f)          '), ('g)          '), ('h)          ')]
    
    ylabs = [
            'vF',
            'vwc1kPa (\u00D710$^3$)',
            'vwc3kPa',
            'vwc6kPa',
            'vwc10kPa',
            'vwc10kPa',
            'vwc10kPa (\u00D710$^{-3}$)',
            'vWc10kPa (\u00D710$^{-3}$)'
            ]
    
    lineqs = {
            "P_vs_Np_1": "y = 4.1e-5x - 0.068, R$^2$ = 0.75",
            "P_vs_Np_2": "y = 1.8e-5x - 0.002, R$^2$ = 0.91",
            "P_vs_Np_3": "y = 2.1e-5x - 0.016, R$^2$ = 0.94"
            }
    
    fig = plt.figure(num=200)
    fig.set_size_inches(7,9)
    plt.clf()
    
    ax1 = fig.add_subplot(4,2,1)
    
    dummyx = df['vwc1kPa'][df['Depth'] == 1]
    dummyy = df['vF'][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax1.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax1.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax1.text(0.07, 0.9, ds[0] + ': ' + lineqs['P_vs_Np_1'] + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax1.transAxes, fontsize=fontsize)
    
    dummyx = df['vwc1kPa'][df['Depth'] == 2]
    dummyy = df['vF'][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax1.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax1.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax1.text(0.07, 0.8, ds[1] + ': ' + lineqs['P_vs_Np_2'] + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax1.transAxes, fontsize=fontsize)
    
    dummyx = df['vwc1kPa'][df['Depth'] == 3]
    dummyy = df['vF'][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax1.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax1.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax1.text(0.07, 0.7, ds[2] + ': ' + lineqs['P_vs_Np_3'] + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax1.transAxes, fontsize=fontsize)
    
    ax1.set_xlim([0, 14000])
    ax1.set_ylim([0.0, 0.9])
    
    ax1.text(0.02, 1.05, subt[0], transform=ax1.transAxes)
    
    #handles, labels = ax1.get_legend_handles_labels()
    #order = [3,4,5]
    #ax1.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax1.legend(fontsize=fontsize)
    
    ax1.set_xlabel('vwc1kPa')
    ax1.set_ylabel('Network porosity')
    
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    ax1.set_position(bppos[0])   
    
    #----------
    
    ax2 = fig.add_subplot(4,2,2)
    
    xcol = 'vF'
    ycol = 'vwc3kPa'
    
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax2.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax2.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax2.text(0.07, 0.9, ds[0] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax2.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax2.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax2.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax2.text(0.07, 0.8, ds[1] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax2.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax2.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax2.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax2.text(0.07, 0.7, ds[2] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax2.transAxes, fontsize=fontsize)
    
    ax2.set_xlim([0, 0.7])
    ax2.set_ylim([2, 12])
    
    ax2.text(0.02, 1.05, subt[1], transform=ax2.transAxes)
    
    #handles, labels = ax2.get_legend_handles_labels()
    #order = [3,4,5]
    #ax2.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax2.legend(fontsize=fontsize)
    
    ax2.set_xlabel('Network porosity')
    ax2.set_ylabel('vwc3kPa')
    
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    ax2.set_position(bppos[1])    
    
    #-------------------
    
    ax3 = fig.add_subplot(4,2,3)
    
    xcol = 'vwc3kPa'
    ycol = 'vwc6kPa'
    
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax3.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax3.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax3.text(0.07, 0.9, ds[0] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax3.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax3.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax3.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax3.text(0.07, 0.8, ds[1] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax3.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax3.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax3.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax3.text(0.07, 0.7, ds[2] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax3.transAxes, fontsize=fontsize)
    
    ax3.set_xlim([2, 10])
    ax3.set_ylim([0, 0.9])
    
    ax3.text(0.02, 1.05, subt[2], transform=ax3.transAxes)
    
    #handles, labels = ax3.get_legend_handles_labels()
    #order = [3,4,5]
    #ax3.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax3.legend(fontsize=fontsize)
    
    ax3.set_xlabel('vwc3kPa')
    ax3.set_ylabel('vwc6kPa')
    
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    
    ax3.set_position(bppos[2])
    
    #----------
    
    ax4 = fig.add_subplot(4,2,4)
    
    xcol = 'vF'
    ycol = 'vwc10kPa'
    
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax4.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax4.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax4.text(0.07, 0.9, ds[0] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax4.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][np.array([0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0]).astype('bool')]
    dummyy = df[ycol][df['Depth'] == 2].dropna()
    lr_1 = stats.linregress(dummyx,dummyy)
    ax4.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax4.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax4.text(0.07, 0.8, ds[1] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax4.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax4.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax4.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax4.text(0.07, 0.7, ds[2] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax4.transAxes, fontsize=fontsize)
    
    ax4.set_xlim([0, 0.7])
    ax4.set_ylim([1, 6.0])
    
    ax4.text(0.02, 1.05, subt[3], transform=ax4.transAxes)
    
    #handles, labels = ax4.get_legend_handles_labels()
    #order = [3,4,5]
    #ax4.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax4.legend(fontsize=fontsize)
    
    ax4.set_xlabel('Network porosity')
    ax4.set_ylabel('vwc10kPa')
    
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    
    ax4.set_position(bppos[3])    
    
    #----------
    
    ax5 = fig.add_subplot(4,2,5)
    
    xcol = 'vwc3kPa'
    ycol = 'vwc10kPa'
    
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax5.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax5.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax5.text(0.07, 0.9, ds[0] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax5.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax5.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax5.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax5.text(0.07, 0.8, ds[1] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax5.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax5.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax5.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax5.text(0.07, 0.7, ds[2] + ': y = ' + str(np.round(lr_1.slope,2)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax5.transAxes, fontsize=fontsize)
    
    ax5.set_xlim([2, 10])
    ax5.set_ylim([15.0, 60.0])
    
    ax5.text(0.02, 1.05, subt[4], transform=ax5.transAxes)
    
    #handles, labels = ax5.get_legend_handles_labels()
    #order = [3,4,5]
    #ax5.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, loc='upper right')
    
    #ax5.legend(fontsize=fontsize)
    
    ax5.set_xlabel('vwc3kPa')
    ax5.set_ylabel('vwc10kPa')
    
    ax5.spines["top"].set_visible(False)
    ax5.spines["right"].set_visible(False)
    
    ax5.set_position(bppos[4])    
    
    #----------
    
    ax6 = fig.add_subplot(4,2,6)
    
    xcol = 'vwc3kPa'
    ycol = 'vwc10kPa'
    
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax6.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax6.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax6.text(0.07, 0.9, ds[0] + ': y = ' + str(np.round(lr_1.slope,4)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax6.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 2]
    dummyy = df[ycol][df['Depth'] == 2]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax6.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax6.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax6.text(0.07, 0.8, ds[1] + ': y = ' + str(np.round(lr_1.slope,4)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax6.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax6.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax6.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax6.text(0.07, 0.7, ds[2] + ': y = ' + str(np.round(lr_1.slope,4)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax6.transAxes, fontsize=fontsize)
    
    ax6.set_xlim([2, 10])
    ax6.set_ylim([0.0, 0.019])
    
    handles, labels = ax6.get_legend_handles_labels()
    order = [3,4,5]
    ax6.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
               fontsize=fontsize,  bbox_to_anchor=(0.68, 0.22),
               bbox_transform=plt.gcf().transFigure)
    
    #ax6.legend(fontsize=fontsize)
    
    ax6.set_xlabel('vwc3kPa')
    ax6.set_ylabel('vwc10kPa')
    
    ax6.text(0.02, 1.05, subt[5], transform=ax6.transAxes)
    
    ax6.spines["top"].set_visible(False)
    ax6.spines["right"].set_visible(False)
    
    ax6.set_position(bppos[5])    
    
    #----------
    
    ax7 = fig.add_subplot(4,2,7)
    
    xcol = 'vwc10kPa'
    ycol = 'vwc10kPa'
    
    dummyx = df[xcol][df['Depth'] == 1]
    dummyy = df[ycol][df['Depth'] == 1]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax7.scatter(dummyx, dummyy, marker='o', c='k', s=6, label=ds[0])
    ax7.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'k-', lw=1, label='fit')
    ax7.text(0.07, 0.9, ds[0] + ': y = ' + str(np.round(lr_1.slope,4)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax7.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 2].dropna()
    dummyy = df[ycol][df['Depth'] == 2].dropna()
    lr_1 = stats.linregress(dummyx,dummyy)
    ax7.scatter(dummyx, dummyy, marker='s', c='b', s=6, label=ds[1])
    ax7.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'b--', lw=1, label='fit')
    ax7.text(0.07, 0.8, ds[1] + ': y = ' + str(np.round(lr_1.slope,4)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax7.transAxes, fontsize=fontsize)
    
    dummyx = df[xcol][df['Depth'] == 3]
    dummyy = df[ycol][df['Depth'] == 3]
    lr_1 = stats.linregress(dummyx,dummyy)
    ax7.scatter(dummyx, dummyy, marker='^', c='r', s=6, label=ds[2])
    ax7.plot(dummyx, lr_1.slope * dummyx + lr_1.intercept, 'r-.', lw=1, label='fit')
    ax7.text(0.07, 0.7, ds[2] + ': y = ' + str(np.round(lr_1.slope,4)) + "x + "
             + str(np.round(lr_1.intercept,3))
             + ", R$^2$ = " + str(np.round(lr_1.rvalue**2,2)) + ' ' + str(np.round(lr_1.pvalue,3)),
             transform=ax7.transAxes, fontsize=fontsize)
    
    ax7.set_xlim([1, 4])
    ax7.set_ylim([0.0, 0.025])
    
    ax7.text(0.02, 1.05, subt[6], transform=ax7.transAxes)
    
    #handles, labels = ax7.get_legend_handles_labels()
    #order = [3,4,5]
    #ax7.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    #           fontsize=fontsize, bbox_to_anchor=(0.8, 0.1),
    #           bbox_transform=plt.gcf().transFigure)
    
    #ax7.legend(fontsize=fontsize)
    
    ax7.set_xlabel('vwc10kPa')
    ax7.set_ylabel('Top-bottom betw. centr.')
    
    ax7.spines["top"].set_visible(False)
    ax7.spines["right"].set_visible(False)
    
    ax7.set_position(bppos[6])    
    
    #----------

if True:
    vals = []
    for h in range(0,1000):
    
        
        
        t_1 = stats.norm.rvs(4, 1, size=500)
        t_2 = stats.norm.rvs(8, 1.15, size=500)
        t_3 = stats.norm.rvs(6, 1.06, size=500)
        
        st = stats.levene(t_1, t_2, t_3, center='mean')[1]
        vals.append(st)
    
        #print(stats.levene(t_1, t_2, t_3, center='mean'))
        
    print(np.mean(np.asarray(vals)))        
    
    #print(stats.levene(t_1, t_2, t_3, center='median'))

'''


comp = mc.MultiComparison(df['vwc10kPa'], df['Depth'])
post_hoc_res = comp.tukeyhsd()
print(post_hoc_res.summary())


#Welch ANOVA https://statisticsbyjim.com/anova/welchs-anova-compared-to-classic-one-way-anova/
#https://www.statsmodels.org/devel/generated/statsmodels.stats.oneway.anova_oneway.html#statsmodels.stats.oneway.anova_oneway
aone = ao.anova_oneway(df['vwc10kPa'], df['Depth'],
                       use_var='unequal', welch_correction=True)
print(aone)
#Standard ANOVA
#aone = ao.anova_oneway(df['vwc10kPa'], df['Depth'],
#                       use_var='equal', welch_correction=False)

#https://pypi.org/project/hypothetical/
#https://aaronschlegel.me/games-howell-post-hoc-multiple-comparisons-test-python.html


aov = hp.aov.AnovaOneWay(df['vwc10kPa'], group=df['Depth'])

print(aov.test_summary)


hp_games = hp.posthoc.GamesHowell(df['vwc10kPa'], group=df['Depth'], alpha=0.05)

print(hp_games.test_result)

hp_tukey = hp.posthoc.TukeysTest(np.asarray(df['vwc10kPa']), group=np.asarray(df['Depth']))

print(hp_tukey.test_result)

'''
