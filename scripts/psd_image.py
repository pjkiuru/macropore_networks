# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 18:13:53 2020

@author: pkiuru
"""

#import openpnm as op
#import porespy as ps
import numpy as np
import matplotlib.pyplot as plt
import more_models as mm
import scipy.stats as stats
from matplotlib.ticker import FuncFormatter#, FormatStrFormatter, MaxNLocator
#from scipy.interpolate import interp1d

from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in' 

fontname = 'Arial'
fontsize = 8

rcParams['font.size'] = fontsize

rcParams['axes.titlesize'] = fontsize
rcParams['axes.labelsize'] = fontsize
rcParams['axes.titleweight'] = 'normal'

rcParams['xtick.labelsize'] = fontsize
rcParams['ytick.labelsize'] = fontsize

rcParams['legend.handlelength'] = 1.0

nbins = 20
volfact = 1e9


a = '20_25'

loading = 0

# Gamma fits for pore size distribution of each network
regular = 0

#Histograms of combined pore sizes of all samples from a single depth
joint = 1

#Calculation of mean and median values of pore and throat sizes
means = 0

#Goodness-of-fit calculations for distribution fits
ppr2 = 0

# -------------------------- LOADING AND CALCULATION --------------------------------

if loading:
    
    a = '0_5'
    
    vol_1_1, vol_all_1_1, int_1_1, int_all_1_1, thr_1_1, thr_int_1_1, thr_all_1_1 = mm.load_network_volumes(a+'_1')
    vol_1_2, vol_all_1_2, int_1_2, int_all_1_2, thr_1_2, thr_int_1_2, thr_all_1_2 = mm.load_network_volumes(a+'_2')
    if a == '40_45':
        vol_1_3, vol_all_1_3, int_1_3, int_all_1_3, thr_1_3, thr_int_1_3, thr_all_1_3 = mm.load_network_volumes(a+'_1_2')
    else:
        vol_1_3, vol_all_1_3, int_1_3, int_all_1_3, thr_1_3, thr_int_1_3, thr_all_1_3 = mm.load_network_volumes(a+'_3')
    vol_1_4, vol_all_1_4, int_1_4, int_all_1_4, thr_1_4, thr_int_1_4, thr_all_1_4 = mm.load_network_volumes(a+'_4')
    if a == '40_45':
        vol_1_5, vol_all_1_5, int_1_5, int_all_1_5, thr_1_5, thr_int_1_5, thr_all_1_5 = mm.load_network_volumes(a+'_5_2')
    else:
        vol_1_5, vol_all_1_5, int_1_5, int_all_1_5, thr_1_5, thr_int_1_5, thr_all_1_5 = mm.load_network_volumes(a+'_5')
    vol_1_6, vol_all_1_6, int_1_6, int_all_1_6, thr_1_6, thr_int_1_6, thr_all_1_6 = mm.load_network_volumes(a+'_6')
    vol_1_7, vol_all_1_7, int_1_7, int_all_1_7, thr_1_7, thr_int_1_7, thr_all_1_7 = mm.load_network_volumes(a+'_7')
    
    diam_1_1 = mm.equiv_diam(volfact, vol_1_1)
    diam_1_2 = mm.equiv_diam(volfact, vol_1_2)
    diam_1_3 = mm.equiv_diam(volfact, vol_1_3)
    diam_1_4 = mm.equiv_diam(volfact, vol_1_4)
    diam_1_5 = mm.equiv_diam(volfact, vol_1_5)
    diam_1_6 = mm.equiv_diam(volfact, vol_1_6)
    diam_1_7 = mm.equiv_diam(volfact, vol_1_7)
    
    diamfit_s_1_1 = stats.lognorm.fit(diam_1_1[int_1_1])
    diamfit_s_1_2 = stats.lognorm.fit(diam_1_2[int_1_2])
    diamfit_s_1_3 = stats.lognorm.fit(diam_1_3[int_1_3])
    diamfit_s_1_4 = stats.lognorm.fit(diam_1_4[int_1_4])
    diamfit_s_1_5 = stats.lognorm.fit(diam_1_5[int_1_5])
    diamfit_s_1_6 = stats.lognorm.fit(diam_1_6[int_1_6])
    diamfit_s_1_7 = stats.lognorm.fit(diam_1_7[int_1_7])
    
    diam_1_all = np.concatenate([diam_1_1[int_1_1], diam_1_2[int_1_2],
                                diam_1_3[int_1_3], diam_1_4[int_1_4],
                                diam_1_5[int_1_5], diam_1_6[int_1_6],
                                diam_1_7[int_1_7]])
    
    diamfit_s_1 = stats.lognorm.fit(diam_1_all)
    
    
    a = '20_25'
    
    vol_2_1, vol_all_2_1, int_2_1, int_all_2_1, thr_2_1, thr_int_2_1, thr_all_2_1 = mm.load_network_volumes(a+'_1')
    vol_2_2, vol_all_2_2, int_2_2, int_all_2_2, thr_2_2, thr_int_2_2, thr_all_2_2 = mm.load_network_volumes(a+'_2')
    if a == '40_45':
        vol_2_3, vol_all_2_3, int_2_3, int_all_2_3, thr_2_3, thr_int_2_3, thr_all_2_3 = mm.load_network_volumes(a+'_2_2')
    else:
        vol_2_3, vol_all_2_3, int_2_3, int_all_2_3, thr_2_3, thr_int_2_3, thr_all_2_3 = mm.load_network_volumes(a+'_3')
    vol_2_4, vol_all_2_4, int_2_4, int_all_2_4, thr_2_4, thr_int_2_4, thr_all_2_4 = mm.load_network_volumes(a+'_4')
    if a == '40_45':
        vol_2_5, vol_all_2_5, int_2_5, int_all_2_5, thr_2_5, thr_int_2_5, thr_all_2_5 = mm.load_network_volumes(a+'_5_2')
    else:
        vol_2_5, vol_all_2_5, int_2_5, int_all_2_5, thr_2_5, thr_int_2_5, thr_all_2_5 = mm.load_network_volumes(a+'_5')
    vol_2_6, vol_all_2_6, int_2_6, int_all_2_6, thr_2_6, thr_int_2_6, thr_all_2_6 = mm.load_network_volumes(a+'_6')
    vol_2_7, vol_all_2_7, int_2_7, int_all_2_7, thr_2_7, thr_int_2_7, thr_all_2_7 = mm.load_network_volumes(a+'_7')
    
    diam_2_1 = mm.equiv_diam(volfact, vol_2_1)
    diam_2_2 = mm.equiv_diam(volfact, vol_2_2)
    diam_2_3 = mm.equiv_diam(volfact, vol_2_3)
    diam_2_4 = mm.equiv_diam(volfact, vol_2_4)
    diam_2_5 = mm.equiv_diam(volfact, vol_2_5)
    diam_2_6 = mm.equiv_diam(volfact, vol_2_6)
    diam_2_7 = mm.equiv_diam(volfact, vol_2_7)
    
    #https://stackoverflow.com/questions/8747761/scipy-lognormal-distribution-parameters?rq=1
    
    diamfit_s_2_1 = stats.lognorm.fit(diam_2_1[int_2_1])
    diamfit_s_2_2 = stats.lognorm.fit(diam_2_2[int_2_2])
    diamfit_s_2_3 = stats.lognorm.fit(diam_2_3[int_2_3])
    diamfit_s_2_4 = stats.lognorm.fit(diam_2_4[int_2_4])
    diamfit_s_2_5 = stats.lognorm.fit(diam_2_5[int_2_5])
    diamfit_s_2_6 = stats.lognorm.fit(diam_2_6[int_2_6])
    diamfit_s_2_7 = stats.lognorm.fit(diam_2_7[int_2_7])
    
    diam_2_all = np.concatenate([diam_2_1[int_2_1], diam_2_2[int_2_2],
                                diam_2_3[int_2_3], diam_2_4[int_2_4],
                                diam_2_5[int_2_5], diam_2_6[int_2_6],
                                diam_2_7[int_2_7]])
    
    diamfit_s_2 = stats.lognorm.fit(diam_2_all)
    
    a = '40_45'
    
    vol_3_1, vol_all_3_1, int_3_1, int_all_3_1, thr_3_1, thr_int_3_1, thr_all_3_1 = mm.load_network_volumes(a+'_1')
    vol_3_2, vol_all_3_2, int_3_2, int_all_3_2, thr_3_2, thr_int_3_2, thr_all_3_2 = mm.load_network_volumes(a+'_2')
    if a == '40_45':
        vol_3_3, vol_all_3_3, int_3_3, int_all_3_3, thr_3_3, thr_int_3_3, thr_all_3_3 = mm.load_network_volumes(a+'_3_2')
    else:
        vol_3_3, vol_all_3_3, int_3_3, int_all_3_3, thr_3_3, thr_int_3_3, thr_all_3_3 = mm.load_network_volumes(a+'_3')
    vol_3_4, vol_all_3_4, int_3_4, int_all_3_4, thr_3_4, thr_int_3_4, thr_all_3_4 = mm.load_network_volumes(a+'_4')
    if a == '40_45':
        vol_3_5, vol_all_3_5, int_3_5, int_all_3_5, thr_3_5, thr_int_3_5, thr_all_3_5 = mm.load_network_volumes(a+'_5_2')
    else:
        vol_3_5, vol_all_3_5, int_3_5, int_all_3_5, thr_3_5, thr_int_3_5, thr_all_3_5 = mm.load_network_volumes(a+'_5')
    vol_3_6, vol_all_3_6, int_3_6, int_all_3_6, thr_3_6, thr_int_3_6, thr_all_3_6 = mm.load_network_volumes(a+'_6')
    vol_3_7, vol_all_3_7, int_3_7, int_all_3_7, thr_3_7, thr_int_3_7, thr_all_3_7 = mm.load_network_volumes(a+'_7')
    
    diam_3_1 = mm.equiv_diam(volfact, vol_3_1)
    diam_3_2 = mm.equiv_diam(volfact, vol_3_2)
    diam_3_3 = mm.equiv_diam(volfact, vol_3_3)
    diam_3_4 = mm.equiv_diam(volfact, vol_3_4)
    diam_3_5 = mm.equiv_diam(volfact, vol_3_5)
    diam_3_6 = mm.equiv_diam(volfact, vol_3_6)
    diam_3_7 = mm.equiv_diam(volfact, vol_3_7)
    
    diamfit_s_3_1 = stats.lognorm.fit(diam_3_1[int_3_1])
    diamfit_s_3_2 = stats.lognorm.fit(diam_3_2[int_3_2])
    diamfit_s_3_3 = stats.lognorm.fit(diam_3_3[int_3_3])
    diamfit_s_3_4 = stats.lognorm.fit(diam_3_4[int_3_4])
    diamfit_s_3_5 = stats.lognorm.fit(diam_3_5[int_3_5])
    diamfit_s_3_6 = stats.lognorm.fit(diam_3_6[int_3_6])
    diamfit_s_3_7 = stats.lognorm.fit(diam_3_7[int_3_7])
    
    
    diam_3_all = np.concatenate([diam_3_1[int_3_1], diam_3_2[int_3_2],
                                diam_3_3[int_3_3], diam_3_4[int_3_4],
                                diam_3_5[int_3_5], diam_3_6[int_3_6],
                                diam_3_7[int_3_7]])
    
    diamfit_s_3 = stats.lognorm.fit(diam_3_all)
    
    thr_1 = np.concatenate([thr_1_1, thr_1_2, thr_1_3, thr_1_4, thr_1_5,
                           thr_1_6, thr_1_7])
    
    thr_2 = np.concatenate([thr_2_1, thr_2_2, thr_2_3, thr_2_4, thr_2_5,
                           thr_2_6, thr_2_7])
    
    thr_3 = np.concatenate([thr_3_1, thr_3_2, thr_3_3, thr_3_4, thr_3_5,
                           thr_3_6, thr_3_7])



# OOOOOOOOOOOOOOOOOOOOOOOOOO PORE SIZE DISTRIBUTIONS OOOOOOOOOOOOOOOOOOOOOOOOOO

#'''

#https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting/33213196

xmin = 0.1
xmax = 6

xminc = 0.08
xmaxc = 4.8

ymina = -0.05
ymaxa = 1.95

yminb = -0.1
ymaxb = 3.9
    
yminc = -0.05
ymaxc = 1.05

lw = 1

spdim = [0.24, 0.25]

sppos = [
            [0.08, 0.70, spdim[0],spdim[1]],
            [0.40, 0.70, spdim[0],spdim[1]],
            [0.72, 0.70, spdim[0],spdim[1]],
            [0.08, 0.38, spdim[0],spdim[1]],
            [0.40, 0.38, spdim[0],spdim[1]],
            [0.72, 0.38, spdim[0],spdim[1]],
            [0.08, 0.06, spdim[0],spdim[1]],
            [0.40, 0.06, spdim[0],spdim[1]],
            [0.72, 0.06, spdim[0],spdim[1]],
            ]
 
subtpos = [0.01,1.06]
   
subt = [('(a)'), ('(b)'), ('(c)'), ('(d)'),
        ('(e)'), ('(f)'), ('(g)'), ('(h)'),
        ('(i)')]

#colors = ['r', 'g', 'b','k', 'c','m','y']

colors = list([(0,0,0), (230/255,159/255,0), (86/255,180/255,233/255),
               (0,158/255,115/255), (240/255,228/255,66/255),
               (0,114/255,178/255), (213/255,94/255,0),(204/255,121/255,167/255)])

labs = ['1', '2', '3', '4', '5', '6', '7']    

if regular:
    
    fig = plt.figure(num=500)
    
    fig.set_size_inches(6,5)
    plt.clf()
    
    ax1 = fig.add_subplot(3,3,1)
    #ax1.hist(diam_3_1[int_3_1], histtype='step', bins=100, density=True,  color = 'r', label='1')
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_1),mm.lognorm_pdf(diamfit_s_1_1),'-', color=colors[0],label=labs[0], lw=lw)
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_2),mm.lognorm_pdf(diamfit_s_1_2),'-', color=colors[1],label=labs[1], lw=lw)
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_3),mm.lognorm_pdf(diamfit_s_1_3),'-', color=colors[2],label=labs[2], lw=lw)
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_4),mm.lognorm_pdf(diamfit_s_1_4),'-', color=colors[3],label=labs[3], lw=lw)
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_5),mm.lognorm_pdf(diamfit_s_1_5),'-', color=colors[4],label=labs[4], lw=lw)
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_6),mm.lognorm_pdf(diamfit_s_1_6),'-', color=colors[5],label=labs[5], lw=lw)
    ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_7),mm.lognorm_pdf(diamfit_s_1_7),'-', color=colors[6],label=labs[6], lw=lw)
    ax1.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax1.xaxis, ax1.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax1.legend()
    #ax1.set_xlabel('Pore equivalent diameter (mm)')
    ax1.set_ylabel('Pore size density')
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymina,ymaxa)
    plt.text(0.10,0.90,'0\u20135 cm',transform=ax1.transAxes, fontsize=fontsize)
    
    ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)
    
    ax1.set_position(sppos[0]) 
    
    
    ax2 = fig.add_subplot(3,3,2)
    dummy = mm.vol_cumsum(volfact*vol_1_1, int_1_1)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[0],label=labs[0], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_2, int_1_2)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[1],label=labs[1], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_3, int_1_3)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[2],label=labs[2], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_4, int_1_4)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[3],label=labs[3],lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_5, int_1_5)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[4],label=labs[4], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_6, int_1_6)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[5],label=labs[5], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_7, int_1_7)
    ax2.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[6],label=labs[6], lw=lw)
    ax2.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax2.xaxis, ax2.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    ax2.legend(fontsize=fontsize, ncol=2, columnspacing=1,borderaxespad=1.0)
    ax2.set_xlim(xmin,xmax)
    ax2.set_ylim(-0.5,17)
    
    #ax2.set_xlabel('Pore equivalent diameter (mm)')
    ax2.set_ylabel(r'Cumulative volume (cm$^3$)')
    
    ax2.text(subtpos[0], subtpos[1], subt[1], transform=ax2.transAxes)
    ax2.set_position(sppos[1]) 
    
    
    ax3 = fig.add_subplot(3,3,3)
    dens = True
    hist_t_1 = ax3.hist(thr_1_1, histtype='step', bins=np.unique(thr_1_1), density=dens, cumulative=-1, color = colors[0], label='1')
    ax3.hist(thr_1_2, histtype='step', bins=np.unique(thr_1_2), density=dens, cumulative=-1, color = colors[1], label='2')
    ax3.hist(thr_1_3, histtype='step', bins=np.unique(thr_1_3), density=dens, cumulative=-1, color = colors[2], label='3')
    ax3.hist(thr_1_4, histtype='step', bins=np.unique(thr_1_4), density=dens, cumulative=-1, color = colors[3], label='4')
    ax3.hist(thr_1_5, histtype='step', bins=np.unique(thr_1_5), density=dens, cumulative=-1, color = colors[4], label='5')
    ax3.hist(thr_1_6, histtype='step', bins=np.unique(thr_1_6), density=dens, cumulative=-1, color = colors[5], label='6')
    ax3.hist(thr_1_7, histtype='step', bins=np.unique(thr_1_7), density=dens, cumulative=-1, color = colors[6], label='7')
    
    ax3.set_xscale('log')
    ax3.set_xticks([0.1, 0.2, 0.5, 1, 2])
    
    for axis in [ax3.xaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax3.legend()
    #ax3.set_xlabel('Throat diameter (mm)')
    
    ax3.set_ylabel('Cumulative frequency')
    ax3.set_xlim(xminc,xmaxc)
    ax3.set_ylim(yminc,ymaxc)
    
    ax2.text(subtpos[0], subtpos[1], subt[2], transform=ax3.transAxes)
    ax3.set_position(sppos[2]) 
    
    
    ax4 = fig.add_subplot(3,3,4)
    #ax4.hist(diam_3_1[int_3_1], histtype='step', bins=100, density=True,  color = 'r', label='1')
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_1),mm.lognorm_pdf(diamfit_s_2_1),'-', color=colors[0],label=labs[0], lw=lw)
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_2),mm.lognorm_pdf(diamfit_s_2_2),'-', color=colors[1],label=labs[1], lw=lw)
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_3),mm.lognorm_pdf(diamfit_s_2_3),'-', color=colors[2],label=labs[2], lw=lw)
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_4),mm.lognorm_pdf(diamfit_s_2_4),'-', color=colors[3],label=labs[3], lw=lw)
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_5),mm.lognorm_pdf(diamfit_s_2_5),'-', color=colors[4],label=labs[4], lw=lw)
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_6),mm.lognorm_pdf(diamfit_s_2_6),'-', color=colors[5],label=labs[5], lw=lw)
    ax4.semilogx(mm.lognorm_xaxis(diamfit_s_2_7),mm.lognorm_pdf(diamfit_s_2_7),'-', color=colors[6],label=labs[6], lw=lw)
    ax4.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax4.xaxis, ax4.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax4.legend()
    #ax4.set_xlabel('Pore equivalent diameter (mm)')
    ax4.set_ylabel('Pore size density')
    ax4.set_xlim(xmin,xmax)
    ax4.set_ylim(ymina,ymaxa)
    plt.text(0.10,0.90,'20\u201325 cm',transform=ax4.transAxes, fontsize=fontsize)
    
    ax4.text(subtpos[0], subtpos[1], subt[3], transform=ax4.transAxes)
    ax4.set_position(sppos[3]) 
    
    
    ax5 = fig.add_subplot(3,3,5)
    dummy = mm.vol_cumsum(volfact*vol_2_1, int_2_1)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[0],label=labs[0], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_2, int_2_2)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[1],label=labs[1], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_3, int_2_3)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[2],label=labs[2], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_4, int_2_4)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[3],label=labs[3], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_5, int_2_5)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[4],label=labs[4], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_6, int_2_6)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[5],label=labs[5], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_7, int_2_7)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[6],label=labs[6], lw=lw)
    ax5.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax5.xaxis, ax5.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax5.legend()
    ax5.set_xlim(xmin,xmax)
    ax5.set_ylim(yminb,ymaxb)
    
    ax5.text(subtpos[0], subtpos[1], subt[4], transform=ax5.transAxes)
    ax5.set_position(sppos[4]) 
    
    
    #ax5.set_xlabel('Pore equivalent diameter (mm)')
    ax5.set_ylabel('Cumulative volume (cm$^3$)')
    
    ax6 = fig.add_subplot(3,3,6)
    dens = True
    hist_t_1 = ax6.hist(thr_2_1, histtype='step', bins=np.unique(thr_2_1), density=dens, cumulative=-1, color = colors[0], label='1')
    ax6.hist(thr_2_2, histtype='step', bins=np.unique(thr_2_2), density=dens, cumulative=-1, color = colors[1], label='2')
    ax6.hist(thr_2_3, histtype='step', bins=np.unique(thr_2_3), density=dens, cumulative=-1, color = colors[2], label='3')
    ax6.hist(thr_2_4, histtype='step', bins=np.unique(thr_2_4), density=dens, cumulative=-1, color = colors[3], label='4')
    ax6.hist(thr_2_5, histtype='step', bins=np.unique(thr_2_5), density=dens, cumulative=-1, color = colors[4], label='5')
    ax6.hist(thr_2_6, histtype='step', bins=np.unique(thr_2_6), density=dens, cumulative=-1, color = colors[5], label='6')
    ax6.hist(thr_2_7, histtype='step', bins=np.unique(thr_2_7), density=dens, cumulative=-1, color = colors[6], label='7')
    
    ax6.set_xscale('log')
    ax6.set_xticks([0.1, 0.2, 0.5, 1, 2])
    
    for axis in [ax6.xaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax6.legend()
    #ax6.set_xlabel('Throat diameter (mm)')
    
    ax6.set_ylabel('Cumulative frequency')
    ax6.set_xlim(xminc,xmaxc)
    ax6.set_ylim(yminc,ymaxc)
    
    ax6.text(subtpos[0], subtpos[1], subt[5], transform=ax6.transAxes)
    ax6.set_position(sppos[5]) 
    
    
    ax7 = fig.add_subplot(3,3,7)
    #ax7.hist(diam_3_1[int_3_1], histtype='step', bins=100, density=True,  color = 'r', label='1')
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_1),mm.lognorm_pdf(diamfit_s_3_1),'-', color=colors[0],label=labs[0], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_2),mm.lognorm_pdf(diamfit_s_3_2),'-', color=colors[1],label=labs[1], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_3),mm.lognorm_pdf(diamfit_s_3_3),'-', color=colors[2],label=labs[2], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_4),mm.lognorm_pdf(diamfit_s_3_4),'-', color=colors[3],label=labs[3], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_5),mm.lognorm_pdf(diamfit_s_3_5),'-', color=colors[4],label=labs[4], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_6),mm.lognorm_pdf(diamfit_s_3_6),'-', color=colors[5],label=labs[5], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_7),mm.lognorm_pdf(diamfit_s_3_7),'-', color=colors[6],label=labs[6], lw=lw)
    ax7.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax7.xaxis, ax7.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax7.legend()
    ax7.set_xlabel('Pore equivalent diameter (mm)')
    ax7.set_ylabel('Pore size density')
    ax7.set_xlim(xmin,xmax)
    ax7.set_ylim(ymina,ymaxa)
    plt.text(0.10,0.90,'40\u201345 cm',transform=ax7.transAxes, fontsize=fontsize)
    
    ax7.text(subtpos[0], subtpos[1], subt[6], transform=ax7.transAxes)
    ax7.set_position(sppos[6]) 
    
    
    ax8 = fig.add_subplot(3,3,8)
    dummy = mm.vol_cumsum(volfact*vol_3_1, int_3_1)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[0],label=labs[0], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_2, int_3_2)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[1],label=labs[1], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_3, int_3_3)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[2],label=labs[2], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_4, int_3_4)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[3],label=labs[3], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_5, int_3_5)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[4],label=labs[4], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_6, int_3_6)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[5],label=labs[5], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_7, int_3_7)
    ax8.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[6],label=labs[6], lw=lw)
    ax8.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax8.xaxis, ax8.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    ax8.set_xlim(xmin,xmax)
    ax8.set_ylim(yminb,ymaxb)
    
    ax8.set_xlabel('Pore equivalent diameter (mm)')
    ax8.set_ylabel('Cumulative volume (cm$^3$)')
    
    ax8.text(subtpos[0], subtpos[1], subt[7], transform=ax8.transAxes)
    ax8.set_position(sppos[7]) 
    
    
    ax9 = fig.add_subplot(3,3,9)
    dens = True
    hist_t_1 = ax9.hist(thr_3_1, histtype='step', bins=np.unique(thr_3_1), density=dens, cumulative=-1, color = colors[0], label='1')
    ax9.hist(thr_3_2, histtype='step', bins=np.unique(thr_3_2), density=dens, cumulative=-1, color = colors[1], label='2')
    ax9.hist(thr_3_3, histtype='step', bins=np.unique(thr_3_3), density=dens, cumulative=-1, color = colors[2], label='3')
    ax9.hist(thr_3_4, histtype='step', bins=np.unique(thr_3_4), density=dens, cumulative=-1, color = colors[3], label='4')
    ax9.hist(thr_3_5, histtype='step', bins=np.unique(thr_3_5), density=dens, cumulative=-1, color = colors[4], label='5')
    ax9.hist(thr_3_6, histtype='step', bins=np.unique(thr_3_6), density=dens, cumulative=-1, color = colors[5], label='6')
    ax9.hist(thr_3_7, histtype='step', bins=np.unique(thr_3_7), density=dens, cumulative=-1, color = colors[6], label='7')
    
    ax9.set_xscale('log')
    ax9.set_xticks([0.1, 0.2, 0.5, 1, 2])
    
    for axis in [ax9.xaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax9.legend()
    ax9.set_xlabel('Throat diameter (mm)')
    
    ax9.set_ylabel('Cumulative frequency')
    ax9.set_xlim(xminc,xmaxc)
    ax9.set_ylim(yminc,ymaxc)
    
    ax9.text(subtpos[0], subtpos[1], subt[8], transform=ax9.transAxes)
    ax9.set_position(sppos[8]) 
    
    #ax8.legend(bbox_to_anchor=(0.88, 0.12), bbox_transform=plt.gcf().transFigure)
    
    ax2.yaxis.set_label_coords(-0.1, 0.5)
    ax5.yaxis.set_label_coords(-0.1, 0.5)
    ax8.yaxis.set_label_coords(-0.1, 0.5)
    
    ax3.yaxis.set_label_coords(-0.20, 0.5)
    ax6.yaxis.set_label_coords(-0.20, 0.5)
    ax9.yaxis.set_label_coords(-0.20, 0.5)
    
    
    ax7.xaxis.set_label_coords(0.5, -0.13)
    ax8.xaxis.set_label_coords(0.5, -0.13)
    ax9.xaxis.set_label_coords(0.5, -0.13)
    
    #plt.savefig('psd_v1.pdf')

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Distributions merged

if joint:
    
    colors_3 = ['k', 'b', 'r']
    
    spdim = [0.26, 0.41]
    spdim_2 = [0.32, 0.35]
    
    sppos = [
            [0.14, 0.63, spdim_2[0],spdim_2[1]], #
            [0.55, 0.63, spdim_2[0],spdim_2[1]], #
            [0.01, 0.01, spdim_2[0],spdim_2[1]],
            [0.06, 0.10, spdim[0],spdim[1]], #
            [0.39, 0.10, spdim[0],spdim[1]], #
            [0.72, 0.10, spdim[0],spdim[1]], #
            [0.08, 0.06, spdim[0],spdim[1]],
            [0.40, 0.06, spdim[0],spdim[1]],
            [0.72, 0.06, spdim[0],spdim[1]],
            ]
    
    subtpos = [0.05,0.88]
    
    fig = plt.figure(num=502)
    
    fig.set_size_inches(6,3.5)
    plt.clf()
    
    ax1 = fig.add_subplot(2,3,1)
    #ax1.hist(diam_3_1[int_3_1], histtype='step', bins=100, density=True,  color = 'r', label='1')
    bins1 =  np.linspace(np.log10(np.min(diam_1_all)), np.log10(np.max(diam_1_all)), num=40)
    testi1=ax1.hist(np.log10(diam_1_all), histtype='step', bins=bins1, density=True,  color = colors_3[0], label='0\u20135 cm')
    bins =  np.linspace(np.log10(np.min(diam_2_all)), np.log10(np.max(diam_2_all)), num=40)
    ax1.hist(np.log10(diam_2_all), histtype='step', bins=bins, density=True,  color = colors_3[1], label='20\u201325 cm')
    bins =  np.linspace(np.log10(np.min(diam_3_all)), np.log10(np.max(diam_3_all)), num=40)
    ax1.hist(np.log10(diam_3_all), histtype='step', bins=bins, density=True,  color = colors_3[2], label='40\u201345 cm')
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1),mm.lognorm_pdf(diamfit_s_1),'-', color=colors[0],label='0-5', lw=lw)
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_2),mm.lognorm_pdf(diamfit_s_2),'-', color=colors[1],label='20-25', lw=lw)
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_3),mm.lognorm_pdf(diamfit_s_3),'-', color=colors[2],label='40-45', lw=lw)
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_4),mm.lognorm_pdf(diamfit_s_1_4),'-', color=colors[3],label=labs[3], lw=lw)
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_5),mm.lognorm_pdf(diamfit_s_1_5),'-', color=colors[4],label=labs[4], lw=lw)
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_6),mm.lognorm_pdf(diamfit_s_1_6),'-', color=colors[5],label=labs[5], lw=lw)
    #ax1.semilogx(mm.lognorm_xaxis(diamfit_s_1_7),mm.lognorm_pdf(diamfit_s_1_7),'-', color=colors[6],label=labs[6], lw=lw)
    
    #ax1.set_xscale('log')
    #ax1.set_xticks([0.2, 0.5, 1, 2, 5])
    ax1.set_xticks([np.log10(0.2), np.log10(0.5), np.log10(1), np.log10(2), np.log10(5)])
    ax1.set_xticklabels(['0.2','0.5','1','2','5'])
    
    for axis in [ax1.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax1.legend()
    ax1.set_xlabel('Pore equivalent diameter (mm)')
    ax1.set_ylabel('Pore size density')
    ax1.set_xlim(np.log10(xmin),np.log10(xmax))
    ax1.set_ylim(ymina,ymaxa)
    #plt.text(0.10,0.90,'0\u20135 cm',transform=ax1.transAxes, fontsize=fontsize)
    
    ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)
    
    ax1.set_position(sppos[0]) 
    
    
    ax2 = fig.add_subplot(2,3,2)
    dens = True
    hist_t_1 = ax2.hist(thr_1, histtype='step', bins=np.unique(thr_1), density=dens, cumulative=-1, color = colors_3[0], label='0\u20135 cm')
    ax2.hist(thr_2, histtype='step', bins=np.unique(thr_2), density=dens, cumulative=-1, color = colors_3[1], label='20\u201325 cm')
    ax2.hist(thr_3, histtype='step', bins=np.unique(thr_3), density=dens, cumulative=-1, color = colors_3[2], label='40\u201345 cm')
    #ax2.hist(thr_1_4, histtype='step', bins=np.unique(thr_1_4), density=dens, cumulative=-1, color = colors[3], label='4')
    #ax2.hist(thr_1_5, histtype='step', bins=np.unique(thr_1_5), density=dens, cumulative=-1, color = colors[4], label='5')
    #ax2.hist(thr_1_6, histtype='step', bins=np.unique(thr_1_6), density=dens, cumulative=-1, color = colors[5], label='6')
    #ax2.hist(thr_1_7, histtype='step', bins=np.unique(thr_1_7), density=dens, cumulative=-1, color = colors[6], label='7')
    
    ax2.set_xscale('log')
    ax2.set_xticks([0.1, 0.2, 0.5, 1, 2])
    
    for axis in [ax2.xaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax2.legend()
    ax2.set_xlabel('Throat diameter (mm)')
    
    ax2.set_ylabel('Cumulative frequency')
    ax2.set_xlim(0.05,xmaxc)
    ax2.set_ylim(yminc,1.2)
    
    ax2.legend(fontsize=fontsize, borderaxespad=1.0)
    
    ax2.text(subtpos[0], subtpos[1], subt[1], transform=ax2.transAxes)
    ax2.set_position(sppos[1]) 
    
    ''' 
    ax3 = fig.add_subplot(2,3,3)
    ax3.semilogx(mm.lognorm_xaxis(diamfit_s_1),mm.lognorm_pdf(diamfit_s_1),'-', color=colors_3[0],label='0\u20135 cm', lw=lw)
    ax3.semilogx(mm.lognorm_xaxis(diamfit_s_2),mm.lognorm_pdf(diamfit_s_2),'-', color=colors_3[1],label='20\u201325 cm', lw=lw)
    ax3.semilogx(mm.lognorm_xaxis(diamfit_s_3),mm.lognorm_pdf(diamfit_s_3),'-', color=colors_3[2],label='40\u201345 cm', lw=lw)
    
    ax3.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax3.xaxis, ax3.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax3.set_xlabel('Pore equivalent diameter (mm)')
    ax3.set_ylabel('Pore size density')
    ax3.set_xlim(xmin,xmax)
    ax3.set_ylim(ymina,ymaxa)
    #plt.text(0.10,0.90,'20\u201325 cm',transform=ax3.transAxes, fontsize=fontsize)
    
    ax3.text(subtpos[0], subtpos[1], subt[2], transform=ax3.transAxes)
    ax3.set_position(sppos[2]) 
    '''
    
    ax4 = fig.add_subplot(2,3,4)
    dummy = mm.vol_cumsum(volfact*vol_1_1, int_1_1)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[0],label=labs[0], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_2, int_1_2)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[1],label=labs[1], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_3, int_1_3)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[2],label=labs[2], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_4, int_1_4)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[3],label=labs[3],lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_5, int_1_5)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[4],label=labs[4], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_6, int_1_6)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[5],label=labs[5], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_1_7, int_1_7)
    ax4.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[6],label=labs[6], lw=lw)
    ax4.set_xticks([0.2, 0.5, 1, 2, 5])
    
    ax4.text(0.16,0.88,'0\u20135 cm',transform=ax4.transAxes, fontsize=fontsize)
    
    for axis in [ax4.xaxis, ax4.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax4.legend(fontsize=fontsize, ncol=2, columnspacing=1,borderaxespad=1.0)
    ax4.set_xlim(xmin,xmax)
    ax4.set_ylim(-0.5,17)
    
    ax4.set_xlabel('Pore equivalent diameter (mm)')
    ax4.set_ylabel(r'Cumulative volume (cm$^3$)')
    
    ax4.text(subtpos[0], subtpos[1], subt[2], transform=ax4.transAxes)
    ax4.set_position(sppos[3]) 
    

    
    ax5 = fig.add_subplot(2,3,5)
    dummy = mm.vol_cumsum(volfact*vol_2_1, int_2_1)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[0],label=labs[0], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_2, int_2_2)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[1],label=labs[1], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_3, int_2_3)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[2],label=labs[2], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_4, int_2_4)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[3],label=labs[3], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_5, int_2_5)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[4],label=labs[4], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_6, int_2_6)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[5],label=labs[5], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_2_7, int_2_7)
    ax5.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[6],label=labs[6], lw=lw)
    ax5.set_xticks([0.2, 0.5, 1, 2, 5])
    
    ax5.text(0.16,0.88,'20\u201325 cm',transform=ax5.transAxes, fontsize=fontsize)
    
    for axis in [ax5.xaxis, ax5.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
        
    ax5.legend(fontsize=fontsize, ncol=2, columnspacing=0.9,#borderaxespad=1.0,
               handlelength=0.7, loc='upper left', bbox_to_anchor=(0.05, 0.71))
    ax5.text(0.08, 0.71, 'Sample number', transform=ax5.transAxes)
    
    ax5.set_xlim(xmin,xmax)
    ax5.set_ylim(yminb,ymaxb)
    
    ax5.text(subtpos[0], subtpos[1], subt[3], transform=ax5.transAxes)
    ax5.set_position(sppos[4]) 
    
    ax5.set_xlabel('Pore equivalent diameter (mm)')
    ax5.set_ylabel('Cumulative volume (cm$^3$)')
    
    
    '''
    ax6 = fig.add_subplot(3,3,6)
    dens = True
    hist_t_1 = ax6.hist(thr_2_1, histtype='step', bins=np.unique(thr_2_1), density=dens, cumulative=-1, color = colors[0], label='1')
    ax6.hist(thr_2_2, histtype='step', bins=np.unique(thr_2_2), density=dens, cumulative=-1, color = colors[1], label='2')
    ax6.hist(thr_2_3, histtype='step', bins=np.unique(thr_2_3), density=dens, cumulative=-1, color = colors[2], label='3')
    ax6.hist(thr_2_4, histtype='step', bins=np.unique(thr_2_4), density=dens, cumulative=-1, color = colors[3], label='4')
    ax6.hist(thr_2_5, histtype='step', bins=np.unique(thr_2_5), density=dens, cumulative=-1, color = colors[4], label='5')
    ax6.hist(thr_2_6, histtype='step', bins=np.unique(thr_2_6), density=dens, cumulative=-1, color = colors[5], label='6')
    ax6.hist(thr_2_7, histtype='step', bins=np.unique(thr_2_7), density=dens, cumulative=-1, color = colors[6], label='7')
    
    ax6.set_xscale('log')
    ax6.set_xticks([0.1, 0.2, 0.5, 1, 2])
    
    for axis in [ax6.xaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax6.legend()
    #ax6.set_xlabel('Throat diameter (mm)')
    
    ax6.set_ylabel('Cumulative frequency')
    ax6.set_xlim(xminc,xmaxc)
    ax6.set_ylim(yminc,ymaxc)
    
    ax6.text(subtpos[0], subtpos[1], subt[5], transform=ax6.transAxes)
    ax6.set_position(sppos[5]) 
    
    
    ax7 = fig.add_subplot(3,3,7)
    #ax7.hist(diam_3_1[int_3_1], histtype='step', bins=100, density=True,  color = 'r', label='1')
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_1),mm.lognorm_pdf(diamfit_s_3_1),'-', color=colors[0],label=labs[0], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_2),mm.lognorm_pdf(diamfit_s_3_2),'-', color=colors[1],label=labs[1], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_3),mm.lognorm_pdf(diamfit_s_3_3),'-', color=colors[2],label=labs[2], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_4),mm.lognorm_pdf(diamfit_s_3_4),'-', color=colors[3],label=labs[3], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_5),mm.lognorm_pdf(diamfit_s_3_5),'-', color=colors[4],label=labs[4], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_6),mm.lognorm_pdf(diamfit_s_3_6),'-', color=colors[5],label=labs[5], lw=lw)
    ax7.semilogx(mm.lognorm_xaxis(diamfit_s_3_7),mm.lognorm_pdf(diamfit_s_3_7),'-', color=colors[6],label=labs[6], lw=lw)
    ax7.set_xticks([0.2, 0.5, 1, 2, 5])
    
    for axis in [ax7.xaxis, ax7.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    #ax7.legend()
    ax7.set_xlabel('Pore equivalent diameter (mm)')
    ax7.set_ylabel('Pore size density')
    ax7.set_xlim(xmin,xmax)
    ax7.set_ylim(ymina,ymaxa)
    plt.text(0.10,0.90,'40\u201345 cm',transform=ax7.transAxes, fontsize=fontsize)
    
    ax7.text(subtpos[0], subtpos[1], subt[6], transform=ax7.transAxes)
    ax7.set_position(sppos[6]) 
    '''
    
    ax6 = fig.add_subplot(2,3,6)
    dummy = mm.vol_cumsum(volfact*vol_3_1, int_3_1)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[0],label=labs[0], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_2, int_3_2)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[1],label=labs[1], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_3, int_3_3)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[2],label=labs[2], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_4, int_3_4)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[3],label=labs[3], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_5, int_3_5)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[4],label=labs[4], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_6, int_3_6)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[5],label=labs[5], lw=lw)
    dummy = mm.vol_cumsum(volfact*vol_3_7, int_3_7)
    ax6.semilogx(dummy.x, 1e-3 * dummy.cumsum,'-', color=colors[6],label=labs[6], lw=lw)
    ax6.set_xticks([0.2, 0.5, 1, 2, 5])
    
    ax6.text(0.16,0.88,'40\u201345 cm',transform=ax6.transAxes, fontsize=fontsize)
    
    for axis in [ax6.xaxis, ax6.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    ax6.set_xlim(xmin,xmax)
    ax6.set_ylim(yminb,ymaxb)
    
    ax6.set_xlabel('Pore equivalent diameter (mm)')
    ax6.set_ylabel('Cumulative volume (cm$^3$)')
    
    ax6.text(subtpos[0], subtpos[1], subt[4], transform=ax6.transAxes)
    ax6.set_position(sppos[5]) 
    
    '''
    ax9 = fig.add_subplot(3,3,9)
    dens = True
    hist_t_1 = ax9.hist(thr_3_1, histtype='step', bins=np.unique(thr_3_1), density=dens, cumulative=-1, color = colors[0], label='1')
    ax9.hist(thr_3_2, histtype='step', bins=np.unique(thr_3_2), density=dens, cumulative=-1, color = colors[1], label='2')
    ax9.hist(thr_3_3, histtype='step', bins=np.unique(thr_3_3), density=dens, cumulative=-1, color = colors[2], label='3')
    ax9.hist(thr_3_4, histtype='step', bins=np.unique(thr_3_4), density=dens, cumulative=-1, color = colors[3], label='4')
    ax9.hist(thr_3_5, histtype='step', bins=np.unique(thr_3_5), density=dens, cumulative=-1, color = colors[4], label='5')
    ax9.hist(thr_3_6, histtype='step', bins=np.unique(thr_3_6), density=dens, cumulative=-1, color = colors[5], label='6')
    ax9.hist(thr_3_7, histtype='step', bins=np.unique(thr_3_7), density=dens, cumulative=-1, color = colors[6], label='7')
    
    ax9.set_xscale('log')
    ax9.set_xticks([0.1, 0.2, 0.5, 1, 2])
    
    for axis in [ax9.xaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
    
    #ax9.legend()
    ax9.set_xlabel('Throat diameter (mm)')
    
    ax9.set_ylabel('Cumulative frequency')
    ax9.set_xlim(xminc,xmaxc)
    ax9.set_ylim(yminc,ymaxc)
    
    ax9.text(subtpos[0], subtpos[1], subt[8], transform=ax9.transAxes)
    ax9.set_position(sppos[8]) 
    '''
    #ax8.legend(bbox_to_anchor=(0.88, 0.12), bbox_transform=plt.gcf().transFigure)
    
    #ax2.yaxis.set_label_coords(-0.1, 0.5)
    #ax5.yaxis.set_label_coords(-0.1, 0.5)
    #ax8.yaxis.set_label_coords(-0.1, 0.5)
    
    #ax3.yaxis.set_label_coords(-0.20, 0.5)
    #ax6.yaxis.set_label_coords(-0.20, 0.5)
    #ax9.yaxis.set_label_coords(-0.20, 0.5)
    
    ax1.xaxis.set_label_coords(0.5, -0.14)
    ax2.xaxis.set_label_coords(0.5, -0.14)
    
    ax4.xaxis.set_label_coords(0.5, -0.125)
    ax5.xaxis.set_label_coords(0.5, -0.125)
    ax6.xaxis.set_label_coords(0.5, -0.125)
    
    ax4.yaxis.set_label_coords(-0.11, 0.5)
    ax5.yaxis.set_label_coords(-0.09, 0.5)
    ax6.yaxis.set_label_coords(-0.09, 0.5)
    
    #plt.savefig('psd_v2.pdf')

if means:

    
    poremeans = []
    
    poremeans.append(np.mean(diam_1_1[int_1_1]))
    poremeans.append(np.mean(diam_1_2[int_1_2]))
    poremeans.append(np.mean(diam_1_3[int_1_3]))
    poremeans.append(np.mean(diam_1_4[int_1_4]))
    poremeans.append(np.mean(diam_1_5[int_1_5]))
    poremeans.append(np.mean(diam_1_6[int_1_6]))
    poremeans.append(np.mean(diam_1_7[int_1_7]))
    poremeans.append(np.mean(diam_2_1[int_2_1]))
    poremeans.append(np.mean(diam_2_2[int_2_2]))
    poremeans.append(np.mean(diam_2_3[int_2_3]))
    poremeans.append(np.mean(diam_2_4[int_2_4]))
    poremeans.append(np.mean(diam_2_5[int_2_5]))
    poremeans.append(np.mean(diam_2_6[int_2_6]))
    poremeans.append(np.mean(diam_2_7[int_2_7]))
    poremeans.append(np.mean(diam_3_1[int_3_1]))
    poremeans.append(np.mean(diam_3_2[int_3_2]))
    poremeans.append(np.mean(diam_3_3[int_3_3]))
    poremeans.append(np.mean(diam_3_4[int_3_4]))
    poremeans.append(np.mean(diam_3_5[int_3_5]))
    poremeans.append(np.mean(diam_3_6[int_3_6]))
    poremeans.append(np.mean(diam_3_7[int_3_7]))
    
    poremeans = np.asarray(poremeans)
    
    print('Pore means:')
    print(poremeans[0:7])
    print(poremeans[7:14])
    print(poremeans[14:21])
    print('')
    
    poremedians = []
    
    poremedians.append(np.median(diam_1_1[int_1_1]))
    poremedians.append(np.median(diam_1_2[int_1_2]))
    poremedians.append(np.median(diam_1_3[int_1_3]))
    poremedians.append(np.median(diam_1_4[int_1_4]))
    poremedians.append(np.median(diam_1_5[int_1_5]))
    poremedians.append(np.median(diam_1_6[int_1_6]))
    poremedians.append(np.median(diam_1_7[int_1_7]))
    poremedians.append(np.median(diam_2_1[int_2_1]))
    poremedians.append(np.median(diam_2_2[int_2_2]))
    poremedians.append(np.median(diam_2_3[int_2_3]))
    poremedians.append(np.median(diam_2_4[int_2_4]))
    poremedians.append(np.median(diam_2_5[int_2_5]))
    poremedians.append(np.median(diam_2_6[int_2_6]))
    poremedians.append(np.median(diam_2_7[int_2_7]))
    poremedians.append(np.median(diam_3_1[int_3_1]))
    poremedians.append(np.median(diam_3_2[int_3_2]))
    poremedians.append(np.median(diam_3_3[int_3_3]))
    poremedians.append(np.median(diam_3_4[int_3_4]))
    poremedians.append(np.median(diam_3_5[int_3_5]))
    poremedians.append(np.median(diam_3_6[int_3_6]))
    poremedians.append(np.median(diam_3_7[int_3_7]))
    
    poremedians = np.asarray(poremedians)
    
    print('Pore medians:')
    print(poremedians[0:7])
    print(poremedians[7:14])
    print(poremedians[14:21])
    print('')
    
    
    throatmeans = []
    
    throatmeans.append(np.mean(thr_1_1))
    throatmeans.append(np.mean(thr_1_2))
    throatmeans.append(np.mean(thr_1_3))
    throatmeans.append(np.mean(thr_1_4))
    throatmeans.append(np.mean(thr_1_5))
    throatmeans.append(np.mean(thr_1_6))
    throatmeans.append(np.mean(thr_1_7))
    throatmeans.append(np.mean(thr_2_1))
    throatmeans.append(np.mean(thr_2_2))
    throatmeans.append(np.mean(thr_2_3))
    throatmeans.append(np.mean(thr_2_4))
    throatmeans.append(np.mean(thr_2_5))
    throatmeans.append(np.mean(thr_2_6))
    throatmeans.append(np.mean(thr_2_7))
    throatmeans.append(np.mean(thr_3_1))
    throatmeans.append(np.mean(thr_3_2))
    throatmeans.append(np.mean(thr_3_3))
    throatmeans.append(np.mean(thr_3_4))
    throatmeans.append(np.mean(thr_3_5))
    throatmeans.append(np.mean(thr_3_6))
    throatmeans.append(np.mean(thr_3_7))
    
    throatmeans = np.asarray(throatmeans)
    
    print('Throats: ')
    print(throatmeans[0:7])
    print(throatmeans[7:14])
    print(throatmeans[14:21])
    print('')
    
    throatmeans_int = []
    
    throatmeans_int.append(np.mean(thr_int_1_1))
    throatmeans_int.append(np.mean(thr_int_1_2))
    throatmeans_int.append(np.mean(thr_int_1_3))
    throatmeans_int.append(np.mean(thr_int_1_4))
    throatmeans_int.append(np.mean(thr_int_1_5))
    throatmeans_int.append(np.mean(thr_int_1_6))
    throatmeans_int.append(np.mean(thr_int_1_7))
    throatmeans_int.append(np.mean(thr_int_2_1))
    throatmeans_int.append(np.mean(thr_int_2_2))
    throatmeans_int.append(np.mean(thr_int_2_3))
    throatmeans_int.append(np.mean(thr_int_2_4))
    throatmeans_int.append(np.mean(thr_int_2_5))
    throatmeans_int.append(np.mean(thr_int_2_6))
    throatmeans_int.append(np.mean(thr_int_2_7))
    throatmeans_int.append(np.mean(thr_int_3_1))
    throatmeans_int.append(np.mean(thr_int_3_2))
    throatmeans_int.append(np.mean(thr_int_3_3))
    throatmeans_int.append(np.mean(thr_int_3_4))
    throatmeans_int.append(np.mean(thr_int_3_5))
    throatmeans_int.append(np.mean(thr_int_3_6))
    throatmeans_int.append(np.mean(thr_int_3_7))
    
    throatmeans_int = np.asarray(throatmeans_int)
    
    print('Internal Throats: ')
    print(throatmeans_int[0:7])
    print(throatmeans_int[7:14])
    print(throatmeans_int[14:21])
    print('')

if ppr2:
    
    pp_r2 = []
    
    #a = stats.probplot(diam_2_6[int_2_6], dist = stats.lognorm, sparams=diamfit_s_2_6)
    #b = a[1][2]**2
    pp_r2.append(mm.prob_plot_r2(diam_1_1, int_1_1, diamfit_s_1_1))
    pp_r2.append(mm.prob_plot_r2(diam_1_2, int_1_2, diamfit_s_1_2))
    pp_r2.append(mm.prob_plot_r2(diam_1_3, int_1_3, diamfit_s_1_3))
    pp_r2.append(mm.prob_plot_r2(diam_1_4, int_1_4, diamfit_s_1_4))
    pp_r2.append(mm.prob_plot_r2(diam_1_5, int_1_5, diamfit_s_1_5))
    pp_r2.append(mm.prob_plot_r2(diam_1_6, int_1_6, diamfit_s_1_6))
    pp_r2.append(mm.prob_plot_r2(diam_1_7, int_1_7, diamfit_s_1_7))
    pp_r2.append(mm.prob_plot_r2(diam_2_1, int_2_1, diamfit_s_2_1))
    pp_r2.append(mm.prob_plot_r2(diam_2_2, int_2_2, diamfit_s_2_2))
    pp_r2.append(mm.prob_plot_r2(diam_2_3, int_2_3, diamfit_s_2_3))
    pp_r2.append(mm.prob_plot_r2(diam_2_4, int_2_4, diamfit_s_2_4))
    pp_r2.append(mm.prob_plot_r2(diam_2_5, int_2_5, diamfit_s_2_5))
    pp_r2.append(mm.prob_plot_r2(diam_2_6, int_2_6, diamfit_s_2_6))
    pp_r2.append(mm.prob_plot_r2(diam_2_7, int_2_7, diamfit_s_2_7))
    pp_r2.append(mm.prob_plot_r2(diam_3_1, int_3_1, diamfit_s_3_1))
    pp_r2.append(mm.prob_plot_r2(diam_3_2, int_3_2, diamfit_s_3_2))
    pp_r2.append(mm.prob_plot_r2(diam_3_3, int_3_3, diamfit_s_3_3))
    pp_r2.append(mm.prob_plot_r2(diam_3_4, int_3_4, diamfit_s_3_4))
    pp_r2.append(mm.prob_plot_r2(diam_3_5, int_3_5, diamfit_s_3_5))
    pp_r2.append(mm.prob_plot_r2(diam_3_6, int_3_6, diamfit_s_3_6))
    pp_r2.append(mm.prob_plot_r2(diam_3_7, int_3_7, diamfit_s_3_7))
    
    pp_r2 = np.asarray(pp_r2)
    