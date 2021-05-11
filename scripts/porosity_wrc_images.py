# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:14:49 2020

@author: pkiuru

The scripts creates images
(a) image porosity vs. measured air-filled porosity @ 10 kPa
(b) simulated vs. measured air intrusion curves (air-filled porosity)

Based on poredist.py (250321 ->)

"""

import numpy as np
import matplotlib.pyplot as plt
import my_models as mm
import scipy.stats as stats
from matplotlib.ticker import FuncFormatter, FormatStrFormatter, MaxNLocator
from scipy.interpolate import interp1d

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

rcParams['legend.handlelength'] = 1.25


#Figure 3 (measured vs. image-derived air-filled porosities)
cylpor = 0

#Figure 4 (water retention curves w/ air-filled porosity)
perc_switch = 1
load_perc = 1
yformat = 'negative'


# -------------- POROSITIES ----------

if cylpor:

    #'''
    
    p_me = np.loadtxt('porosity_meas_tot_volume.txt')
    p_im = np.loadtxt('porosity_image.txt')
    
    #https://stackoverflow.com/questions/6148207/linear-regression-with-matplotlib-numpy
    
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
    slope, intercept, r_value, p_value, std_err = stats.linregress(p_me, p_im)
    
    fig = plt.figure(num=88)
    fig.set_size_inches(3.2,3.2)
    plt.clf()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(p_me[0:7], p_im[0:7], marker='o', s=14, c='k', label='0\u20135 cm',zorder=3)
    ax.scatter(p_me[7:14], p_im[7:14], marker='s', s=14, c='b', label='20\u201325 cm',zorder=4)
    ax.scatter(p_me[14:21], p_im[14:21], marker='^', s=16, c='r', label='40\u201345 cm',zorder=5)
    ax.plot(np.linspace(np.min(p_me),np.max(p_me),100), np.linspace(np.min(p_me),
                        np.max(p_me),100), 'k:', lw=1, label='1:1 line',zorder=2)
    ax.plot(p_me, intercept + slope * p_me, 'k', lw=1, label = 'Linear fit',zorder=1)
    
    formatter = FormatStrFormatter('%1.1f')
    locator = MaxNLocator(nbins=4, steps=[2])
    ax.yaxis.set_major_locator(locator)
    ax.yaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    
    ax.set_xlabel('Measured air-filled porosity')
    ax.set_ylabel('Image porosity')
    
    
    ax.text(0.53 ,0.27, '$y = $' + np.format_float_positional(slope, 3) + '$x + $'
            + np.format_float_positional(intercept, 3), transform=ax.transAxes)
    ax.text(0.53, 0.21, '$R^2 = $' + np.format_float_positional(r_value**2, 2),
            transform=ax.transAxes)
    
    lims = [0.09,0.81]
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    
    handles, labels = ax.get_legend_handles_labels()
    order = [2,3,4,1,0]
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
               fontsize=fontsize, loc='upper left', borderaxespad=1.5)
    
    ax.yaxis.set_label_coords(-0.085, 0.5)
    ax.xaxis.set_label_coords(0.5, -0.065)
    
    ax.set_position([0.12, 0.12, 0.85, 0.85])
    
    #plt.savefig('im_meas_poros.pdf')
    
    #'''
# ooooooooooooooooooooooooooooooo PERCOLATION CURVES ooooooooooooooooooooooooooooooooooo

#'''

if perc_switch:
    
    if load_perc:
        
        caf = 0.785321
        cacorr = np.cos(np.radians(180)) / np.cos(np.radians(120))
        
        a = '0_5'
        press_1_1, afp_1_1 = mm.load_percolation_curve(a+'_1')
        press_1_2, afp_1_2 = mm.load_percolation_curve(a+'_2')
        press_1_3, afp_1_3 = mm.load_percolation_curve(a+'_3')
        press_1_4, afp_1_4 = mm.load_percolation_curve(a+'_4')
        press_1_5, afp_1_5 = mm.load_percolation_curve(a+'_5')
        press_1_6, afp_1_6 = mm.load_percolation_curve(a+'_6')
        press_1_7, afp_1_7 = mm.load_percolation_curve(a+'_7')
        
        a = '20_25'
        press_2_1, afp_2_1 = mm.load_percolation_curve(a+'_1')
        press_2_2, afp_2_2 = mm.load_percolation_curve(a+'_2')
        press_2_3, afp_2_3 = mm.load_percolation_curve(a+'_3')
        press_2_4, afp_2_4 = mm.load_percolation_curve(a+'_4')
        press_2_5, afp_2_5 = mm.load_percolation_curve(a+'_5')
        press_2_6, afp_2_6 = mm.load_percolation_curve(a+'_6')
        press_2_7, afp_2_7 = mm.load_percolation_curve(a+'_7')
        
        a = '40_45'
        press_3_1, afp_3_1 = mm.load_percolation_curve(a+'_1')
        press_3_2, afp_3_2 = mm.load_percolation_curve(a+'_2')
        press_3_3, afp_3_3 = mm.load_percolation_curve(a+'_3_2')
        press_3_4, afp_3_4 = mm.load_percolation_curve(a+'_4')
        press_3_5, afp_3_5 = mm.load_percolation_curve(a+'_5_2')
        press_3_6, afp_3_6 = mm.load_percolation_curve(a+'_6')
        press_3_7, afp_3_7 = mm.load_percolation_curve(a+'_7')
        
        
        #Contact angle correction: 120 degrees -> 180 degrees
        
        p_max = 4400
        
        press_1_1 = cacorr * np.concatenate([press_1_1, np.array([4400/cacorr])])
        press_1_2 = cacorr * np.concatenate([press_1_2, np.array([4400/cacorr])])
        press_1_3 = cacorr * np.concatenate([press_1_3, np.array([4400/cacorr])])
        press_1_4 = cacorr * np.concatenate([press_1_4, np.array([4400/cacorr])])
        press_1_5 = cacorr * np.concatenate([press_1_5, np.array([4400/cacorr])])
        press_1_6 = cacorr * np.concatenate([press_1_6, np.array([4400/cacorr])])
        press_1_7 = cacorr * np.concatenate([press_1_7, np.array([4400/cacorr])])
        
        press_2_1 = cacorr * np.concatenate([press_2_1, np.array([4400/cacorr])])
        press_2_2 = cacorr * np.concatenate([press_2_2, np.array([4400/cacorr])])
        press_2_3 = cacorr * np.concatenate([press_2_3, np.array([4400/cacorr])])
        press_2_4 = cacorr * np.concatenate([press_2_4, np.array([4400/cacorr])])
        press_2_5 = cacorr * np.concatenate([press_2_5, np.array([4400/cacorr])])
        press_2_6 = cacorr * np.concatenate([press_2_6, np.array([4400/cacorr])])
        press_2_7 = cacorr * np.concatenate([press_2_7, np.array([4400/cacorr])])
        
        press_3_1 = cacorr * np.concatenate([press_3_1, np.array([4400/cacorr])])
        press_3_2 = cacorr * np.concatenate([press_3_2, np.array([4400/cacorr])])
        press_3_3 = cacorr * np.concatenate([press_3_3, np.array([4400/cacorr])])
        press_3_4 = cacorr * np.concatenate([press_3_4, np.array([4400/cacorr])])
        press_3_5 = cacorr * np.concatenate([press_3_5, np.array([4400/cacorr])])
        press_3_6 = cacorr * np.concatenate([press_3_6, np.array([4400/cacorr])])
        press_3_7 = cacorr * np.concatenate([press_3_7, np.array([4400/cacorr])])
        
        
        #Area correction: square -> cylinder
        
        afp_1_1 = np.concatenate([afp_1_1, np.array([afp_1_1[-1]])]) / caf
        afp_1_2 = np.concatenate([afp_1_2, np.array([afp_1_2[-1]])]) / caf
        afp_1_3 = np.concatenate([afp_1_3, np.array([afp_1_3[-1]])]) / caf
        afp_1_4 = np.concatenate([afp_1_4, np.array([afp_1_4[-1]])]) / caf
        afp_1_5 = np.concatenate([afp_1_5, np.array([afp_1_5[-1]])]) / caf
        afp_1_6 = np.concatenate([afp_1_6, np.array([afp_1_6[-1]])]) / caf
        afp_1_7 = np.concatenate([afp_1_7, np.array([afp_1_7[-1]])]) / caf
        
        afp_2_1 = np.concatenate([afp_2_1, np.array([afp_2_1[-1]])]) / caf
        afp_2_2 = np.concatenate([afp_2_2, np.array([afp_2_2[-1]])]) / caf
        afp_2_3 = np.concatenate([afp_2_3, np.array([afp_2_3[-1]])]) / caf
        afp_2_4 = np.concatenate([afp_2_4, np.array([afp_2_4[-1]])]) / caf
        afp_2_5 = np.concatenate([afp_2_5, np.array([afp_2_5[-1]])]) / caf
        afp_2_6 = np.concatenate([afp_2_6, np.array([afp_2_6[-1]])]) / caf
        afp_2_7 = np.concatenate([afp_2_7, np.array([afp_2_7[-1]])]) / caf
        
        afp_3_1 = np.concatenate([afp_3_1, np.array([afp_3_1[-1]])]) / caf
        afp_3_2 = np.concatenate([afp_3_2, np.array([afp_3_2[-1]])]) / caf
        afp_3_3 = np.concatenate([afp_3_3, np.array([afp_3_3[-1]])]) / caf
        afp_3_4 = np.concatenate([afp_3_4, np.array([afp_3_4[-1]])]) / caf
        afp_3_5 = np.concatenate([afp_3_5, np.array([afp_3_5[-1]])]) / caf
        afp_3_6 = np.concatenate([afp_3_6, np.array([afp_3_6[-1]])]) / caf
        afp_3_7 = np.concatenate([afp_3_7, np.array([afp_3_7[-1]])]) / caf

    

    
    colors = list([(0,0,0), (230/255,159/255,0), (86/255,180/255,233/255),
                   (0,158/255,115/255), (240/255,228/255,66/255),
                   (0,114/255,178/255), (213/255,94/255,0),(204/255,121/255,167/255)]) 
    
    phi = np.array([0.1, 1, 3, 6, 10])
    
    

    # ==================== WRC w.r.t. air-filled porosity ========================
    if True:
        #totpor = np.loadtxt('total_porosity_meas.txt')
        
        thetas = np.loadtxt('wrc_meas_tot_volume_with_zero.txt')
        
        totpor = thetas[:,0]
        
        lw = 1
        ylims = [0.06, 15]
        
        depthpos = [0.08, 0.92]
        
        spdim = [0.27, 0.79]
    
        sppos = [
                [0.07, 0.13, spdim[0],spdim[1]],
                [0.39, 0.13, spdim[0],spdim[1]],
                [0.71, 0.13, spdim[0],spdim[1]]
                ]
        
        subtpos = [0.0,1.04]
        subt = [('(a)'), ('(b)'), ('(c)')]
    
        
        fig = plt.figure(num=309)
        fig.set_size_inches(6,2.25)
        plt.clf()
        
        ax1 = fig.add_subplot(1,3,1)
        
        k = 1
        j = 0
        for i in range(7):
        
            theta = totpor[i+j] - thetas[i+j,:]
            
            f2 = interp1d(np.flip(theta), np.flip(np.log10(phi)), kind='linear')
            theta_fit_spline = np.linspace(theta[-1], theta[0],10000)
            
            ax1.plot(theta, phi, 's', c=colors[i], label = str(k), markersize=3)
        
            ax1.plot(theta_fit_spline, 10**(f2(theta_fit_spline)), ':', c=colors[i], lw=lw)
            print(i)
            if i==0:
                ax1.plot(afp_1_1, 0.001*press_1_1, '-', c=colors[i], lw=lw)
            if i==1:
                ax1.plot(afp_1_2, 0.001*press_1_2, '-', c=colors[i], lw=lw)
            if i==2:
                ax1.plot(afp_1_3, 0.001*press_1_3, '-', c=colors[i], lw=lw)
            if i==3:
                ax1.plot(afp_1_4, 0.001*press_1_4, '-', c=colors[i], lw=lw)
            if i==4:
                ax1.plot(afp_1_5, 0.001*press_1_5, '-', c=colors[i], lw=lw)
            if i==5:
                ax1.plot(afp_1_6, 0.001*press_1_6, '-', c=colors[i], lw=lw)
            if i==6:
                ax1.plot(afp_1_7, 0.001*press_1_7, '-', c=colors[i], lw=lw)
            
            k += 1
            
        ax1.set_yscale('log')
        #ax1.legend(handletextpad=0.2, ncol=2, columnspacing=0.5, borderaxespad=1)
        ax1.set_xlabel('Air-filled porosity')# (m$^3$/m$^3$)')
        ax1.set_ylabel('Matric potential (kPa)')
        ax1.set_ylim(ylims)
        
        if yformat == 'positive':
        
            for axis in [ax1.yaxis]:
                formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
                axis.set_major_formatter(formatter)
        else:
            ax1.set_yticks([0.1,1, 10])
            ax1.set_yticklabels(['-0.1','-1', '-10'])
        
        ax1.set_xlim(-0.08,0.80)
        ax1.text(depthpos[0], depthpos[1], '0\u20135 cm', transform=ax1.transAxes)
        #plt.text(-0.01,9,'0\u20135 cm')
        
        ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)
        
        ax1.set_position(sppos[0]) 
        
        ax2 = fig.add_subplot(1,3,2)
        
        k = 1
        j = 7
        for i in range(7):
        
            theta = totpor[i+j] - thetas[i+j,:]
            
            f2 = interp1d(np.flip(theta), np.flip(np.log10(phi)), kind='linear')
            theta_fit_spline = np.linspace(theta[-1], theta[0],10000)
            
            ax2.plot(theta, phi, 's', c=colors[i], label = str(k), markersize=3)
        
            ax2.plot(theta_fit_spline, 10**(f2(theta_fit_spline)), ':', c=colors[i], lw=lw)
            print(i)
            if i==0:
                ax2.plot(afp_2_1, 0.001*press_2_1, '-', c=colors[i], lw=lw)
            if i==1:
                ax2.plot(afp_2_2, 0.001*press_2_2, '-', c=colors[i], lw=lw)
            if i==2:
                ax2.plot(afp_2_3, 0.001*press_2_3, '-', c=colors[i], lw=lw)
            if i==3:
                ax2.plot(afp_2_4, 0.001*press_2_4, '-', c=colors[i], lw=lw)
            if i==4:
                ax2.plot(afp_2_5, 0.001*press_2_5, '-', c=colors[i], lw=lw)
            if i==5:
                ax2.plot(afp_2_6, 0.001*press_2_6, '-', c=colors[i], lw=lw)
            if i==6:
                ax2.plot(afp_2_7, 0.001*press_2_7, '-', c=colors[i], lw=lw)
            
            k += 1
            
        ax2.set_yscale('log')
        #ax2.legend(handletextpad=0.2, ncol=2, columnspacing=0.5, borderaxespad=1)
        ax2.set_xlabel('Air-filled porosity')# (m$^3$/m$^3$)')
        #ax2.set_ylabel('Pressure [kPa]')
        ax2.set_ylim(ylims)
        
        if yformat == 'positive':
        
            for axis in [ax2.yaxis]:
                formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
                axis.set_major_formatter(formatter)
        else:
            ax2.set_yticks([0.1,1, 10])
            ax2.set_yticklabels(['-0.1','-1', '-10'])
        
        ax2.set_xlim(-0.03,0.37)
        
        ax2.text(depthpos[0], depthpos[1], '20\u201325 cm', transform=ax2.transAxes)
        #plt.text(-0.03,9,'20\u201325 cm')
        
        ax2.text(subtpos[0], subtpos[1], subt[1], transform=ax2.transAxes)
        
        ax2.set_position(sppos[1]) 
        
        
        ax3 = fig.add_subplot(1,3,3)
        
        k = 1
        j = 14
        for i in range(7):
        
            theta = totpor[i+j] - thetas[i+j,:]
            
            f2 = interp1d(np.flip(theta), np.flip(np.log10(phi)), kind='linear')
            theta_fit_spline = np.linspace(theta[-1], theta[0],10000)
            
            ax3.plot(theta, phi, 's', c=colors[i], label = str(k), markersize=3)
        
            ax3.plot(theta_fit_spline, 10**(f2(theta_fit_spline)), ':', c=colors[i], lw=lw)
            print(i)
            if i==0:
                ax3.plot(afp_3_1, 0.001*press_3_1, '-', c=colors[i], lw=lw)
            if i==1:
                ax3.plot(afp_3_2, 0.001*press_3_2, '-', c=colors[i], lw=lw)
            if i==2:
                ax3.plot(afp_3_3, 0.001*press_3_3, '-', c=colors[i], lw=lw)
            if i==3:
                ax3.plot(afp_3_4, 0.001*press_3_4, '-', c=colors[i], lw=lw)
            if i==4:
                ax3.plot(afp_3_5, 0.001*press_3_5, '-', c=colors[i], lw=lw)
            if i==5:
                ax3.plot(afp_3_6, 0.001*press_3_6, '-', c=colors[i], lw=lw)
            if i==6:
                ax3.plot(afp_3_7, 0.001*press_3_7, '-', c=colors[i], lw=lw)
            
            k += 1
            
        ax3.set_yscale('log')
        ax3.legend(handletextpad=0.2, ncol=2, columnspacing=0.5, borderaxespad=1)
        ax3.set_xlabel('Air-filled porosity')# (m$^3$/m$^3$)')
        #ax3.set_ylabel('Pressure [kPa]')
        ax3.set_ylim(ylims)
        
        if yformat == 'positive':
        
            for axis in [ax3.yaxis]:
                formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
                axis.set_major_formatter(formatter)
        else:
            ax3.set_yticks([0.1,1, 10])
            ax3.set_yticklabels(['-0.1','-1', '-10'])
        
        ax3.set_xlim(-0.03,0.31)
        ax3.text(depthpos[0], depthpos[1], '40\u201345 cm', transform=ax3.transAxes)
        #plt.text(-0.03,9,'40\u201345 cm')
        
        ax3.text(subtpos[0], subtpos[1], subt[2], transform=ax3.transAxes)
        
        ax3.set_position(sppos[2]) 
        
        ax1.yaxis.set_label_coords(-0.14, 0.5)
    
        ax1.xaxis.set_label_coords(0.5, -0.09)
        ax2.xaxis.set_label_coords(0.5, -0.09)
        ax3.xaxis.set_label_coords(0.5, -0.09)
        
        #plt.savefig("wrcs.pdf")
