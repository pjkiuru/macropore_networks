# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:18:16 2020

@author: pkiuru

The script creates the hysteresis image.

Similar to script hystfig.py (<- 250321) but does not need poredist.py

"""

import numpy as np
import matplotlib.pyplot as plt
import more_models as mm
from matplotlib.ticker import FuncFormatter

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

load_perc = True

datafile = '../IO/wrc_meas_tot_volume_with_zero.txt'

# Retention & imbibition

if load_perc:
        
        caf = 0.785321
        cacorr = np.cos(np.radians(180)) / np.cos(np.radians(120))
        
        a = '0_5'
        press_1_6, afp_1_6 = mm.load_percolation_curve(a+'_6')
        
        a = '20_25'
        press_2_6, afp_2_6 = mm.load_percolation_curve(a+'_6')
        
        a = '40_45'
        press_3_5, afp_3_5 = mm.load_percolation_curve(a+'_5_2')
        
        a = '0_5'
        #Pressure sign is correct in press_w & air_filled porosity in afp_w! 
        press_1_6_w, afp_1_6_w = mm.load_water_percolation_curve(a+'_6')

        a = '20_25'
        press_2_6_w, afp_2_6_w = mm.load_water_percolation_curve(a+'_6')
        
        a = '40_45'
        press_3_5_w, afp_3_5_w = mm.load_water_percolation_curve(a+'_5')
        
        
        #Contact angle correction: 120 degrees -> 180 degrees
        
        p_max = 4400
        
        press_1_6 = cacorr * np.concatenate([press_1_6, np.array([4400/cacorr])])
        
        press_2_6 = cacorr * np.concatenate([press_2_6, np.array([4400/cacorr])])     

        press_3_5 = cacorr * np.concatenate([press_3_5, np.array([4400/cacorr])])
        
        press_1_6_w = np.flip(cacorr * np.concatenate([1.001*np.array([press_1_6_w[0]]), press_1_6_w]))
        
        press_2_6_w = np.flip(cacorr * np.concatenate([1.001*np.array([press_2_6_w[0]]), press_2_6_w]))
        
        press_3_5_w = np.flip(cacorr * np.concatenate([1.001*np.array([press_3_5_w[0]]), press_3_5_w]))
        
        #Area correction: square -> cylinder
        

        afp_1_6 = np.concatenate([afp_1_6, np.array([afp_1_6[-1]])]) / caf
       
        afp_2_6 = np.concatenate([afp_2_6, np.array([afp_2_6[-1]])]) / caf

        afp_3_5 = np.concatenate([afp_3_5, np.array([afp_3_5[-1]])]) / caf

        afp_1_6_w = np.flip(np.concatenate([np.array([afp_1_6_w[0]]), afp_1_6_w]) / caf)
        
        afp_2_6_w = np.flip(np.concatenate([np.array([afp_2_6_w[0]]), afp_2_6_w]) / caf)
        
        afp_3_5_w = np.flip(np.concatenate([np.array([afp_3_5_w[0]]), afp_3_5_w]) / caf)


spdim = [0.26, 0.67]

sppos = [
            [0.07, 0.20, spdim[0],spdim[1]],
            [0.38, 0.20, spdim[0],spdim[1]],
            [0.69, 0.20, spdim[0],spdim[1]]
            ]
 
subtpos = [0.01,1.06]
   
subt = [('(a)'), ('(b)'), ('(c)')]

phi = np.array([0.1, 1, 3, 6, 10])

colors = list([(0,0,0), (230/255,159/255,0), (86/255,180/255,233/255),
                   (0,158/255,115/255), (240/255,228/255,66/255),
                   (0,114/255,178/255), (213/255,94/255,0),(204/255,121/255,167/255)])

thetas = np.loadtxt(datafile)      
totpor = thetas[:,0]
 
lw = 1
y_min = 0.04
yformat = 'negative'

fig = plt.figure(num=600)
fig.set_size_inches(6,1.5)
plt.clf()
ax1 = fig.add_subplot(1,3,1)

k = 1
a = '0_5'
i = 5

afp_drain = afp_1_6
press_drain = press_1_6
afp_imb = afp_1_6_w
press_imb = np.copy(press_1_6_w)
press_imb[press_imb>p_max] = p_max

    
theta = totpor[i+0] - thetas[i+0,:]

#f2 = interp1d(np.flip(theta), np.flip(np.log10(phi)), kind='linear')
#theta_fit_spline = np.linspace(theta[-1], theta[0],10000)

ax1.plot(theta, phi, 's', c=colors[0], label = 'sample ' + str(i+1), markersize=3)

#ax1.plot(theta_fit_spline, 10**(f2(theta_fit_spline)), ':', c=colors[0])

ax1.plot(afp_drain, 0.001*press_drain, '-', c=colors[0], lw=lw,label='Drainage')
ax1.plot(afp_imb, 0.001*press_imb, '--', c=colors[0], lw=lw,label='Imbibition')
ax1.fill_between(afp_imb, 0.001*press_imb, 0.001*np.interp(afp_imb, afp_drain, press_drain), alpha=0.2)


ax1.set_yscale('log')
#ax1.legend(handletextpad=0.2)
ax1.set_xlabel('Air-filled porosity')# (m$^3$/m$^3$)')
ax1.set_ylabel('Matric potential (kPa)')
ax1.set_ylim(y_min,15)

if yformat == 'positive':
    
    for axis in [ax1.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
else:
    ax1.set_yticks([0.1,1, 10])
    ax1.set_yticklabels(['-0.1','-1', '-10'])

ax1.set_xlim(-0.08,0.65)
ax1.text(0.08,0.88,'0\u20135 cm',transform=ax1.transAxes, fontsize=fontsize)

ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)

ax1.set_position(sppos[0]) 


ax2 = fig.add_subplot(1,3,2)

k = 1
a = '20_25'
i = 5

 
afp_drain = afp_2_6
press_drain = press_2_6
afp_imb = afp_2_6_w
press_imb = np.copy(press_2_6_w)
press_imb[press_imb>p_max] = p_max
   
theta = totpor[i+7] - thetas[i+7,:]

#f2 = interp1d(np.flip(theta), np.flip(np.log10(phi)), kind='linear')
#theta_fit_spline = np.linspace(theta[-1], theta[0],10000)

ax2.plot(theta, phi, 's', c=colors[0], label = 'sample ' + str(i+1), markersize=3)

#ax2.plot(theta_fit_spline, 10**(f2(theta_fit_spline)), ':', c=colors[0])

ax2.plot(afp_drain, 0.001*press_drain, '-', c=colors[0], lw=lw,label='Drainage')
ax2.plot(afp_imb, 0.001*press_imb, '--', c=colors[0], lw=lw,label='Imbibition')
ax2.fill_between(afp_imb, 0.001*press_imb, 0.001*np.interp(afp_imb, afp_drain, press_drain), alpha=0.2)
 
ax2.set_yscale('log')
#ax2.legend(handletextpad=0.2)
ax2.set_xlabel('Air-filled porosity')# (m$^3$/m$^3$)')
#ax2.set_ylabel('Pressure (kPa)')
ax2.set_ylim(y_min,15)

if yformat == 'positive':
    
    for axis in [ax2.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
else:
    ax2.set_yticks([0.1,1, 10])
    ax2.set_yticklabels(['-0.1','-1', '-10'])


ax2.set_xlim(-0.03,0.27)
ax2.text(0.08,0.88,'20\u201325 cm',transform=ax2.transAxes, fontsize=fontsize)

ax2.text(subtpos[0], subtpos[1], subt[1], transform=ax2.transAxes)

ax2.set_position(sppos[1]) 


ax3 = fig.add_subplot(1,3,3)

k = 1
a = '40_45'
i = 4

 
afp_drain = afp_3_5
press_drain = press_3_5
afp_imb = afp_3_5_w
press_imb = np.copy(press_3_5_w)
press_imb[press_imb>p_max] = p_max
   
theta = totpor[i+14] - thetas[i+14,:]

#f2 = interp1d(np.flip(theta), np.flip(np.log10(phi)), kind='linear')
#theta_fit_spline = np.linspace(theta[-1], theta[0],10000)

ax3.plot(theta, phi, 's', c=colors[0], label = 'Measured', markersize=3)

#ax3.plot(theta_fit_spline, 10**(f2(theta_fit_spline)), ':', c=colors[0])

ax3.plot(afp_drain, 0.001*press_drain, '-', c=colors[0], lw=lw,label='Drainage')
ax3.plot(afp_imb, 0.001*press_imb, '--', c=colors[0], lw=lw,label='Imbibition')
ax3.fill_between(afp_imb, 0.001*press_imb, 0.001*np.interp(afp_imb, afp_drain, press_drain), alpha=0.2)
 
ax3.set_yscale('log')
ax3.legend(handletextpad=0.6, borderaxespad=0.5)
ax3.set_xlabel('Air-filled porosity')# (m$^3$/m$^3$)')
#ax3.set_ylabel('Pressure (kPa)')
ax3.set_ylim(y_min,15)

if yformat == 'positive':
    
    for axis in [ax3.yaxis]:
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)
else:
    ax3.set_yticks([0.1,1, 10])
    ax3.set_yticklabels(['-0.1','-1', '-10'])


ax3.set_xlim(-0.02,0.16)
ax3.text(0.08,0.88,'40\u201345 cm',transform=ax3.transAxes, fontsize=fontsize)

ax3.text(subtpos[0], subtpos[1], subt[2], transform=ax3.transAxes)

ax3.set_position(sppos[2]) 


ax1.yaxis.set_label_coords(-0.16, 0.5)

ax1.xaxis.set_label_coords(0.5, -0.16)
ax2.xaxis.set_label_coords(0.5, -0.16)
ax3.xaxis.set_label_coords(0.5, -0.16)

#plt.savefig('hysteresis.pdf')
