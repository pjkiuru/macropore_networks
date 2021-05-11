# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:41:17 2020

@author: pkiuru

Script for drawing result graphs for the development of air volume fraction
in the cubical-domain pore spaces.

"""

import numpy as np
import matplotlib.pyplot as plt

import my_models as mm
import pandas as pd

from matplotlib import rcParams

fontname = 'Arial'
fontsize = 8

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

rcParams['font.size'] = fontsize

rcParams['axes.titlesize'] = fontsize + 2
rcParams['axes.labelsize'] = fontsize
rcParams['axes.titleweight'] = 'normal'

rcParams['xtick.labelsize'] = fontsize
rcParams['ytick.labelsize'] = fontsize

#rcParams['legend.handlelength'] = 1.0

#Pie charts for a single depth
draw_figure = 0

#Generates 2D tables for pasting to Excel-file
metrics = 0

#Box plots
boxplots = 0

if boxplots:
    #Box plots for imbibition
    boxplots_imb = 0

compound_bp = 1

if draw_figure:

    colors = ['black', 'dimgray', 'gray','darkgray', 'silver','lightgray','gainsboro','whitesmoke']
    
    # 0.30  0.50  0.75  1.00 1.50 2.00 2.50 3.00
    mask = np.array([1, 1, 0, 1, 0, 1, 1, 1]).astype('bool')
    
    a = '0_5'
    
    if a == '0_5':
        fignum = 970
    elif a == '20_25':
        fignum = 971
    elif a == '40_45':
        fignum = 972
        
    fig = plt.figure(num=fignum)
    fig.set_size_inches(12,6)
    plt.clf()
    
    for w in range(1,8):
        
        ax = fig.add_subplot(2,4,w)
        try:
            
            points_1, afp_1, vfcn_1 = mm.load_volumefraction_imb(a+'_'+str(w))
            
            points_sub = points_1[mask]
            afp_sub = afp_1[mask]
            vfcn_sub = vfcn_1[mask]
            
            fractions = np.diff(np.concatenate([np.array([0.0]), vfcn_sub]))
            
            labels = ["%1.1f" % number for number in points_sub]
            
            
            patches, texts = ax.pie(fractions, colors=colors, startangle=90, radius=0.8)
            pielabels = ['{0} kPa - {1:1.2f}%'.format(i,j) for i,j in zip(np.asarray(labels), 100.*vfcn_sub)]
            
            ax.legend(patches, pielabels, loc='upper center', bbox_to_anchor=(-0.1, 1.),
                       fontsize=8)
        except Exception:
            
            pass
        
    fig.tight_layout()


if metrics:
    
    #Generate a 2D array of result for saving to Excel
    
    if True:
        
        a = '0_5'
        array_of_arrays = np.ndarray(shape = (7,21))
        
        for w in range(1,8):
        
           points_print, afp_print, vfcn_print = mm.load_volumefraction(a+'_'+str(w))
           
           array_of_arrays[w-1,:] = 100*vfcn_print
       
   
#Boxplots
   
if boxplots:
    
    
   bbdim = [0.21, 0.36]

   bppos = [
            [0.07, 0.565, bbdim[0],bbdim[1]],
            [0.30, 0.565, bbdim[0],bbdim[1]],
            [0.53, 0.565, bbdim[0],bbdim[1]],
            [0.76, 0.565, bbdim[0],bbdim[1]],
            [0.07, 0.095, bbdim[0],bbdim[1]],
            [0.30, 0.095, bbdim[0],bbdim[1]],
            [0.53, 0.095, bbdim[0],bbdim[1]],
            [0.76, 0.095, bbdim[0],bbdim[1]],
            ]
   
   subtpos = [0.01,1.06]
   
   subt = [('(a)'), ('(b)'), ('(c)'), ('(d)'),
           ('(e)'), ('(f)'), ('(g)'), ('(h)')]
   
   medianprops = dict(linestyle='-', linewidth=1, color='blue')
   flierprops = dict(marker='o', markersize=4)#, linestyle='none')
   
   if boxplots_imb:
       df=pd.read_excel('volfracs.xlsx', sheet_name='Sheet5', header=0)
   else:
       df=pd.read_excel('volfracs.xlsx', sheet_name='Sheet1', header=0)
   
   fig = plt.figure(num=310)
   fig.set_size_inches(5.5,3)
   plt.clf()
      
   for i in range(1,9):
        
       
       ax = fig.add_subplot(2,4,i)
       
       data = [df[df.columns[i]][df['Depth'] == '0-5'],
               df[df.columns[i]][df['Depth'] == '20-25'].dropna(),
               df[df.columns[i]][df['Depth'] == '40-45']]
       
       data = 0.01 * np.asarray(data)
       
       ax.boxplot(data, labels= ['0\u20135', '20\u201325', '40\u201345'],
                  showmeans= False, medianprops=medianprops,
                  flierprops=flierprops)
       
       ax.set_ylim([-0.1, 1.1])#([-10, 110])
       if i >= 5:
           ax.set_xlabel("Depth (cm)")   
           ax.xaxis.set_label_coords(0.5, -0.15)
       
       if i ==1:
           ax.set_ylabel("Volume fraction of connected pore space")# (%)")
       
       #ax.set_title(subt[i-1],loc='left')
       
       ax.text(subtpos[0], subtpos[1], subt[i-1], transform=ax.transAxes)
       
       ax.text(0.5, 0.87, df.columns[i], transform=ax.transAxes,
               horizontalalignment='center', fontsize=fontsize)
       
       ax.xaxis.set_tick_params(length=0)
       
       if i != 1 and i != 5:
           ax.yaxis.set_ticklabels([])
       
       if i == 1 or i == 5: 
           ax.yaxis.set_label_coords(-0.2, -0.15)
 
       
       ax.set_position(bppos[i-1])


#Compound boxplots
   
if compound_bp:
    
    
   bbdim = [0.21, 0.36]

   bppos = [
            [0.07, 0.565, bbdim[0],bbdim[1]],
            [0.30, 0.565, bbdim[0],bbdim[1]],
            [0.53, 0.565, bbdim[0],bbdim[1]],
            [0.76, 0.565, bbdim[0],bbdim[1]],
            [0.07, 0.095, bbdim[0],bbdim[1]],
            [0.30, 0.095, bbdim[0],bbdim[1]],
            [0.53, 0.095, bbdim[0],bbdim[1]],
            [0.76, 0.095, bbdim[0],bbdim[1]],
            ]
   
   subtpos = [0.01,1.06]
   
   subt = [('(a)'), ('(b)'), ('(c)'), ('(d)'),
           ('(e)'), ('(f)'), ('(g)'), ('(h)')]
   
   pos_left = [0.65, 2.65, 4.65]
   pos_right = [1.35, 3.35, 5.35]
   xticks = [1, 3, 5]
   xlabels= ['0\u20135', '20\u201325', '40\u201345']
   
   boxprops = dict(linestyle='-', linewidth=1)
   medianprops = dict(linestyle='-', linewidth=1)
   flierprops = dict(marker='o', markersize=3)#, linestyle='none')
   
   df_d = pd.read_excel('volfracs.xlsx', sheet_name='Sheet1', header=0)
   df_i = pd.read_excel('volfracs.xlsx', sheet_name='Sheet5', header=0)
   
   def set_box_color(bp, color, color2):
       plt.setp(bp['boxes'], color=color)
       plt.setp(bp['boxes'], facecolor=color)
       plt.setp(bp['whiskers'], color=color2)
       plt.setp(bp['caps'], color=color2)
       plt.setp(bp['medians'], color=color2)
       plt.setp(bp['fliers'], markeredgecolor=color2)
   
   fig = plt.figure(num=310)
   fig.set_size_inches(5.5,3)
   plt.clf()
      
   for i in range(1,9):
        
       
       ax = fig.add_subplot(2,4,i)
       
       data_d = [df_d[df_d.columns[i]][df_d['Depth'] == '0-5'],
               df_d[df_d.columns[i]][df_d['Depth'] == '20-25'].dropna(),
               df_d[df_d.columns[i]][df_d['Depth'] == '40-45']]
       
       data_d = 0.01 * np.asarray(data_d)
       
       data_i = [df_i[df_i.columns[i]][df_i['Depth'] == '0-5'],
               df_i[df_i.columns[i]][df_i['Depth'] == '20-25'].dropna(),
               df_i[df_i.columns[i]][df_i['Depth'] == '40-45']]
       
       data_i = 0.01 * np.asarray(data_i)
       
       bpl_d = ax.boxplot(data_d, positions = pos_left, widths = 0.40,
                  showmeans= False, patch_artist=True, medianprops=medianprops,
                  flierprops=flierprops, boxprops=boxprops)
       
       bpl_i = ax.boxplot(data_i, positions = pos_right, widths = 0.40,
                  showmeans= False, patch_artist=True, medianprops=medianprops,
                  flierprops=flierprops, boxprops=boxprops)
       
       # ['black', 'dimgray', 'gray','darkgray', 'silver','lightgray','gainsboro','whitesmoke']
       set_box_color(bpl_d, 'gray', 'k')
       set_box_color(bpl_i, 'silver', 'k')
       
       ax.set_xticks(xticks)
       ax.set_xticklabels(xlabels)
       
       ax.set_ylim([-0.12, 1.12])
       ax.set_xlim([0.0, 6.0])
       
       if i >= 5:
           ax.set_xlabel("Depth (cm)")   
           ax.xaxis.set_label_coords(0.5, -0.15)
       
       if i == 1:
           ax.set_ylabel("Volume fraction of connected pore space")# (%)")
           ax.legend([bpl_d["boxes"][0], bpl_i["boxes"][0]],
                     ['Drainage', 'Imbibition'], handlelength=0.6, 
                     loc='upper left', bbox_to_anchor=(0.30, 0.85))
       
       #ax.set_title(subt[i-1],loc='left')
       
       ax.text(subtpos[0], subtpos[1], subt[i-1], transform=ax.transAxes)
       
       ax.text(0.5, 0.87, df_d.columns[i], transform=ax.transAxes,
               horizontalalignment='center', fontsize=fontsize)
       
       ax.xaxis.set_tick_params(length=0)
       
       if i != 1 and i != 5:
           ax.yaxis.set_ticklabels([])
       
       if i == 1 or i == 5: 
           ax.yaxis.set_label_coords(-0.2, -0.15)
 
       
       ax.set_position(bppos[i-1])
  
       #plt.savefig('connvolfrac_di.pdf')