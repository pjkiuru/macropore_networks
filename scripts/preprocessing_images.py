# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 13:04:28 2020

@author: pkiuru

Script fof generating figures of tomography image preprocessing

"""

import numpy as np
import matplotlib.pyplot as plt
import my_models as mm
import porespy as ps

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

subtpos = [0.01,1.06]
subt = ['(a)', '(b)', '(c)', '(d)']

# Filtered + binary image, 40_45_4 (scan = 19)

if True:

    print('Minimum and maximum intensities in slice:')
    
    im_filt = filtered[200:800,475,200:800]
    bin_img_slice = bin_img_nofloat[200:800,475,200:800]
    
    if False:
        np.random.seed(42)
        psnet = ps.networks.snow(bin_img_nofloat[200:800,200:800,200:800], voxel_size = 50e-6, boundary_faces=['left','right'])
        regions_slice = psnet.regions[:,275,:].astype('float')
        regions_slice[regions_slice==0] = np.NaN
         
    
    #im_filt = filtered[:,275,:]
    #bin_img_slice = bin_img_nofloat[:,275,:]
    #np.random.seed(42)
    #psnet = ps.networks.snow(bin_img_nofloat, voxel_size = 50e-6, boundary_faces=['left','right'])
    #regions_slice = psnet.regions[:,275,:].astype('float')
    #regions_slice[regions_slice==0] = np.NaN
    
    #------------------------------
    
    hist_slice, bins_center_slice = exposure.histogram(im_filt)
    print(np.min(im_filt),np.max(im_filt))
    
    fig = plt.figure(num=2)
    fig.set_size_inches(6,1.8)
    plt.clf()
    
    ax1 = fig.add_subplot(1,4,1)
    ax1.imshow(im_filt, vmin=0, vmax=255, cmap='gray', interpolation='nearest')
    #plt.axis('off')
    
    #ax1.set_title('a)', loc='left', fontsize = fontsize)
    ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)
    #ax1.patch.set_edgecolor('black')  
    #ax1.patch.set_linewidth('0.5') 
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks([]) 
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])
    
    ax1.spines['left'].set_linewidth(0.5)
    ax1.spines['right'].set_linewidth(0.5)
    ax1.spines['top'].set_linewidth(0.5)
    ax1.spines['bottom'].set_linewidth(0.5)
    
    ax2 = fig.add_subplot(1,4,3)
    ax2.imshow(~bin_img_slice, cmap='gray', interpolation='nearest')
    #plt.axis('off')
    #ax2.set_title('b', loc='left', fontsize = fontsize)
    ax2.text(subtpos[0], subtpos[1], subt[2], transform=ax2.transAxes)
    ax2.patch.set_edgecolor('black')  
    ax2.patch.set_linewidth('0.5') 
    ax2.xaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks([]) 
    ax2.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticklabels([])
    
    ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['right'].set_linewidth(0.5)
    ax2.spines['top'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)
    
    ax3 = fig.add_subplot(1,4,4)
    ax3.imshow(regions_slice, cmap=plt.cm.jet)
    #plt.axis('off')
    #ax3.set_title('c', loc='left', fontsize = fontsize)
    ax3.text(subtpos[0], subtpos[1], subt[3], transform=ax3.transAxes)
    ax3.patch.set_edgecolor('black')  
    ax3.patch.set_linewidth('0.5') 
    ax3.xaxis.set_ticks_position('none')
    ax3.yaxis.set_ticks([]) 
    ax3.xaxis.set_ticklabels([])
    ax3.yaxis.set_ticklabels([])
    
    ax3.spines['left'].set_linewidth(0.5)
    ax3.spines['right'].set_linewidth(0.5)
    ax3.spines['top'].set_linewidth(0.5)
    ax3.spines['bottom'].set_linewidth(0.5)
    
    ax4 = fig.add_subplot(1,4,2)
    ax4.plot(bins_centerf, 1e-6*histf, c = 'k', lw=1)
    ax4.set_xlabel('Grayscale intensity', fontname=fontname, fontsize = fontsize, labelpad=2)
    #ax4.set_ylabel(r'Number of voxels ($ \times 10^6$)', fontname=fontname, fontsize = fontsize, labelpad=2)
    ax4.set_ylabel('Number of voxels (\u00D710$^6$)', fontsize = fontsize, labelpad=0)
    
    ax4.set_xlim([0,255])
    
    ax4.tick_params(labelsize=fontsize)
    #ax4.text(10,20,'d')
    #ax4.set_title('d', loc='left', fontsize = fontsize)
    ax4.text(subtpos[0], subtpos[1], subt[1], transform=ax4.transAxes)
    
    #ax1.set_position([0.02, 0.10, 0.21, 0.81])
    #ax4.set_position([0.30, 0.20, 0.20, 0.66])
    #ax2.set_position([0.52, 0.10, 0.21, 0.81])
    #ax3.set_position([0.75, 0.10, 0.21, 0.81])
    
    ax1.set_position([0.02, 0.11, 0.20, 0.81])
    ax4.set_position([0.30, 0.19, 0.20, 0.66])
    ax2.set_position([0.535, 0.11, 0.20, 0.81])
    ax3.set_position([0.77, 0.11, 0.20, 0.81])
    
    #plt.savefig("preprocess.pdf")
    
    #ax4 = fig.add_subplot(2,3,6)
    #ax4.plot(bins_center_slice, hist_slice, c = 'k', lw=1.5)
    
    #plt.tight_layout()
    



#ps.io.to_vtk(bin_img_nofloat[200:800,200:800,200:800], 
#             path='C:/Users/pkiuru/Pictures/metnet_ct_reco/binary_40_45_4_600', vox=True)



#op.io.VTK.save(network=pn_5, filename='C:/Users/pkiuru/Pictures/metnet_ct_reco/pn_40_45_4_600')























