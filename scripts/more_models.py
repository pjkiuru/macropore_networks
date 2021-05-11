# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:11:07 2020

@author: pkiuru
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import cumsum
import scipy.stats as stats
from collections import namedtuple
import matplotlib.animation as anim 
import networkx as nx
import pdb
import csv

#Units: exchange: mol(m3/s)
#       exchange_coefficient: 1/s


#phys_air['pore.exchC'] = 0.95
#phys_air['pore.Camb'] = 0.0
#phys_air.add_model(propname='pore.exchange', model=my_models.atmospheric_exchange,
               #exchange_coefficient=phys_air['pore.exchC'], ambient_concentration=phys_air['pore.Camb'], X='pore.concentration')

def atmospheric_exchange(target, X, exchange_coefficient, ambient_concentration):
    C_pore = target['pore.concentration']
    exchange = exchange_coefficient * (ambient_concentration - C_pore)
    return exchange

#def small_

def osakuva(target, numb, cmap='gray'):
    #import matplotlib.pyplot as plt
    plt.figure(num=numb, figsize=(7, 7))
    plt.clf()
    plt.imshow(target, cmap=cmap, interpolation='nearest')
    plt.title('')
    #plt.axis('off')
    plt.tight_layout()
    
def osakuva_8bit(target, numb, cmap='gray'):
    plt.figure(num=numb, figsize=(7, 7))
    plt.clf()
    plt.imshow(target, cmap=cmap, interpolation='nearest', vmin=0, vmax=255)
    plt.tight_layout()
    
def osakuvat_z(target1, target2, target3, coord, numb, cmap='gray', cmap2=plt.cm.jet):
    #import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(25, 7),
                           sharex=True, sharey=True)
    #plt.figure(num=numb, figsize=(25, 7))
    #plt.clf()
    #plt.subplot(1,3,1)
    ax[0].imshow(target1[:, :, coord], cmap=cmap, interpolation='nearest')
    #ax[0].axis('off')
    #plt.subplot(1,3,2)
    ax[1].imshow(target2[:, :, coord], cmap=cmap2, interpolation='nearest')
    #ax[1].axis('off')
    #plt.subplot(1,3,3)
    ax[2].imshow(target3[:, :, coord], cmap=cmap2, interpolation='nearest')
    #ax[2].axis('off')
    fig.tight_layout()

def anim3D(target, numb, cmap='gray'):
    
    counte = 1
    ims = []
    
    fig = plt.figure(num=numb)
    fig.clf()
    imi = plt.imshow(target[:,:,0], cmap=cmap, interpolation='nearest', animated = True)
    ims.append([imi])
    
    while counte < np.shape(target)[2]:
    
        imi = plt.imshow(target[:,:,counte], cmap=cmap, interpolation='nearest', animated = True)
        plt.title(counte)
        ims.append([imi])
        counte = counte + 1
    
    anim.ArtistAnimation(fig, ims, interval=200, blit=True, repeat_delay=1000)

def hajonta(xs, ys, numb):
    #import matplotlib.pyplot as plt
    #import numpy as np
    plt.figure(num=numb, figsize=(7, 7))
    plt.clf()
    plt.scatter(xs, ys)
    plt.grid()
    plt.xlim([0,1.05*np.max(xs)])
    plt.ylim([0,1.05*np.max(ys)])
    
def viiva(ys, numb):
    #import matplotlib.pyplot as plt
    #import numpy as np
    plt.figure(num=numb, figsize=(7, 7))
    plt.clf()
    plt.plot(np.arange(len(ys)), ys,'b-')
    plt.ylim([0,1.05*np.max(ys)])

def plotti(xs, ys, numb):
    #import matplotlib.pyplot as plt
    #import numpy as np
    plt.figure(num=numb, figsize=(7, 7))
    plt.clf()
    plt.plot(xs, ys,'bo-')
    plt.xlim([0,1.05*np.max(xs)])
    plt.ylim([0,1.05*np.max(ys)])

def plotti_piste(xs, ys, numb):
    #import matplotlib.pyplot as plt
    #import numpy as np
    plt.figure(num=numb, figsize=(7, 7))
    plt.clf()
    plt.plot(xs, ys,'bo')
    plt.xlim([0,1.05*np.max(xs)])
    plt.ylim([0,1.05*np.max(ys)])
  
def viipale(target, zcoord, numb):
    #import matplotlib.pyplot as plt
    fig, ax = plt.subplots(num=numb, figsize=(8, 8))
    #reg = snow_out.regions.astype(float) - 1
    #reg[reg == -1] = np.nan
    target_slice = target[:, :, zcoord] - 1
    #mask = target_slice >= 0
    plt.imshow(target_slice.T)
    plt.tight_layout()
    
def viipale_vapaa(target, index, coord, numb):
    #import matplotlib.pyplot as plt
    fig, ax = plt.subplots(num=numb, figsize=(8, 8))
    plt.clf()
    #reg = snow_out.regions.astype(float) - 1
    #reg[reg == -1] = np.nan
    if index == 0:
        target_slice = target[coord, :, :] - 1
    elif index == 1:
        target_slice = target[:, coord, :] - 1
    elif index == 2:
        target_slice = target[:, :, coord] - 1
    #mask = target_slice >= 0
    plt.imshow(target_slice.T)
    plt.tight_layout()
      
def pylvas(target, numb):
   #import matplotlib.pyplot as plt
   plt.figure(num=numb, figsize=[7, 7])
   plt.clf()
   plt.hist(target)
   
def regionslaissi(target, coord, numb):
    #import matplotlib.pyplot as plt
    #import numpy as np
    nsize = np.ndim(target)
    fig, ax = plt.subplots(num=numb, figsize=(8, 8))
    plt.clf()
    reg = target.astype(float) - 1
    reg[reg == -1] = np.nan
    if nsize == 3:
        region_slice = target[:, :, coord] - 1
    else:
        region_slice = target - 1
    plt.imshow(region_slice.T)
    #pdb.set_trace()
   
def pot_fit(x, a, b):
    return a * x**(b)

def local_maxima(histogram, bins_center, num=70, rang=100):
    
    ww = np.int(np.size(histogram)/rang)
    
    avg_dir = np.zeros(ww+1)
    
    for k in range(ww):
            avg_dir[k] = np.average(histogram[rang*k:rang*(k+1)])
    avg_dir[-1] = np.average(histogram[rang*(k+1):])
    
    bins_av = rang // 2 + np.min(bins_center) + np.arange(0, np.max(bins_center)-np.min(bins_center)+1, rang) 
    
    #Fluctuation in the long tails results in negative gradients
    avg_dir_truncated = np.copy(avg_dir)
    avg_dir_truncated[avg_dir_truncated<1] = 0
    
    grad = np.diff(avg_dir_truncated)
    
    boo = np.zeros(np.size(grad)).astype('bool')
    
    for k in range(np.size(grad)):
        try:
            boo[k] = grad[k] < 0 and grad [k-1] < 0 and grad[k+1] > 0 and grad[k+2] > 0
        except:
            continue
    
    if np.sum(boo)==0:
        
        for k in range(np.size(grad)):
            try:
                boo[k] = grad[k] < 0  and grad[k+1] > 0
            except Exception:
                continue
        
    #pdb.set_trace()      
    loclower = np.where(avg_dir == np.max(avg_dir[0:np.int(np.where(boo)[0][0])]))
    lochigher = np.where(avg_dir == np.max(avg_dir[np.int(np.where(boo)[0][0])+1:]))
    
    lower = bins_av[loclower[0][0]]
    higher = bins_av[lochigher[0][0]]
    #pdb.set_trace()
    
    plt.figure(num=num,figsize=(7,7))
    plt.clf()
    plt.plot(bins_center, histogram, linewidth=2)
    plt.plot(bins_av, avg_dir, 'r*')
    plt.plot([lower,higher], [avg_dir[loclower[0][0]],avg_dir[lochigher[0][0]]], 'k*',markersize=10)
    
    return bins_av, avg_dir, loclower, lochigher, lower, higher

def parse_histogram(h, voxel_size=1):
    delta_x = h[1]
    P = h[0]
    temp = P*(delta_x[1:] - delta_x[:-1])
    C = cumsum(temp[-1::-1])[-1::-1]
    S = P*(delta_x[1:] - delta_x[:-1])
    bin_edges = delta_x * voxel_size
    bin_widths = (delta_x[1:] - delta_x[:-1]) * voxel_size
    bin_centers = ((delta_x[1:] + delta_x[:-1])/2) * voxel_size
    psd = namedtuple('histogram', ('pdf', 'cdf', 'relfreq',
                                   'bin_centers', 'bin_edges', 'bin_widths'))
    return psd(P, C, S, bin_centers, bin_edges, bin_widths)

def trim_small_clusters(im, size=1):
    r"""
    Remove isolated voxels or clusters smaller than a given size

    Parameters
    ----------
    im : ND-array
        The binary image from which voxels are to be removed
    size : scalar
        The threshold size of clusters to trim.  As clusters with this many
        voxels or fewer will be trimmed.  The default is 1 so only single
        voxels are removed.

    Returns
    -------
    im : ND-image
        A copy of ``im`` with clusters of voxels smaller than the given
        ``size`` removed.

    """
    import scipy as sp
    import scipy.ndimage as spim
    from skimage.morphology import ball, disk
        
    if im.ndim == 2:
        strel = disk(1)
    elif im.ndim == 3:
        strel = ball(1)
    else:
        raise Exception('Only 2D or 3D images are accepted')
    filtered_array = sp.copy(im)
    labels, N = spim.label(filtered_array, structure=strel)
    id_sizes = sp.array(spim.sum(im, labels, range(N + 1)))
    area_mask = (id_sizes <= size)
    filtered_array[area_mask[labels]] = 0
    return filtered_array

def read_params(filename):
    
    crop_coords = []
    rot_params = []
    scale_params = []
    
    with open(filename, newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1   
            else:
                crop_coords.append([int(row[2]),int(row[3]),int(row[4]),int(row[5]),int(row[6]),int(row[7])])
                rot_params.append([float(row[8]),float(row[9])])
                scale_params.append([int(row[11]),int(row[12])])
                line_count += 1
    
    crop_coords = np.asarray(crop_coords)
    rot_params = np.asarray(rot_params)
    scale_params = np.asarray(scale_params)
    
    return crop_coords, rot_params, scale_params

def load_network_volumes(filename):
    
    #import os
    import openpnm as op
    
    topfolder = '../IO/pnms/'
    suffix = '_600.pnm'

    
    ws = op.Workspace()
    ws.clear()
    ws.load_workspace(filename=topfolder+filename+suffix, overwrite=True)
    
    proj5 = ws['image']
    vols = proj5['net_01']['pore.volume']
    internal = proj5['net_01']['pore.internal']
    thrdiams = 1000 * proj5['net_01']['throat.diameter']
    thrdiams_int = 1000 * proj5['net_01']['throat.diameter'][proj5['net_01']['throat.internal']]
    
    proj5w = ws['whole_network']
    vols_all = proj5w['net_01']['pore.volume']
    internal_all = proj5w['net_01']['pore.internal']
    thrdiams_all = 1000 * proj5w['net_01']['throat.diameter']

    return vols, vols_all, internal, internal_all, thrdiams, thrdiams_int, thrdiams_all

def load_percolation_curve(filename):
    
    topfolder = '../IO/results/perc_'
    suffix = '.npz'
    
    dummy = np.load(topfolder+filename+suffix) 
    
    return dummy['press'], dummy['afporos'] 

def load_water_percolation_curve(filename):
    
    topfolder = '../IO/results/perc_water_'
    suffix = '.npz'
    
    dummy = np.load(topfolder+filename+suffix) 
    
    return dummy['press_w'], dummy['afporos_w'] 

def load_volumefraction(filename):
    
    topfolder = '../IO/results/volfrac_'
    suffix = '_wide.npz'
    
    dummy = np.load(topfolder+filename+suffix) 
    
    return dummy['prpoints'], dummy['afporos'], dummy['volfraction'] 

def load_volumefraction_imb(filename):
    
    topfolder = 'C:/Users/pkiuru/Pictures/metnet_ct_reco/volfrac_'
    suffix = '_imb.npz'
    
    dummy = np.load(topfolder+filename+suffix) 
    
    return dummy['prpoints'], dummy['afporos'], dummy['volfraction'] 
    
def volumehist(vols, internal, nbins, w=True):
    
    if w:
        weights=vols[internal]
    else:
        weights = None
    
    return np.histogram(np.log10(vols[internal]), bins=nbins, density=False, weights=weights)

def diamhist(vols, internal, nbins, w=True):
    
    if w:
        weights=vols[internal]
    else:
        weights = None
    
    diams = equiv_diam(volfact=1, vols=vols)
    
    return np.histogram(np.log10(diams[internal]), bins=nbins, density=False, weights=weights)

def equiv_diam(volfact, vols):
    
    return 2 * (0.75 * volfact * vols / np.pi)**(1/3)

def diamdistr_2(diamfit, diamhist, internal):
    
    fig = plt.figure(figsize=(13,5))
    fig.clf()
    ax1 = fig.add_subplot(1,3,1)
    xax = lognorm_xaxis(diamfit)
    ax1.semilogx(10**xax, stats.lognorm.pdf(xax, diamfit[0], diamfit[1], diamfit[2]))
    ax2 = fig.add_subplot(1,3,2)
    ax2.bar(np.array(10**(diamhist[1][0:-1])), np.array(diamhist[0]), width = np.diff(10**(diamhist[1])), align='edge')
    ax2.set_xscale('log')
    ax3 = fig.add_subplot(1,3,3)
    ax3.semilogx(np.array(10**(diamhist[1][0:-1])), np.array(diamhist[0])/np.sum(internal))
    ax3.set_xscale('log')

def diamdistr_3(diamfit, diamhist, internal):

    fig = plt.figure(figsize=(13,5))
    fig.clf()
    ax1 = fig.add_subplot(1,3,1)
    xax = lognorm_xaxis(diamfit)
    cdf = stats.lognorm.cdf(xax, diamfit[0], loc=diamfit[1], scale=diamfit[2])
    ax1.semilogx(10**xax[0:-1], np.diff(cdf)*np.sum(internal))
    ax3= fig.add_subplot(1,3,2)
    ax3.bar(np.array(10**xax[0:-1]), np.diff(cdf)*np.sum(internal), width = np.diff(10**xax), align='edge')
    ax3.set_xscale('log')
    

def lognorm_xaxis(fit, xmin=0.0001, xmax=0.999, n=200):

    xax = np.linspace(stats.lognorm.ppf(xmin, fit[0],fit[1],fit[2]), stats.lognorm.ppf(xmax, fit[0],fit[1],fit[2]), n)

    return xax

def lognorm_pdf(fit, xmin=0.0001, xmax=0.999, n=200): 
    
    return stats.lognorm.pdf(lognorm_xaxis(fit, xmin, xmax, n), fit[0], fit[1], fit[2])

def gamma_xaxis(fit, xmin=0.0001, xmax=0.999, n=200):

    xax = np.linspace(stats.gamma.ppf(xmin, fit[0],fit[1],fit[2]), stats.gamma.ppf(xmax, fit[0],fit[1],fit[2]), n)

    return xax

def gamma_pdf(fit, xmin=0.0001, xmax=0.999, n=200): 
    
    return stats.gamma.pdf(gamma_xaxis(fit, xmin, xmax, n), fit[0], fit[1], fit[2])

def weibull_xaxis(fit, xmin=0.0001, xmax=0.999, n=200):

    xax = np.linspace(stats.weibull_min.ppf(xmin, fit[0],fit[1],fit[2]), stats.weibull_min.ppf(xmax, fit[0],fit[1],fit[2]), n)

    return xax

def weibull_pdf(fit, xmin=0.0001, xmax=0.999, n=200): 
    
    return stats.weibull_min.pdf(weibull_xaxis(fit, xmin, xmax, n), fit[0], fit[1], fit[2])

def expon_xaxis(fit, xmin=0.0001, xmax=0.999, n=200):

    xax = np.linspace(stats.expon.ppf(xmin, fit[0],fit[1]), stats.expon.ppf(xmax, fit[0],fit[1]), n)

    return xax

def expon_pdf(fit, xmin=0.0001, xmax=0.999, n=200): 
    
    return stats.expon.pdf(expon_xaxis(fit, xmin, xmax, n), fit[0], fit[1])

def vol_cumsum(vols,internal):
    
    tup = namedtuple('cs', field_names=['x', 'cumsum'])
    
    diams = equiv_diam(volfact=1, vols=vols[internal])
    uniq_diam = np.unique(diams)
    intvols = vols[internal]
    cumsum_a = np.zeros(np.size(uniq_diam))
    for v, w in enumerate(uniq_diam):
        cumsum_a[v] = np.sum(intvols[diams <= uniq_diam[v]])
    
    tup.x = uniq_diam
    tup.cumsum = cumsum_a
    
    return tup
    
def tortuosity(pnx, space=25):
    
    print('Calculating tortuosity')
    
    coords = nx.get_node_attributes(pnx,'coords')
    #lengths = nx.get_edge_attributes(pnx,'total_length')
    vcoords = []
    for k, kk in enumerate(coords):
        vcoords.append(coords[k][0])
    
    #Tähänkin näyttää olevan useampi keino. (Yksinkertaisemmalla filter-funktion
    #käytöllä saa sitten haluttujen solujen arvot kun nyt tulee indeksit.)
    # https://www.geeksforgeeks.org/python-ways-to-find-indices-of-value-in-list/
    top = list(filter(lambda x: vcoords[x] == min(vcoords), range(len(vcoords))))
    print(len(top))
    bottom = list(filter(lambda x: vcoords[x] == max(vcoords), range(len(vcoords))))
    print(len(bottom))
    
    le = []
    for i in range(0, len(top),space):
        for j in range(0, len(bottom),space):
            le.append(nx.dijkstra_path_length(pnx,top[i],bottom[j],'conduitlength'))
        if i%10==0:
            print(i)
    path_lengths = np.asarray(le)/(max(vcoords)-min(vcoords))
    
    return path_lengths

def tortuosity_n(pnx, space=25, axis=0):
    
    print('Calculating tortuosity in a specified direction:')
    print(str(axis))
    
    coords = nx.get_node_attributes(pnx,'coords')
    #lengths = nx.get_edge_attributes(pnx,'total_length')
    vcoords = []
    for k, kk in enumerate(coords):
        vcoords.append(coords[k][axis])
    
    #Tähänkin näyttää olevan useampi keino. (Yksinkertaisemmalla filter-funktion
    #käytöllä saa sitten haluttujen solujen arvot kun nyt tulee indeksit.)
    # https://www.geeksforgeeks.org/python-ways-to-find-indices-of-value-in-list/
    top = list(filter(lambda x: vcoords[x] == min(vcoords), range(len(vcoords))))
    print(len(top))
    bottom = list(filter(lambda x: vcoords[x] == max(vcoords), range(len(vcoords))))
    print(len(bottom))
    
    le = []
    for i in range(0, len(top),space):
        for j in range(0, len(bottom),space):
            le.append(nx.dijkstra_path_length(pnx,top[i],bottom[j],'conduitlength'))
        if i%10==0:
            print(i)
    path_lengths = np.asarray(le)/(max(vcoords)-min(vcoords))
    
    return path_lengths

def betweenness_centrality_subset(pnx,norm=True,space=25):
    
    internal_nodes = [x for x,y in pnx.nodes(data=True) if y['internal']==True]
    
    print('Calculating top-bottom betweenness centrality')
    
    coords = nx.get_node_attributes(pnx,'coords')
    vcoords = []
    for k, kk in enumerate(coords):
        vcoords.append(coords[k][0])
    
    top = list(filter(lambda x: vcoords[x] == min(vcoords), range(len(vcoords))))
    print(len(top))
    bottom = list(filter(lambda x: vcoords[x] == max(vcoords), range(len(vcoords))))
    print(len(bottom))
    
    poredata = np.asarray(list([len(internal_nodes), len(top), len(bottom)]))
    
    bcs = nx.algorithms.centrality.betweenness_centrality_subset(pnx, sources=top, targets=bottom, normalized=norm, weight='conduitlength')
    
    return np.asarray(list(bcs.values())), poredata

def betweenness_centrality_subset_n(pnx,norm=True,space=25,axis=0):
    
    internal_nodes = [x for x,y in pnx.nodes(data=True) if y['internal']==True]
    
    print('Calculating throughflow betweenness centrality')
    
    coords = nx.get_node_attributes(pnx,'coords')
    vcoords = []
    for k, kk in enumerate(coords):
        vcoords.append(coords[k][axis])
    
    top = list(filter(lambda x: vcoords[x] == min(vcoords), range(len(vcoords))))
    print(len(top))
    bottom = list(filter(lambda x: vcoords[x] == max(vcoords), range(len(vcoords))))
    print(len(bottom))
    
    poredata = np.asarray(list([len(internal_nodes), len(top), len(bottom)]))
    
    bcs = nx.algorithms.centrality.betweenness_centrality_subset(pnx, sources=top, targets=bottom, normalized=norm, weight='conduitlength')
    
    return np.asarray(list(bcs.values())), poredata

def probability_plot(fignum, data_1, data_2, data_3):

    fig = plt.figure(num=fignum)
    fig.set_size_inches(10,4)
    plt.clf()
    ax1 = fig.add_subplot(131)
    stats.probplot(data_1, plot= plt, rvalue= True)
    ax1.set_title(str(fignum))
    ax2 = fig.add_subplot(132)
    stats.probplot(data_2, plot= plt, rvalue= True)
    ax2.set_title(str(fignum))
    ax2.set_ylabel('')
    ax3 = fig.add_subplot(133)
    stats.probplot(data_3, plot= plt, rvalue= True)
    ax3.set_title(str(fignum))
    ax3.set_ylabel('')
    plt.tight_layout()
    
def probability_plot_ols(fignum, data_1):

    fig = plt.figure(num=fignum)
    fig.set_size_inches(4,4)
    plt.clf()
    ax1 = fig.add_subplot(111)
    stats.probplot(data_1, plot= plt, rvalue= True)
    ax1.set_title(str(fignum))
    plt.tight_layout()

def prob_plot_r2(diams, internal, params):   
    
    #params = stats.gamma.fit(diams[internal]))
    a = stats.probplot(diams[internal], dist = stats.lognorm, sparams=params)
    b = a[1][2]**2
    
    return b
    
def to_networkx(network, geometry):
    r"""
    Write OpenPNM Network to a NetworkX object.
    
    Parameters
    ----------
    network : OpenPNM Network Object
        The OpenPNM Network to be converted to a NetworkX object
    
    Returns
    -------
    A NetworkX object with all pore/throat properties attached to it
    """

    G = nx.Graph()
    
    # Extracting node list and connectivity matrix from Network
    nodes = map(int, network.Ps)
    conns = network['throat.conns']
    
    # Explicitly add nodes and connectivity matrix
    G.add_nodes_from(nodes)
    G.add_edges_from(conns)
    
    # Attach Network properties to G
    for prop in network.props(deep=False) + network.labels():
        if 'pore.' in prop:
            if len(network[prop].shape) > 1:
                val = {i: list(network[prop][i]) for i in network.Ps}
            else:
                val = {i: network[prop][i] for i in network.Ps}
            nx.set_node_attributes(G, name=prop[5:], values=val)
        if 'throat.' in prop:
            val = {tuple(conn): network[prop][i] for i, conn
                   in enumerate(conns)}
            nx.set_edge_attributes(G, name=prop[7:], values=val)
    
    for prop in geometry.props():
        if 'pore.' in prop:
            if len(geometry[prop].shape) > 1:
                val = {i: list(geometry[prop][i]) for i in network.Ps}
            else:
                val = {i: geometry[prop][i] for i in network.Ps}
            nx.set_node_attributes(G, name=prop[5:], values=val)
        if 'throat.' in prop:
            val = {tuple(conn): geometry[prop][i] for i, conn
                   in enumerate(conns)}
            nx.set_edge_attributes(G, name=prop[7:], values=val)
    
    return G
