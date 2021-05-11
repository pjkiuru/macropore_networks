# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:50:10 2020

@author: pkiuru

Script for calculating network metrics for the cubical domain networks

Required input:
    OpenPNM network object from a file or from memory    
    NOTE: If a network is created during the same session and not loaded,
    perc_diff_simulations.py with domain_switch = cube, percolation_switch = True and
    diffusion_switch = True must be run in advance!

Output:
    Network metrics calculation results (printed in the console)


"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import networkx as nx
import openpnm as op
import more_models as mm

# True: Load the OpenPNM workspace (network object) from file; False: use the 
# workspace / network object in memory 
load_switch = True
ws_file = '../IO/pnms/40_45_1_600.pnm'

#Tortuosity and top-bottom betweenness centrality uses the boundary pores at the
#top and bottom of the network

# Calculation of tortuosity
tortuosity_switch = True

#Calculation of top-bottom betweenness centrality
tb_betw_centr_switch = True

#Other network metrics use only the internal network (boundary pores excluded)

#Calculation of pore coordination number
degree_switch = True

#Calculation of pore clustering coefficient
clust_switch = True

#Calculation of pore centralities
centrality_switch = True

#Reduction of nodes included in tortuosity / centrality calculations by
#a factor 'space'
space = 1

#Calculation of horizontal tortuosity (horizontal boundary pores have to be
#determined)
horiz = False

# -----------------------------------------------------------------------------

# Load network object from file
if load_switch:
    
    ws = op.Workspace()
    ws.clear()
    ws.load_workspace(filename=ws_file, overwrite=True)
    pn_5 = ws['image']['net_01']

# Create a NetworkX object from the total pore network and from the network of
# internal pores

pn_5['throat.conduitlength'] = pn_5['throat.conduit_lengths.pore1'] + \
pn_5['throat.conduit_lengths.pore2'] +pn_5['throat.conduit_lengths.throat']
pnx = op.io.NetworkX.to_networkx(pn_5)

#Remove the boundary pores from the NetworkX object
internal_nodes = [x for x,y in pnx.nodes(data=True) if y['internal']==True]        
pnx_int = pnx.subgraph(nodes=internal_nodes)

 
if tortuosity_switch:

    path_lengths = mm.tortuosity(pnx,space)
    
    if horiz:
        path_lengths_h1 = mm.tortuosity_n(pnx,space = space, axis = 1)
        path_lengths_h2 = mm.tortuosity_n(pnx,space = space, axis = 2)
        
    his = mm.parse_histogram(np.histogram(path_lengths, bins=25, density=True))

    fit_alpha, fit_loc, fit_beta = stats.gamma.fit(path_lengths)
        
    xax = np.linspace(stats.gamma.ppf(0.005, fit_alpha, fit_loc, fit_beta),
                          stats.gamma.ppf(0.995, fit_alpha, fit_loc, fit_beta), 100)
    
    print('')
    print('Vertical tortuosity:')
    print(np.mean(path_lengths))
    #print(stats.gamma.stats(fit_alpha, fit_loc, fit_beta,'mvsk'))        
            
    plt.figure(num=572,figsize=(5,5))
    plt.clf()
    plt.subplot(1,1,1)
    plt.bar(his.bin_centers, his.pdf, label='Calculated', width = his.bin_centers[1]-his.bin_centers[0])
    plt.plot(xax, stats.gamma.pdf(xax, fit_alpha, fit_loc, fit_beta), c='sandybrown', label='Fitted gamma distribution')
    plt.xlabel('Tortuosity')
    plt.ylabel('Probability density')
    plt.legend()
        
    if horiz:
        
        print('')
        print('H_1 tortuosity:')
        fit_alpha, fit_loc, fit_beta = stats.gamma.fit(path_lengths_h1)
        print(stats.gamma.stats(fit_alpha, fit_loc, fit_beta,'mvsk'))
        
        
        print('')
        print('H_2 tortuosity:')
        fit_alpha, fit_loc, fit_beta = stats.gamma.fit(path_lengths_h2)
        print(stats.gamma.stats(fit_alpha, fit_loc, fit_beta,'mvsk'))
        
        print('')
        print('Average h1 tortuosity')
        print(np.mean(path_lengths_h1))
        print('Average h2 tortuosity')
        print(np.mean(path_lengths_h2))
        print('Average total h tortuosity')
        print(np.mean(np.concatenate([path_lengths_h1, path_lengths_h2])))


if tb_betw_centr_switch:

    bcs, nodecount = mm.betweenness_centrality_subset(pnx,False,1)

    #nodecount: internal pores, top pores, bottom pores
    n_internal, n_top, n_bottom = nodecount[0], nodecount[1], nodecount[2]
    n_all = np.sum(nodecount)
    
    #https://github.com/networkx/networkx/issues/3481
    print('Average betweenness centrality for subset')
    print(2*np.mean(bcs))
    print('Average normalized betweenness centrality for subset')
    print(2*np.mean(bcs)/((n_all-1)*(n_all-2)))
    print('Average betweenness centrality for subset normalized for throughflow')
    print(2*np.mean(bcs[0:n_internal])/(n_top*n_bottom))

    if horiz:
    
        bcs_h1, nodecount_h1 = mm.betweenness_centrality_subset_n(pnx,False, space=1, axis=1)
        n_internal_h1, n_top_h1, n_bottom_h1 = nodecount_h1[0], nodecount_h1[1], nodecount_h1[2]
        n_all_h1 = np.sum(nodecount_h1)
        
        print('Average betweenness centrality for subset h1')
        print(2*np.mean(bcs_h1))
        print('Average normalized betweenness centrality for subset h1')
        print(2*np.mean(bcs_h1)/((n_all_h1-1)*(n_all_h1-2)))
        print('Average betweenness centrality for subset normalized for throughflow h1')
        print(2*np.mean(bcs_h1[0:n_internal_h1])/(n_top_h1*n_bottom_h1))
        
        bcs_h2, nodecount_h2 = mm.betweenness_centrality_subset_n(pnx,False, space=1, axis=2)
        n_internal_h2, n_top_h2, n_bottom_h2 = nodecount_h2[0], nodecount_h2[1], nodecount_h2[2]
        n_all_h2 = np.sum(nodecount_h2)
        
        print('Average betweenness centrality for subset h2')
        print(2*np.mean(bcs_h2))
        print('Average normalized betweenness centrality for subset h2')
        print(2*np.mean(bcs_h2)/((n_all_h2-1)*(n_all_h2-2)))
        print('Average betweenness centrality for subset normalized for throughflow h2')
        print(2*np.mean(bcs_h2[0:n_internal_h2])/(n_top_h2*n_bottom_h2))
    
    plt.figure(num=571,figsize=(5,5))
    plt.clf()
    plt.subplot(1,1,1)
    tdbc_hist = plt.hist(2*bcs[0:n_internal]/(n_top*n_bottom), bins=20, density=False)
    plt.xlabel('Top-down betweenness centrality')
    plt.ylabel('Number')


if degree_switch:
    #Pore degree or coordination number
    
    nnodes = pnx_int.number_of_nodes()
    
    degrees = np.asarray([val for (node, val) in pnx_int.degree()])
    av_degree = np.sum(degrees)/nx.classes.function.number_of_nodes(pnx_int)
    
    print('Average pore coordination number')
    print(av_degree)
    
    degreehist = np.histogram(degrees, bins=np.arange(np.min(degrees),np.max(degrees)+2),density=False)
    
    
    fig1 = plt.figure(num=75, figsize=(5,5))
    plt.clf()
    ax2 = fig1.add_subplot(1,1,1)
    ax2.bar(degreehist[1][0:-1], degreehist[0], width = np.diff(degreehist[1]),align='edge')
    ax2.set_yscale('log')
    ax2.set_xlabel('Pore coordination number')
    ax2.set_ylabel('Number of pores')
    
    fig1 = plt.figure(num=85, figsize=(5,5))
    plt.clf()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.plot(pn_5['pore.coords'][pn_5['pore.internal'],0], degrees, 'b.')
    ax1.set_xlabel('Pore vertical coordinate')
    ax1.set_ylabel('Coordination number')


if clust_switch:

    #Clustering coefficient
    clust = np.asarray(list(nx.algorithms.cluster.clustering(pnx_int).values()))
    
    #Network average clustering coefficient
    clust_av = np.mean(clust)
    
    print('Network average clustering coefficient, single network')
    print(clust_av)
    
    clusthist = np.histogram(clust, bins=10,density=False)

    fig1 = plt.figure(num=76, figsize=(5,5))
    ax2 = fig1.add_subplot(1,1,1)
    ax2.bar(clusthist[1][0:-1], clusthist[0], width = np.diff(clusthist[1]),align='edge')
    ax2.set_yscale('log')
    ax2.set_xlabel('Pore clustering coefficient')
    ax2.set_ylabel('Number of pores')


    fig1 = plt.figure(num=6)
    fig1.set_size_inches(14,7)
    plt.clf()
    cnvscc_hist = []
    for i in range(2,17):
        try:
            ax1 = fig1.add_subplot(3,5,i-1)
            qwer = np.unique(clust[degrees==i])
            qwert = np.asarray(qwer[-1] + np.diff(qwer)[0])
            np.concatenate([qwer,np.array([qwert])])
            dummy1 =  ax1.hist(clust[degrees==i], bins=np.concatenate([qwer,np.array([qwert])]), density=False, align='left')
            cnvscc_hist.append(dummy1)
            ax1.set_ylim(0.9, 1.1*np.max(dummy1[0]))
            ax1.set_yscale('log')
            if i-1 == 1 or i-1 == 6 or i-1 == 11:
                ax1.set_ylabel('Pore clustering coefficient')
            if i-1 > 10:
                ax1.set_xlabel('Pore coordination number')
            ax1.set_title('CN = '+str(i))
        except Exception:
            pass
        
    fig1 = plt.figure(num=86, figsize=(5,5))
    plt.clf()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.plot(pn_5['pore.coords'][pn_5['pore.internal'],0], clust, 'b.')
    ax1.set_xlabel('Pore vertical coordinate')
    ax1.set_ylabel('Clustering coefficient')
    
    # Network transitivity (global clustering coefficient)
    transitiv = nx.algorithms.cluster.transitivity(pnx_int)
    
    print('Network transitivity')
    print(transitiv)


if centrality_switch:    
    
    print('Calculating closeness centrality')
    
    ccentr = []
    for i in range(0, pnx_int.number_of_nodes(), space):
        ccentr.append(nx.algorithms.centrality.closeness_centrality(pnx_int, u=i, distance='conduitlength'))
        if i%1000==0:
            print(i)
    
    #Closeness centrality
    ccentr = np.asarray(ccentr)
    ccentrhist = np.histogram(ccentr, bins=10,density=False)
    
    try:
        print(pn_5['pore.equivalent_diameter'][0])
    except:
        pn_5['pore.equivalent_diameter']=pn_5['pore.diameter']
        
    
    fig1 = plt.figure(num=777)
    fig1.set_size_inches(7,7)
    plt.clf()
    ax1 = fig1.add_subplot(2,2,1)
    ax1.hist(ccentr,bins=20)
    ax1.set_xlabel('Pore closeness centrality')
    ax1.set_ylabel('Number of pores')
    if space == 1: 
        ax2 = fig1.add_subplot(2,2,2)
        ax2.plot(pn_5['pore.coords'][pn_5['pore.internal'],0], ccentr, 'b.')
        ax2.set_xlabel('Vertical coordinate')
        ax2.set_ylabel('Closeness centrality')
        ax2.set_ylim([0,1.05*np.max(ccentr)])
        
        ax3 = fig1.add_subplot(2,2,3)
        ax3.scatter(pn_5['pore.equivalent_diameter'][pn_5['pore.internal']], ccentr, c='b')
        ax3.set_xlabel('Pore equivalent diameter')
        ax3.set_ylabel('Closeness centrality')
        ax3.set_ylim([0,1.05*np.max(ccentr)])
        ax3.set_xlim([-0.0005,0.005])
        
        ax4 = fig1.add_subplot(2,2,4)
        ax4.scatter(pn_5['pore.coords'][pn_5['pore.internal'],0],
                    pn_5['pore.equivalent_diameter'][pn_5['pore.internal']], c='b')
        ax4.set_xlabel('Vertical coordinate')
        ax4.set_ylabel('Pore equivalent diameter')
        ax4.set_xlim([-0.001,0.031])
        ax4.set_ylim([-0.0005,0.005])
    plt.tight_layout()
    
    print('Average closeness centrality')
    print(np.mean(ccentr))
    
    #Betweenness centrality
    bc_all = np.asarray(list(nx.algorithms.centrality.betweenness_centrality(pnx_int, k=np.int(pnx_int.number_of_nodes()/space), weight='conduitlength').values()))
    
    print('Average betweenness centrality for the whole network')
    print(np.mean(bc_all))        
    
    bcentrhist = np.histogram(bc_all, bins=20,density=False)
    
    fig78 = plt.figure(num=788)
    fig78.set_size_inches(7,7)
    plt.clf()
    ax1 = fig78.add_subplot(2,2,1)
    ax1.bar(bcentrhist[1][0:-1], bcentrhist[0], width = np.diff(bcentrhist[1]),align='edge')
    ax1.set_yscale('log')
    ax1.set_xlabel('Pore betweenness centrality')
    ax1.set_ylabel('Number of pores')
    if space == 1: 
        ax2 = fig78.add_subplot(2,2,2)
        ax2.plot(pn_5['pore.coords'][pn_5['pore.internal'],0], bc_all, 'b.')
        ax2.set_xlabel('Vertical coordinate')
        ax2.set_ylabel('Betweenness centrality')
        ax2.set_ylim([0,1.05*np.max(bc_all)])
        
        ax3 = fig78.add_subplot(2,2,3)
        ax3.scatter(pn_5['pore.equivalent_diameter'][pn_5['pore.internal']], bc_all, c='b')
        ax3.set_xlabel('Pore equivalent diameter')
        ax3.set_ylabel('Betweenness centrality')
        ax3.set_ylim([0,1.05*np.max(bc_all)])
        ax3.set_xlim([-0.0005,0.005])
    plt.tight_layout()
    
    fig1 = plt.figure(num=88, figsize=(11,4))
    plt.clf()
    ax1 = fig1.add_subplot(1,3,1)
    ax1.plot(pn_5['pore.coords'][pn_5['pore.internal'],0], bc_all, 'b.')
    ax2 = fig1.add_subplot(1,3,2)
    ax2.plot(pn_5['pore.coords'][pn_5['pore.internal'],1], bc_all, 'b.')
    ax3 = fig1.add_subplot(1,3,3)
    ax3.plot(pn_5['pore.coords'][pn_5['pore.internal'],2], bc_all, 'b.')    
    plt.suptitle('Betweenness centrality vs. directions')
