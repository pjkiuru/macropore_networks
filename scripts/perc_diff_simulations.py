# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 12:40:29 2020

@author: pkiuru

Script for performing water retention (percolation) simulations (drainage
and/or imbibition) and relative diffusion simulations

Required input:
    OpenPNM network object from a file or from memory                                      

Output:
    Percolation result data (pressure, air-filled porosity), can be saved
    Relative diffusion results (printed in the console)
    NOTE: Complete version of OpenPNM network object, ready to be saved for
    later use


"""

import openpnm as op
import porespy as ps
import numpy as np
import matplotlib.pyplot as plt
import more_models as mm


# True: OpenPNM workspace (network object) from file; False: use the workspace
# in memory 
load_switch = True
# Path for the image array .npy file
ws_file = '../IO/pnms/40_45_1_600.pnm'


# Shape of the network domain ('cyl' or 'cube'):
# cyl: 800 x 1000 cylinder (height x diameter)
# cube: centered 600 x 600 x 600 cube
domain_switch = 'cube'

# NOTE: percolation and diffusion simulations in this script are performed for cylindrical
# for cylindrical domain networks only!
# Cubical domain percolation simulation is needed before the volume fraction analysis
# simulation in a following script!

if domain_switch == 'cyl':
    areafact = 0.785321
else:
    areafact = 1


# Perform percolation simulation (drainage and imbibition, i.e. air intrusion
# and water intrusion)
percolation_switch = True

# Save the percolation simulation results (pressure points, air-filled
# porosities) to a .npz file
save_percolation = False
#airpercresultfile = 'C:/Users/pkiuru/Pictures/metnet_ct_reco/perc_20_25_2'
#waterpercresultfile = 'C:/Users/pkiuru/Pictures/metnet_ct_reco/perc_water_0_5_1'

# Percolation result units:
# Pressure: Pa, calculated with contact angle 120 degrees
# air-filled porosity: volume of air-filled fraction of pore network divided by
# the total volume of the cuboid-shaped domain (the cylindrical shape is not taken
# into account here)

# Perform steady-state relative diffusion simulation
# NOTE: percolation simulation must be performed in advance!
diffusion_switch = True

# Save the OpenPNM network object to a .pnm file
save_network = False
#networkfile = 'C:/Users/pkiuru/Pictures/metnet_ct_reco/subs_20_25_2.pnm'

# -----------------------------------------------------------------------------

# Arbitrary boundary values for gas concentrations in the diffusion
# simulation; for the calculation of effective diffusivity
D_min = 1
D_max = 2


if load_switch:
    
    voxel_size = 50e-6
    
    if domain_switch == 'cube':
    
        z_top = 200
        z_bot = 800
        h_min = 200
        h_med = 500
        h_max = 800
        
    else:
        
        # Sample 1: z_top = 150, z_bot = 750 (z = 600)
        # Sample 2: 200...900 (z = 700)
        # Sample 7: 100...800 (z = 700)
        # Sample 9: 150...950 (z = 800)
        # Otherwise: 100...900 (z = 800)
        
        z_top = 100
        z_bot = 900
        h_min = 0
        h_med = 500
        h_max = 1000 
    
    subvol = (z_bot - z_top) * (h_med - h_min) * (h_med - h_min) * voxel_size**3
    totvol = (z_bot - z_top) * (h_max - h_min) * (h_max - h_min) * voxel_size**3
    
    ws = op.Workspace()
    ws.clear()
    ws.load_workspace(filename=ws_file, overwrite=True)
    
    proj1 = ws['subimage_1']
    proj2 = ws['subimage_2']
    proj3 = ws['subimage_3']
    proj4 = ws['subimage_4']
    
    pn_1 = proj1['net_01']
    pn_2 = proj2['net_01']
    pn_3 = proj3['net_01']
    pn_4 = proj4['net_01']
    
    geom_1 = proj1['geo_01']
    geom_2 = proj2['geo_01']
    geom_3 = proj3['geo_01']
    geom_4 = proj4['geo_01']
    
    proj1w = ws['whole_subnetwork_1']
    proj2w = ws['whole_subnetwork_2']
    proj3w = ws['whole_subnetwork_3']
    proj4w = ws['whole_subnetwork_4']
    
    
    if domain_switch == 'cube':

        proj5 = ws['image']
        pn_5 = proj5['net_01']
        geom_5 = proj5['geo_01']
        proj5w = ws['whole_network']
        pn_5w = proj5w['net_01']
        
        fig = plt.figure(num=22,figsize=(9,9))
        plt.clf()
        ax = fig.gca(projection='3d')
        ax.scatter(1000*pn_5['pore.coords'][:,0], 1000*pn_5['pore.coords'][:,1],
                   1000*pn_5['pore.coords'][:,2], c='#000080', s = 50e6*(pn_5['pore.equivalent_diameter'])**2,
                   cmap=plt.cm.jet)

    
if percolation_switch:
    
    # Porosities of the "cubical" subdomains - total volumes include the space
    # outside the cylinder
    
    poros_sub = np.array([np.sum(pn_1['pore.volume'][pn_1['pore.internal']]) / subvol,
                          np.sum(pn_2['pore.volume'][pn_2['pore.internal']]) / subvol,
                          np.sum(pn_3['pore.volume'][pn_3['pore.internal']]) / subvol,
                          np.sum(pn_4['pore.volume'][pn_4['pore.internal']]) / subvol])
    if domain_switch == 'cube':
        
        poros_tot = np.array([np.sum(pn_5['pore.volume'][pn_5['pore.internal']]) / totvol])
    
    # If no boundary pores exist, pore(s) nearest the boundary are "artificially"
    # labeled as boundary pore(s).
    
    if ~np.any(pn_1['pore.left']):
        
        print('No left boundary in net 1')
        leftmost = np.where(pn_1['pore.coords'][:,0] == np.min(pn_1['pore.coords'][:,0]))
        pn_1['pore.left'][leftmost] = True
        pn_1['pore.boundary'][leftmost] = True
        pn_1['pore.internal'][leftmost] = False
    
    if ~np.any(pn_2['pore.left']):
        
        print('No left boundary in net 2')
        leftmost = np.where(pn_2['pore.coords'][:,0] == np.min(pn_2['pore.coords'][:,0]))
        pn_2['pore.left'][leftmost] = True
        pn_2['pore.boundary'][leftmost] = True
        pn_2['pore.internal'][leftmost] = False
        
    if ~np.any(pn_3['pore.left']):
        
        print('No left boundary in net 3')
        leftmost = np.where(pn_3['pore.coords'][:,0] == np.min(pn_3['pore.coords'][:,0]))
        pn_3['pore.left'][leftmost] = True
        pn_3['pore.boundary'][leftmost] = True
        pn_3['pore.internal'][leftmost] = False
        
    if ~np.any(pn_4['pore.left']):
        
        print('No left boundary in net 4')
        leftmost = np.where(pn_4['pore.coords'][:,0] == np.min(pn_4['pore.coords'][:,0]))
        pn_4['pore.left'][leftmost] = True
        pn_4['pore.boundary'][leftmost] = True
        pn_4['pore.internal'][leftmost] = False
    
    if ~np.any(pn_1['pore.right']):
        
        print('No right boundary in net 1')
        rightmost = np.where(pn_1['pore.coords'][:,0] == np.max(pn_1['pore.coords'][:,0]))
        pn_1['pore.right'][rightmost] = True
        pn_1['pore.boundary'][rightmost] = True
        pn_1['pore.internal'][rightmost] = False
    
    if ~np.any(pn_2['pore.right']):
        
        print('No right boundary in net 2')
        rightmost = np.where(pn_2['pore.coords'][:,0] == np.max(pn_2['pore.coords'][:,0]))
        pn_2['pore.right'][rightmost] = True
        pn_2['pore.boundary'][rightmost] = True
        pn_2['pore.internal'][rightmost] = False
        
    if ~np.any(pn_3['pore.right']):
        
        print('No right boundary in net 3')
        rightmost = np.where(pn_3['pore.coords'][:,0] == np.max(pn_3['pore.coords'][:,0]))
        pn_3['pore.right'][rightmost] = True
        pn_3['pore.boundary'][rightmost] = True
        pn_3['pore.internal'][rightmost] = False
        
    if ~np.any(pn_4['pore.right']):
        
        print('No right boundary in net 4')
        rightmost = np.where(pn_4['pore.coords'][:,0] == np.max(pn_4['pore.coords'][:,0]))
        pn_4['pore.right'][rightmost] = True
        pn_4['pore.boundary'][rightmost] = True
        pn_4['pore.internal'][rightmost] = False
    
    
    if domain_switch == 'cube':
        
        if ~np.any(pn_5['pore.left']):
            
            print('No left boundary in net 5')
            leftmost = np.where(pn_5['pore.coords'][:,0] == np.min(pn_5['pore.coords'][:,0]))
            pn_5['pore.left'][leftmost] = True
            pn_5['pore.boundary'][leftmost] = True
            pn_5['pore.internal'][leftmost] = False
            
        if ~np.any(pn_5['pore.right']):
            
            print('No right boundary in net 5')
            rightmost = np.where(pn_5['pore.coords'][:,0] == np.max(pn_5['pore.coords'][:,0]))
            pn_5['pore.right'][rightmost] = True
            pn_5['pore.boundary'][rightmost] = True
            pn_5['pore.internal'][rightmost] = False

    
    # Set the properties of air phase and water phase for the subnetworks
    
    air_1 = op.phases.Air(network=pn_1)
    air_1['pore.contact_angle'] = 120.0
    phys_air_1 = op.physics.Standard(network=pn_1, phase=air_1, geometry=geom_1)
    
    water_1 = op.phases.Water(network=pn_1)
    water_1['pore.contact_angle'] = 60.0
    phys_water_1 = op.physics.Standard(network=pn_1, phase=water_1, geometry=geom_1)
    phys_water_1.add_model(propname='pore.entry_pressure',
                 model=op.models.physics.capillary_pressure.washburn,
                 diameter='pore.diameter')
        
    qw_1 = np.copy(pn_1['pore.volume'])
    
    
    air_2 = op.phases.Air(network=pn_2)
    air_2['pore.contact_angle'] = 120.0
    phys_air_2 = op.physics.Standard(network=pn_2, phase=air_2, geometry=geom_2)
    
    water_2 = op.phases.Water(network=pn_2)
    water_2['pore.contact_angle'] = 60.0
    phys_water_2 = op.physics.Standard(network=pn_2, phase=water_2, geometry=geom_2) 
    phys_water_2.add_model(propname='pore.entry_pressure',
                 model=op.models.physics.capillary_pressure.washburn,
                 diameter='pore.diameter')   
    
    qw_2 = np.copy(pn_2['pore.volume'])
    
    
    air_3 = op.phases.Air(network=pn_3)
    air_3['pore.contact_angle'] = 120.0
    phys_air_3 = op.physics.Standard(network=pn_3, phase=air_3, geometry=geom_3)
    
    water_3 = op.phases.Water(network=pn_3)
    water_3['pore.contact_angle'] = 60.0
    phys_water_3 = op.physics.Standard(network=pn_3, phase=water_3, geometry=geom_3)
    phys_water_3.add_model(propname='pore.entry_pressure',
                 model=op.models.physics.capillary_pressure.washburn,
                 diameter='pore.diameter')
    
    qw_3 = np.copy(pn_3['pore.volume'])
    
    
    air_4 = op.phases.Air(network=pn_4)
    air_4['pore.contact_angle'] = 120.0
    phys_air_4 = op.physics.Standard(network=pn_4, phase=air_4, geometry=geom_4)
    
    water_4 = op.phases.Water(network=pn_4)
    water_4['pore.contact_angle'] = 60.0
    phys_water_4 = op.physics.Standard(network=pn_4, phase=water_4, geometry=geom_4)
    phys_water_4.add_model(propname='pore.entry_pressure',
                 model=op.models.physics.capillary_pressure.washburn,
                 diameter='pore.diameter')

    qw_4 = np .copy(pn_4['pore.volume'])
    
    # Find total minimum and maximum entry pressures in the four subnetworks
    # and create 100 pressure points for simulation
    
    pres_min = np.min([np.min(phys_air_1['throat.entry_pressure']), np.min(phys_air_2['throat.entry_pressure']),
            np.min(phys_air_3['throat.entry_pressure']), np.min(phys_air_4['throat.entry_pressure'])])
    
    pres_max = np.max([np.max(phys_air_1['throat.entry_pressure']), np.max(phys_air_2['throat.entry_pressure']),
            np.max(phys_air_3['throat.entry_pressure']), np.max(phys_air_4['throat.entry_pressure'])])
    
    points = np.logspace(start=np.log10(max(1, 0.75*pres_min)),
                                 stop=np.log10(1.5*pres_max), num=100)
    
    
    pres_min_w = np.min([np.min(phys_water_1['pore.entry_pressure']), np.min(phys_water_2['pore.entry_pressure']),
            np.min(phys_water_3['pore.entry_pressure']), np.min(phys_water_4['pore.entry_pressure'])])
    
    pres_max_w = np.max([np.max(phys_water_1['pore.entry_pressure']), np.max(phys_water_2['pore.entry_pressure']),
            np.max(phys_water_3['pore.entry_pressure']), np.max(phys_water_4['pore.entry_pressure'])])
        
    points_w = -np.flip(np.logspace(start=np.log10(max(1, 0.75*(-pres_max_w))),
                                     stop=np.log10(-1.5*pres_min_w), num=100))
    
    # Run percolation simulation for air
    
    perc_1 = op.algorithms.OrdinaryPercolation(network=pn_1)
    perc_1.setup(phase=air_1)
    perc_1.set_inlets(pores=pn_1['pore.left'])
    perc_1.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_1.run(points = points)
    perc_1_data = perc_1.get_intrusion_data()
    
    
    perc_2 = op.algorithms.OrdinaryPercolation(network=pn_2)
    perc_2.setup(phase=air_2)
    perc_2.set_inlets(pores=pn_2['pore.left'])
    perc_2.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_2.run(points = points)
    perc_2_data = perc_2.get_intrusion_data()
    
    
    perc_3 = op.algorithms.OrdinaryPercolation(network=pn_3)
    perc_3.setup(phase=air_3)
    perc_3.set_inlets(pores=pn_3['pore.left'])
    perc_3.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_3.run(points = points)
    perc_3_data = perc_3.get_intrusion_data()
   
    
    perc_4 = op.algorithms.OrdinaryPercolation(network=pn_4)
    perc_4.setup(phase=air_4)
    perc_4.set_inlets(pores=pn_4['pore.left'])
    perc_4.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_4.run(points = points)
    perc_4_data = perc_4.get_intrusion_data()
    
    # Run percolation simulation for water
    
    perc_1_w = op.algorithms.OrdinaryPercolation(network=pn_1)
    perc_1_w.setup(phase=water_1, mode='site')
    perc_1_w.set_inlets(pores=pn_1['pore.right'])
    perc_1_w.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_1_w.run(points = points_w)
    perc_1_w_data = perc_1_w.get_intrusion_data()
    
    
    perc_2_w = op.algorithms.OrdinaryPercolation(network=pn_2)
    perc_2_w.setup(phase=water_2, mode='site')
    perc_2_w.set_inlets(pores=pn_2['pore.right'])
    perc_2_w.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_2_w.run(points = points_w)
    perc_2_w_data = perc_2_w.get_intrusion_data()
    
    
    perc_3_w = op.algorithms.OrdinaryPercolation(network=pn_3)
    perc_3_w.setup(phase=water_3, mode='site')
    perc_3_w.set_inlets(pores=pn_3['pore.right'])
    perc_3_w.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_3_w.run(points = points_w)
    perc_3_w_data = perc_3_w.get_intrusion_data()
   
    
    perc_4_w = op.algorithms.OrdinaryPercolation(network=pn_4)
    perc_4_w.setup(phase=water_4, mode='site')
    perc_4_w.set_inlets(pores=pn_4['pore.right'])
    perc_4_w.setup(pore_volume='pore.volume', throat_volume='throat.volume')
    perc_4_w.run(points = points_w)
    perc_4_w_data = perc_4_w.get_intrusion_data()
        
    
    sat_tot = (np.asarray(perc_1_data.Snwp) * np.sum(qw_1) \
            + np.asarray(perc_2_data.Snwp) * np.sum(qw_2) \
            + np.asarray(perc_3_data.Snwp) * np.sum(qw_3) \
            + np.asarray(perc_4_data.Snwp) * np.sum(qw_4)) \
            / (np.sum(qw_1) + np.sum(qw_2) +np.sum(qw_3) +np.sum(qw_4))
    
    sat_tot_w = (np.asarray(perc_1_w_data.Snwp) * np.sum(qw_1) \
            + np.asarray(perc_2_w_data.Snwp) * np.sum(qw_2) \
            + np.asarray(perc_3_w_data.Snwp) * np.sum(qw_3) \
            + np.asarray(perc_4_w_data.Snwp) * np.sum(qw_4)) \
            / (np.sum(qw_1) + np.sum(qw_2) +np.sum(qw_3) +np.sum(qw_4))
    
    # sat_tot is the degree of saturation (volume of pores filled with invading
    # fluid divided by volume of total network volume)
    
    if domain_switch == 'cube':
        
        air_5 = op.phases.Air(network=pn_5)
        air_5['pore.contact_angle'] = 120.0
        
        phys_air_5 = op.physics.Standard(network=pn_5, phase=air_5, geometry=geom_5)

        
        qw_5 = np.copy(pn_5['pore.volume'])
        
        perc_5 = op.algorithms.OrdinaryPercolation(network=pn_5)
        perc_5.setup(phase=air_5)
        perc_5.set_inlets(pores=pn_5['pore.left'])
        perc_5.setup(pore_volume='pore.volume', throat_volume='throat.volume')
        perc_5.run(points = points)
    
        perc_5_data = perc_5.get_intrusion_data()
        
        #------------
    
    #subplot 1: water saturation as a function of pressure
    
    #subplot 2: air-filled porosity as a function of pressure.
    # air-filled porosity = porosity * air saturation 
    # = (total network volume/total domain volume) * (air-filled volume/total network volume)
    
    plt.figure(num=942, figsize=[11,5])
    plt.clf()
    plt.subplot(1,2,1)
    plt.semilogy(1-np.asarray(perc_1_data.Snwp), perc_1_data.Pcap, 'b--', label='Subvolume 1',lw=0.5)
    plt.semilogy(1-np.asarray(perc_2_data.Snwp), perc_2_data.Pcap, 'r--', label='Subvolume 2',lw=0.5)
    plt.semilogy(1-np.asarray(perc_3_data.Snwp), perc_3_data.Pcap, 'g--', label='Subvolume 3',lw=0.5)
    plt.semilogy(1-np.asarray(perc_4_data.Snwp), perc_4_data.Pcap, 'k--', label='Subvolume 4',lw=0.5)
    
    plt.semilogy(np.asarray(perc_1_w_data.Snwp), -perc_1_w_data.Pcap, 'b-.', label='Imbibition subvolume 1',lw=0.5)
    plt.semilogy(np.asarray(perc_2_w_data.Snwp), -perc_2_w_data.Pcap, 'r-.', label='Imbibition subvolume 2',lw=0.5)
    plt.semilogy(np.asarray(perc_3_w_data.Snwp), -perc_3_w_data.Pcap, 'g-.', label='Imbibition subvolume 3',lw=0.5)
    plt.semilogy(np.asarray(perc_4_w_data.Snwp), -perc_4_w_data.Pcap, 'k-.', label='Imbibition subvolume 4',lw=0.5)
    
    plt.semilogy(1-sat_tot, perc_4_data.Pcap, 'k-', label='Combined volume')
    plt.semilogy(sat_tot_w, -perc_4_w_data.Pcap, 'y-', label='Imbibition combined volume')
    #if domain_switch == 'cube':
        #plt.semilogy(1-np.asarray(perc_5_data.Snwp), perc_5_data.Pcap, 'r-', label='Single volume')
        #plt.semilogy(1-inj_data.S_tot, inj_data.Pcap, 'g-', label='MIP')
        #plt.semilogy(np.asarray(perc_5_w_data.Snwp), -perc_5_w_data.Pcap, 'c-', label='Single volume imbibition')
        #plt.semilogy(wtd_data.S_tot, -wtd_data.Pcap, 'm-', label='MIP imbibition')
    plt.xlabel('Water saturation (water volume / total void space volume)')
    plt.ylabel('Pressure (soil matric potential)')
    plt.legend()
    plt.subplot(1,2,2)
    plt.semilogy(poros_sub[0]*np.asarray(perc_1_data.Snwp)/areafact, perc_1_data.Pcap, 'b--', label='Subvolume 1',lw=0.5)
    plt.semilogy(poros_sub[1]*np.asarray(perc_2_data.Snwp)/areafact, perc_2_data.Pcap, 'r--', label='Subvolume 2',lw=0.5)
    plt.semilogy(poros_sub[2]*np.asarray(perc_3_data.Snwp)/areafact, perc_3_data.Pcap, 'g--', label='Subvolume 3',lw=0.5)
    plt.semilogy(poros_sub[3]*np.asarray(perc_4_data.Snwp)/areafact, perc_4_data.Pcap, 'k--', label='Subvolume 4',lw=0.5)
    
    plt.semilogy(poros_sub[0]*(1-np.asarray(perc_1_w_data.Snwp))/areafact, -perc_1_w_data.Pcap, 'b-.', label='Imbibition subvolume 1',lw=0.5)
    plt.semilogy(poros_sub[1]*(1-np.asarray(perc_2_w_data.Snwp))/areafact, -perc_2_w_data.Pcap, 'r-.', label='Imbibition subvolume 2',lw=0.5)
    plt.semilogy(poros_sub[2]*(1-np.asarray(perc_3_w_data.Snwp))/areafact, -perc_3_w_data.Pcap, 'g-.', label='Imbibition subvolume 3',lw=0.5)
    plt.semilogy(poros_sub[3]*(1-np.asarray(perc_4_w_data.Snwp))/areafact, -perc_4_w_data.Pcap, 'k-.', label='Imbibition subvolume 4',lw=0.5)
    
    plt.semilogy(np.mean(poros_sub)*sat_tot/areafact, perc_4_data.Pcap, 'k-', label='Combined volume')
    plt.semilogy(np.mean(poros_sub)*(1-sat_tot_w)/areafact, -perc_4_w_data.Pcap, 'y-', label='Combined volume')
    #if domain_switch == 'cube':
        #plt.semilogy(np.asarray(poros_tot*perc_5_data.Snwp), perc_5_data.Pcap, 'r-', label='Single volume')
        #plt.semilogy(poros_tot*(1-np.asarray(perc_5_w_data.Snwp)), -perc_5_w_data.Pcap, 'c-', label='Single volume imbibition')
    plt.xlabel('Air-filled porosity (air volume / total sample volume)')
    plt.ylabel('Pressure (soil matric potential)')
    plt.legend()
        
    
    if save_percolation:
        
        np.savez(airpercresultfile, press=perc_4_data.Pcap, afporos=np.mean(poros_sub)*sat_tot)
        np.savez(waterpercresultfile, press_w=-perc_4_data.Pcap, afporos_w=np.mean(poros_sub)*(1-sat_tot_w))
    
    if diffusion_switch:
        
        #Relative diffusion
        
        
        #Information on pore occupancy is required for the determination of
        #multiphase conduit conductance
        #This initial occupancy is random and will be updated and corrected
        #as the simulation starts
        air_1.update(perc_1.results(Pc=1000))
        
        pn_1.add_model(propname='throat.length_piecewise',
                 model=op.models.geometry.throat_length.piecewise,
                 throat_endpoints='throat.endpoints',
                 throat_centroid='throat.centroid')
            
        pn_1.add_model(propname='throat.conduit_lengths',
                 model=op.models.geometry.throat_length.conduit_lengths,
                 throat_length='throat.length_piecewise')
        
        phys_air_1.add_model(propname='throat.diffusive_conductance',
                             model=op.models.physics.diffusive_conductance.ordinary_diffusion,
                             regen_mode='normal')
        
        phys_air_1.add_model(model=op.models.physics.multiphase.conduit_conductance,
                       propname='throat.conduit_diffusive_conductance',
                       throat_conductance='throat.diffusive_conductance')    
        
        diff_air_1 = []
        sat_1 = []
        tot_vol_1 = np.sum(pn_1["pore.volume"]) + np.sum(pn_1["throat.volume"])
        rate_1 = []
        
        for Pc in points:#np.unique(perc_1['pore.invasion_pressure']):
            air_1.update(perc_1.results(Pc=Pc))
            phys_air_1.regenerate_models()
            this_sat = 0
            this_sat += np.sum(pn_1["pore.volume"][air_1["pore.occupancy"] == 1])
            sat_1.append(this_sat)
            BC1_pores = pn_1['pore.left']
            BC2_pores = pn_1['pore.right']
            FD_1 = op.algorithms.FickianDiffusion(network=pn_1)
            FD_1.setup(phase=air_1, conductance='throat.conduit_diffusive_conductance')
            FD_1.set_value_BC(values=D_max, pores=BC1_pores)
            FD_1.set_value_BC(values=D_min, pores=BC2_pores)
            FD_1.run()
            rate_1.append(FD_1.rate(pores=BC1_pores))
            eff_diff = FD_1.calc_effective_diffusivity(domain_area=areafact*((h_med-h_min)*voxel_size)**2, domain_length= (z_bot-z_top)*voxel_size)
            if np.size(eff_diff)==0:
                eff_diff = 0
            diff_air_1.append(eff_diff)
            pn_1.project.purge_object(FD_1)
        
        sat_1 = np.asarray(sat_1)
        sats_1 = np.copy(sat_1)
        sat_1 /= tot_vol_1
        
        rel_diff_air_1 = np.asarray(diff_air_1)
        if rel_diff_air_1[-1] > 0:
            rel_diff_air_1 /= rel_diff_air_1[-1]
        
        
        air_2.update(perc_2.results(Pc=1000))
        
        pn_2.add_model(propname='throat.length_piecewise',
                 model=op.models.geometry.throat_length.piecewise,
                 throat_endpoints='throat.endpoints',
                 throat_centroid='throat.centroid')
            
        pn_2.add_model(propname='throat.conduit_lengths',
                 model=op.models.geometry.throat_length.conduit_lengths,
                 throat_length='throat.length_piecewise')
        
        phys_air_2.add_model(propname='throat.diffusive_conductance',
                         model=op.models.physics.diffusive_conductance.ordinary_diffusion,
                         regen_mode='normal')
        
        phys_air_2.add_model(model=op.models.physics.multiphase.conduit_conductance,
                       propname='throat.conduit_diffusive_conductance',
                       throat_conductance='throat.diffusive_conductance')    
        
        diff_air_2 = []
        sat_2 = []
        tot_vol_2 = np.sum(pn_2["pore.volume"]) + np.sum(pn_2["throat.volume"])
        rate_2 = []
        
        for Pc in points:#np.unique(perc_2['pore.invasion_pressure']):
            air_2.update(perc_2.results(Pc=Pc))
            phys_air_2.regenerate_models()
            this_sat = 0
            this_sat += np.sum(pn_2["pore.volume"][air_2["pore.occupancy"] == 1])
            sat_2.append(this_sat)
            BC1_pores = pn_2['pore.left']
            BC2_pores = pn_2['pore.right']
            FD_2 = op.algorithms.FickianDiffusion(network=pn_2)
            FD_2.setup(phase=air_2, conductance='throat.conduit_diffusive_conductance')
            FD_2.set_value_BC(values=D_max, pores=BC1_pores)
            FD_2.set_value_BC(values=D_min, pores=BC2_pores)
            FD_2.run()
            rate_2.append(FD_2.rate(pores=BC1_pores))
            eff_diff = FD_2.calc_effective_diffusivity(domain_area=areafact*((h_med-h_min)*voxel_size)**2, domain_length= (z_bot-z_top)*voxel_size)
            if np.size(eff_diff)==0:
                eff_diff = 0
            diff_air_2.append(eff_diff)
            pn_2.project.purge_object(FD_2)
        
        sat_2 = np.asarray(sat_2)
        sats_2 = np.copy(sat_2)
        sat_2 /= tot_vol_2
        
        rel_diff_air_2 = np.asarray(diff_air_2)
        if rel_diff_air_2[-1] > 0:
            rel_diff_air_2 /= rel_diff_air_2[-1]
        
        
        air_3.update(perc_3.results(Pc=1000))
        
        pn_3.add_model(propname='throat.length_piecewise',
                 model=op.models.geometry.throat_length.piecewise,
                 throat_endpoints='throat.endpoints',
                 throat_centroid='throat.centroid')
            
        pn_3.add_model(propname='throat.conduit_lengths',
                 model=op.models.geometry.throat_length.conduit_lengths,
                 throat_length='throat.length_piecewise')
        
        phys_air_3.add_model(propname='throat.diffusive_conductance',
                                 model=op.models.physics.diffusive_conductance.ordinary_diffusion,
                                 regen_mode='normal')
    
        phys_air_3.add_model(model=op.models.physics.multiphase.conduit_conductance,
                       propname='throat.conduit_diffusive_conductance',
                       throat_conductance='throat.diffusive_conductance')    
        
        diff_air_3 = []
        sat_3 = []
        tot_vol_3 = np.sum(pn_3["pore.volume"]) + np.sum(pn_3["throat.volume"])
        rate_3 = []
        print('3 starts')
        for Pc in points:#np.unique(perc_3['pore.invasion_pressure']):
            air_3.update(perc_3.results(Pc=Pc))
            phys_air_3.regenerate_models()
            this_sat = 0
            this_sat += np.sum(pn_3["pore.volume"][air_3["pore.occupancy"] == 1])
            sat_3.append(this_sat)
            BC1_pores = pn_3['pore.left']
            BC2_pores = pn_3['pore.right']
            FD_3 = op.algorithms.FickianDiffusion(network=pn_3)
            FD_3.setup(phase=air_3, conductance='throat.conduit_diffusive_conductance')
            FD_3.set_value_BC(values=D_max, pores=BC1_pores)
            FD_3.set_value_BC(values=D_min, pores=BC2_pores)
            FD_3.run()
            rate_3.append(FD_3.rate(pores=BC1_pores))
            eff_diff = FD_3.calc_effective_diffusivity(domain_area=areafact*((h_med-h_min)*voxel_size)**2, domain_length= (z_bot-z_top)*voxel_size)
            if np.size(eff_diff)==0:
                eff_diff = 0
            diff_air_3.append(eff_diff)
            pn_3.project.purge_object(FD_3)
        
        sat_3 = np.asarray(sat_3)
        sats_3 = np.copy(sat_3)
        sat_3 /= tot_vol_3
        print('3 ends')
        rel_diff_air_3 = np.asarray(diff_air_3)
        if rel_diff_air_3[-1] > 0:
            rel_diff_air_3 /= rel_diff_air_3[-1]
        
        
        air_4.update(perc_4.results(Pc=1000))
        
        pn_4.add_model(propname='throat.length_piecewise',
                 model=op.models.geometry.throat_length.piecewise,
                 throat_endpoints='throat.endpoints',
                 throat_centroid='throat.centroid')
            
        pn_4.add_model(propname='throat.conduit_lengths',
                 model=op.models.geometry.throat_length.conduit_lengths,
                 throat_length='throat.length_piecewise')
        
        phys_air_4.add_model(propname='throat.diffusive_conductance',
                             model=op.models.physics.diffusive_conductance.ordinary_diffusion,
                             regen_mode='normal')
     
        phys_air_4.add_model(model=op.models.physics.multiphase.conduit_conductance,
                       propname='throat.conduit_diffusive_conductance',
                       throat_conductance='throat.diffusive_conductance')    
        
        diff_air_4 = []
        sat_4 = []
        tot_vol_4 = np.sum(pn_4["pore.volume"]) + np.sum(pn_4["throat.volume"])
        rate_4 = []
        
        for Pc in points:#np.unique(perc_4['pore.invasion_pressure']):
            air_4.update(perc_4.results(Pc=Pc))
            phys_air_4.regenerate_models()
            this_sat = 0
            this_sat += np.sum(pn_4["pore.volume"][air_4["pore.occupancy"] == 1])
            sat_4.append(this_sat)
            BC1_pores = pn_4['pore.left']
            BC2_pores = pn_4['pore.right']
            FD_4 = op.algorithms.FickianDiffusion(network=pn_4)
            FD_4.setup(phase=air_4, conductance='throat.conduit_diffusive_conductance')
            FD_4.set_value_BC(values=D_max, pores=BC1_pores)
            FD_4.set_value_BC(values=D_min, pores=BC2_pores)
            FD_4.run()
            rate_4.append(FD_4.rate(pores=BC1_pores))
            eff_diff = FD_4.calc_effective_diffusivity(domain_area=areafact*((h_med-h_min)*voxel_size)**2, domain_length= (z_bot-z_top)*voxel_size)
            if np.size(eff_diff)==0:
                eff_diff = 0
            diff_air_4.append(eff_diff)
            pn_4.project.purge_object(FD_4)
        
        sat_4 = np.asarray(sat_4)
        sats_4 = np.copy(sat_4)
        sat_4 /= tot_vol_4
        
        rel_diff_air_4 = np.asarray(diff_air_4)
        if rel_diff_air_4[-1] > 0:
            rel_diff_air_4 /= rel_diff_air_4[-1]
        
        # Combined flow rate through all the subnetworks is obtained by summing up
        # the individual flow rates
        
        rate_comb = np.asarray(rate_1) + np.asarray(rate_2) + np.asarray(rate_3) + \
                         np.asarray(rate_4)
        
        # Effective diffusion is calculated by inverting the bulk flow rate equation
        
        diff_air_comb = (np.asarray(rate_1) + np.asarray(rate_2) + np.asarray(rate_3) + 
                         np.asarray(rate_4)) * (z_bot-z_top)*voxel_size \
                         / ((D_max - D_min) * areafact* ((h_max-h_min)*voxel_size)**2)
        
        # Diffusion rate at a certain pressure relative to maximum diffusion rate
        
        rel_diff_air_comb = diff_air_comb / diff_air_comb[-1]
        
        #Air saturation of the combined network
        sat_diff_tot = (sats_1 + sats_2 + sats_3 + sats_4) / (tot_vol_1 + tot_vol_2 + tot_vol_3 + tot_vol_4)
        
        if domain_switch == 'cube':
            
            air_5.update(perc_5.results(Pc=1000))
            
            pn_5.add_model(propname='throat.length_piecewise',
                 model=op.models.geometry.throat_length.piecewise,
                 throat_endpoints='throat.endpoints',
                 throat_centroid='throat.centroid')
            
            pn_5.add_model(propname='throat.conduit_lengths',
                 model=op.models.geometry.throat_length.conduit_lengths,
                 throat_length='throat.length_piecewise')
            
            phys_air_5.add_model(propname='throat.diffusive_conductance',
                             model=op.models.physics.diffusive_conductance.ordinary_diffusion,
                             regen_mode='normal')
        
            phys_air_5.add_model(model=op.models.physics.multiphase.conduit_conductance,
                           propname='throat.conduit_diffusive_conductance',
                           throat_conductance='throat.diffusive_conductance')    
            
            diff_air_5 = []
            sat_5 = []
            tot_vol_5 = np.sum(pn_5["pore.volume"]) + np.sum(pn_5["throat.volume"])
            rate_5 = []
            
            for Pc in points:#np.unique(perc_5['pore.invasion_pressure']):
                air_5.update(perc_5.results(Pc=Pc))
                phys_air_5.regenerate_models()
                this_sat = 0
                this_sat += np.sum(pn_5["pore.volume"][air_5["pore.occupancy"] == 1])
                sat_5.append(this_sat)
                BC1_pores = pn_5['pore.left']
                BC2_pores = pn_5['pore.right']
                FD_5 = op.algorithms.FickianDiffusion(network=pn_5)
                FD_5.setup(phase=air_5, conductance='throat.conduit_diffusive_conductance')
                FD_5.set_value_BC(values=D_max, pores=BC1_pores)
                FD_5.set_value_BC(values=D_min, pores=BC2_pores)
                FD_5.run()
                rate_5.append(FD_5.rate(pores=BC1_pores))
                eff_diff = FD_5.calc_effective_diffusivity(domain_area=areafact*((h_max-h_min)*voxel_size)**2, domain_length= (z_bot-z_top)*voxel_size)
                diff_air_5.append(eff_diff)
                pn_5.project.purge_object(FD_5)
            
            sat_5 = np.asarray(sat_5)
            sat_5 /= tot_vol_5
            
            rel_diff_air_5 = np.asarray(diff_air_5)
            rel_diff_air_5 /= rel_diff_air_5[-1]
        
        plt.figure(188)
        plt.clf()
        plt.plot(sat_1, np.asarray(rel_diff_air_1), '^-r', label='1')
        plt.plot(sat_2, np.asarray(rel_diff_air_2), '^-g', label='2')
        plt.plot(sat_3, np.asarray(rel_diff_air_3), '^-b', label='3')
        plt.plot(sat_4, np.asarray(rel_diff_air_4), '^-k', label='4')
        plt.plot(sat_diff_tot, rel_diff_air_comb, '^-c', label='1...4')
        if domain_switch == 'cube':
            plt.plot(sat_5, np.asarray(rel_diff_air_5), '^-m', label='5')
        plt.legend()
        plt.xlabel('Air saturation (Air volume / total void space volume)')
        plt.ylabel('Relative diffusivity (relative to totally dry pore conditions)')
        
        plt.figure(189)
        plt.clf()
        plt.plot(points, np.asarray(rate_1), '^-r', label='Subvolume 1')
        plt.plot(points, np.asarray(rate_2), '^-g', label='Subvolume 2')
        plt.plot(points, np.asarray(rate_3), '^-b', label='Subvolume 3')
        plt.plot(points, np.asarray(rate_4), '^-k', label='Subvolume 4')
        plt.plot(points, rate_comb, '^-c', label='Combined volume')
        if domain_switch == 'cube':
            plt.plot(points, np.asarray(rate_5), '^-m', label='Single volume')
        plt.legend()
        plt.xlabel('Pressure)')
        plt.ylabel('Flow')
        
        plt.figure(190,figsize=(13,5))
        plt.clf()
        plt.subplot(1,2,1)
        plt.plot(sat_1, np.asarray(rate_1), '-r', label='Subvolume 1',lw=0.5)
        plt.plot(sat_2, np.asarray(rate_2), '-g', label='Subvolume 2',lw=0.5)
        plt.plot(sat_3, np.asarray(rate_3), '-b', label='Subvolume 3',lw=0.5)
        plt.plot(sat_4, np.asarray(rate_4), '-k', label='Subvolume 4',lw=0.5)
        plt.plot(sat_diff_tot, rate_comb, '^-c', label='Combined volume')
        if domain_switch == 'cube':
            plt.plot(sat_5, np.asarray(rate_5), '^-m', label='Single volume')
        plt.legend()
        plt.xlabel('Air saturation (air volume / total void space volume)')
        plt.ylabel('Flow rate [units per second]')
        
        plt.subplot(1,2,2)
        plt.plot(poros_sub[0]*sat_1, np.asarray(rate_1), '-r', label='Subvolume 1',lw=0.5)
        plt.plot(poros_sub[1]*sat_2, np.asarray(rate_2), '-g', label='Subvolume 2',lw=0.5)
        plt.plot(poros_sub[2]*sat_3, np.asarray(rate_3), '-b', label='Subvolume 3',lw=0.5)
        plt.plot(poros_sub[3]*sat_4, np.asarray(rate_4), '-k', label='Subvolume 4',lw=0.5)
        plt.plot(np.mean(poros_sub)*sat_diff_tot, rate_comb, '^-c', label='Combined volume')
        if domain_switch == 'cube':
            plt.plot(sat_5, np.asarray(rate_5), '^-m', label='Single volume')
        plt.legend()
        plt.xlabel('Air-filled porosity (air volume / total sample volume)')
        plt.ylabel('Flow rate [units per second]')
        
        
        plt.figure(191,figsize=(13,5))
        plt.clf()
        plt.subplot(1,2,1)
        plt.plot(sat_1, 1e4*np.asarray(diff_air_1), '-r', label='Subvolume 1',lw=0.5)
        plt.plot(sat_2, 1e4*np.asarray(diff_air_2), '-g', label='Subvolume 2',lw=0.5)
        plt.plot(sat_3, 1e4*np.asarray(diff_air_3), '-b', label='Subvolume 3',lw=0.5)
        plt.plot(sat_4, 1e4*np.asarray(diff_air_4), '-k', label='Subvolume 4',lw=0.5)
        plt.plot(sat_diff_tot, 1e4*diff_air_comb, '^-c', label='Combined volume')
        if domain_switch == 'cube':
            plt.plot(sat_5, 1e4*np.asarray(diff_air_5), '^-m', label='Single volume')
        plt.legend()
        plt.xlabel('Air saturation (air volume / total void space volume)')
        plt.ylabel('Diffusivity D [cm2/s]')
        plt.subplot(1,2,2)
        plt.plot(poros_sub[0]*sat_1, 1e4*np.asarray(diff_air_1), '-r', label='Subvolume 1',lw=0.5)
        plt.plot(poros_sub[1]*sat_2, 1e4*np.asarray(diff_air_2), '-g', label='Subvolume 2',lw=0.5)
        plt.plot(poros_sub[2]*sat_3, 1e4*np.asarray(diff_air_3), '-b', label='Subvolume 3',lw=0.5)
        plt.plot(poros_sub[3]*sat_4, 1e4*np.asarray(diff_air_4), '-k', label='Subvolume 4',lw=0.5)
        plt.plot(np.mean(poros_sub)*sat_diff_tot, 1e4*diff_air_comb, '^-c', label='Combined volume')
        if domain_switch == 'cube':
            plt.plot(poros_tot*sat_5, 1e4*np.asarray(diff_air_5), '^-m', label='Single volume')
        plt.legend()
        plt.xlabel('Air-filled porosity (air volume / total sample volume)')
        plt.ylabel('Diffusivity D [cm2/s]')
        
        print('Top boundary pores')
        print( np.sum(pn_1['pore.left']), np.sum(pn_2['pore.left']),
               np.sum(pn_3['pore.left']), np.sum(pn_4['pore.left']))
        
        print('Bottom boundary pores')
        print( np.sum(pn_1['pore.right']), np.sum(pn_2['pore.right']),
               np.sum(pn_3['pore.right']), np.sum(pn_4['pore.right']))
        
        print('Maximum diffusivities')
        print( diff_air_1[-1], diff_air_2[-1],
               diff_air_3[-1], diff_air_4[-1])
        
        print('Maximum relative diffusivities')
        print( diff_air_1[-1]/air_1['pore.diffusivity'][0], diff_air_2[-1]/air_1['pore.diffusivity'][0],
               diff_air_3[-1]/air_1['pore.diffusivity'][0], diff_air_4[-1]/air_1['pore.diffusivity'][0])
        
        print('Pore counts')
        print( np.sum(~pn_1['pore.boundary']),np.sum(~pn_2['pore.boundary']),
               np.sum(~pn_3['pore.boundary']),np.sum(~pn_4['pore.boundary']))
        
        print('Combined pore count')
        print( np.sum(~pn_1['pore.boundary']) + np.sum(~pn_2['pore.boundary']) +
               np.sum(~pn_3['pore.boundary']) + np.sum(~pn_4['pore.boundary']))
        
        print('Throat counts')
        print( np.sum(pn_1['throat.internal']),np.sum(pn_2['throat.internal']),
               np.sum(pn_3['throat.internal']),np.sum(pn_4['throat.internal']))
        
        print('Combined throat count')
        print( np.sum(pn_1['throat.internal']) + np.sum(pn_2['throat.internal']) +
               np.sum(pn_3['throat.internal']) + np.sum(pn_4['throat.internal']))
        
        if domain_switch == 'cube':
            
            print('Single pore count')
            print(np.sum(pn_5['pore.internal']))
            
            print('Single throat count')
            print(np.sum(pn_5['throat.internal']))
            
            print('Top boundary pores, single')
            print( np.sum(pn_5['pore.left']))
            
            print('Bottom boundary pores, single')
            print( np.sum(pn_5['pore.right']))
                
        print('Total pore volumes')
        print( np.sum(pn_1['pore.volume'][~pn_1['pore.boundary']])/voxel_size**3,
               np.sum(pn_2['pore.volume'][~pn_2['pore.boundary']])/voxel_size**3,
               np.sum(pn_3['pore.volume'][~pn_3['pore.boundary']])/voxel_size**3,
               np.sum(pn_4['pore.volume'][~pn_4['pore.boundary']])/voxel_size**3)
        
        print('Combined pore volume')
        print( np.sum(np.sum(pn_1['pore.volume'][~pn_1['pore.boundary']]) +
               np.sum(pn_2['pore.volume'][~pn_2['pore.boundary']]) +
               np.sum(pn_3['pore.volume'][~pn_3['pore.boundary']]) +
               np.sum(pn_4['pore.volume'][~pn_4['pore.boundary']]))/voxel_size**3)
    
        
        print('Network porosities')
        print(poros_sub/areafact)
        
        print('Combined porosity')
        print( np.sum(np.sum(pn_1['pore.volume'][~pn_1['pore.boundary']]) +
               np.sum(pn_2['pore.volume'][~pn_2['pore.boundary']]) +
               np.sum(pn_3['pore.volume'][~pn_3['pore.boundary']]) +
               np.sum(pn_4['pore.volume'][~pn_4['pore.boundary']]))/(totvol*areafact))
        
        # Effective diffusivity of the totally air-filled combination of the four subnetworks
        
        total_diffusivity = (np.asarray(rate_1[-1])+np.asarray(rate_2[-1])+np.asarray(rate_3[-1])+np.asarray(rate_4[-1])) \
                        * ((z_bot-z_top)*voxel_size)/((D_max - D_min) * areafact * ((h_max-h_min)*voxel_size)**2)
        
        print('Combined diffusivity of four subnetworks')
        print(total_diffusivity)
        
        print('Total flow through four subnetworks')
        print(rate_comb[-1])
        
        print('Subnetwork flows')
        print(rate_1[-1], rate_2[-1], rate_3[-1], rate_4[-1])
        
        
        if domain_switch == 'cube':
        
            print('Diffusivity of total network')
            print(diff_air_5[-1])
            
            print('Total network pore volume')
            print(np.sum(pn_5['pore.volume'][~pn_5['pore.boundary']]))
            
            print('Total network porosity - cubical domain')
            print(np.sum(pn_5['pore.volume'][~pn_5['pore.boundary']])/totvol)
        
    
if save_network:
    ws.save_workspace(filename=networkfile)
    
    