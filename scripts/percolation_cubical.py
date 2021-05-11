# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 15:17:06 2020

@author: pkiuru

Script for performing percolation simulations (drainage & imbibition) for a cubical-domain network
for volume fraction analysis

perc_diff_simulations.py with domain_switch = 'cube' & percolation_switch = 'True' must 
to be run in advance!

Required input:
    Output of perc_diff_simulations.py                                      

Output:
    Volume fraction of connected network as a function of external pressure
    (as NumPy arrays)

"""
# Save the percolation simulation results (pressure points, air-filled
# porosities) to a .npz file for drainage and imbibition
save_results = False
drainage_file = "../IO/results/volfrac_dra_0_5_2"
imbibition_file = "../IO/results/volfrac_imb_0_5_2" 

#Pressure scaling for output purposes
anglefactor = np.cos(np.radians(180)) / np.cos(np.radians(120))

#Applied pressures (Pa) in drainage percolation simulation with a contact angle of 180 degrees
points_vf = np.array([100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,
                       600.0, 700.0, 750.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 2000.0, 3000.0])
points_vf_orig = 0.001 * np.copy(points_vf) #Pa -> kPa

#Pressures with a contact angle of 120 degrees (value determined in the properties of Air)
points_vf /= anglefactor

#Applied pressures (Pa) in drainage percolation simulation with a contact angle of 180 degrees
points_imb = np.array([100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,
                       600.0, 700.0, 750.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 2000.0, 3000.0])
points_imb_orig = 0.001 * np.copy(points_imb) #Pa -> kPa

points_imb /= anglefactor

#Set the boundary pore volumes to zero
pn_5['pore.volume'][pn_5['pore.boundary']] = 0

'''
phys_air_5.add_model(propname='throat.entry_pressure',
         model=op.models.physics.capillary_pressure.washburn,
         diameter='throat.diameter')
phys_air_5['pore.entry_pressure'] = 0
'''

# Ordinary percolation simulation for air intrusion (drainage)
perc_5_vf = op.algorithms.OrdinaryPercolation(network=pn_5)
perc_5_vf.setup(phase=air_5)
perc_5_vf.set_inlets(pores=pn_5['pore.left'])
perc_5_vf.setup(pore_volume='pore.volume', throat_volume='throat.volume')
perc_5_vf.run(points = points_vf)

perc_5_vf_data = perc_5_vf.get_intrusion_data()

#Water phase
water_5 = op.phases.Water(network=pn_5)
water_5['pore.contact_angle'] = 60.0 #(180 degrees minus air contact angle)

phys_water_5 = op.physics.Standard(network=pn_5, phase=water_5, geometry=geom_5)

phys_water_5.add_model(propname='pore.entry_pressure',
         model=op.models.physics.capillary_pressure.washburn,
         diameter='pore.diameter')
phys_water_5['throat.entry_pressure'] = 100*np.min(phys_water_5['pore.entry_pressure'])
#a large negative value

# Ordinary percolation simulation for water intrusion (imbibition)
perc_5_w = op.algorithms.OrdinaryPercolation(network=pn_5)
perc_5_w.setup(phase=water_5, mode='site')
perc_5_w.set_inlets(pores=pn_5['pore.right'])
perc_5_w.setup(pore_volume='pore.volume', throat_volume='throat.volume')
perc_5_w.run(points = -np.flip(points_imb))

perc_5_w_data = perc_5_w.get_intrusion_data()


plt.figure(num=599, figsize=[11,5])
plt.clf()
plt.subplot(1,2,1)
plt.semilogy(1-np.asarray(perc_5_vf_data.Snwp), anglefactor * perc_5_vf_data.Pcap,
             'r-',label='Drainage')
plt.semilogy(np.asarray(perc_5_w_data.Snwp), -anglefactor * perc_5_w_data.Pcap,
             'k-',label='Imbibition')
plt.xlabel('Water saturation (water volume / total void space volume)')
plt.ylabel('Pressure (soil matric potential)')

plt.subplot(1,2,2)
plt.semilogy(np.asarray(poros_tot*perc_5_vf_data.Snwp), anglefactor * perc_5_vf_data.Pcap,
             'r-', label='Drainage')
plt.semilogy(poros_tot*(1-np.asarray(perc_5_w_data.Snwp)),
             -anglefactor * perc_5_w_data.Pcap, 'k-', label='Imbibition')
plt.xlabel('Air-filled porosity (air volume / total sample volume)')
plt.ylabel('Pressure (soil matric potential)')

# Air-filled porosity at each pressure in drainage (pressure = initial points_vf @ 180 degree contact angle)

afp = np.asarray(poros_tot*perc_5_vf_data.Snwp)

# Air-filled porosity at each pressure in imbibition

afp_imb = poros_tot*(1-np.asarray(perc_5_w_data.Snwp))

# Volume fraction of connected network as a function of external pressure
# (afp divided by the volume of total pore space)

vfcn = afp / (np.sum(pn_5w['pore.volume'][~pn_5w['pore.boundary']])/totvol)

vfcn_imb = afp_imb / (np.sum(pn_5w['pore.volume'][~pn_5w['pore.boundary']])/totvol)


if save_results:
    np.savez(drainage_file, prpoints=points_vf_orig, afporos=afp, 
             volfraction = vfcn)
    np.savez(imbibition_file, prpoints=points_imb_orig, 
             afporos=np.flip(afp_imb), volfraction = np.flip(vfcn_imb))
