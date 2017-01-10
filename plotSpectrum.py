"""
 plotSpectrum.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the output of the XMM THELSim simulation 
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotSpectrum.py N_file source N_in type E_bin E_bin_part E_min E_max
 ---------------------------------------------------------------------------------
 Parameters:
 - N_file = number of files
 - source = background source [0 = CXB]
 - N_in = number of simulated particles
 - type = 0[Photons], 1[Electrons], 2[Positrons], 3[Detector]
 - E_bin = energy bin spectrum
 - E_min = energy min
 - E_max = energy max
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2015/02/12: creation date
"""

from astropy.io import fits
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
from collections import Counter

# Import the input parameters
arg_list = sys.argv
N_file = int(arg_list[1])
source = int(arg_list[2])
N_in = int(arg_list[3])
type = int(arg_list[4])
E_bin = float(arg_list[5])
E_bin_part = float(arg_list[6])
E_min = float(arg_list[7])
E_max = float(arg_list[8])

# set-up
if source == 0: source_name = 'CXB'
if type == 0: type_name = 'Photon'
if type == 1: type_name = 'Electron'
if type == 2: type_name = 'Positron'

# Normalization
flux_part = 197.2  #part/cm2/s
det_side = 7.67  # cm
det_area = np.pi*((det_side/2.)**2.)  #cm2
flux_part_sr = flux_part/(4.*np.pi)  # part/cm2 s sr

R_big = 900.  # Sphere radius - cm
q = 0.067  # Half angle of aperture - rad

A_sphere = 4.*np.pi*(R_big**2.)   
omega = np.pi*((np.sin(q))**2.)  # sr
F_int = flux_part_sr*A_sphere*omega  #part/s
Time = (N_in*1000000.)/F_int   #s 

vecParticleEnergy = []
vecParticleProcName = []
vecDetEnergy = []

countPart = 0
countDet = 0

N_each = N_in/N_file
for jfits in xrange(N_file):

	path_output = './'+source_name+'/'+str(N_each*(jfits+1))+'mil/'
	hdulist_part = fits.open(path_output+'/'+source_name+'_'+str(N_each*(jfits+1))+'mil_XMMSim_'+type_name+'s.fits.gz')
	
	tbdata_part = hdulist_part[1].data
	
	evt_id_part = tbdata_part.field('EVENT_ID')
	energy_part = tbdata_part.field('ENERGY')
	e_dep_part = tbdata_part.field('E_DEP')
	pos_x_part = tbdata_part.field('POS_X')
	pos_y_part = tbdata_part.field('POS_Y')
	proc_id_part = tbdata_part.field('PROC_ID')

	hdulist_det = fits.open(path_output+'/'+source_name+'_'+str(N_each*(jfits+1))+'mil_XMMSim_Detector.fits.gz')	
	tbdata_det = hdulist_det[1].data
	
	evt_id_det = tbdata_det.field('EVENT_ID')
	e_dep_det = tbdata_det.field('E_DEP')

	where_esel = np.where((e_dep_det >= E_min) & (e_dep_det <= E_max))  
	evt_id_det_sel = evt_id_det[where_esel[0]]
	e_dep_det_sel = e_dep_det[where_esel[0]]
	
	for jev in xrange(len(evt_id_det_sel)):
		vecDetEnergy.append(e_dep_det_sel[jev])
		countDet += 1
		where_evsel = np.where(evt_id_part == evt_id_det_sel[jev])  		
		if where_evsel[0].size:
			countPart += 1
			temp_enpart = energy_part[where_evsel[0]]
			temp_proc_id = proc_id_part[where_evsel[0]]
			for jpart in xrange(len(temp_enpart)):
				vecParticleEnergy.append(temp_enpart[jpart])
				if (temp_proc_id[jpart] == 1): vecParticleProcName.append('Primary')
				if (temp_proc_id[jpart] == 2): vecParticleProcName.append('Photoelectric')
				if (temp_proc_id[jpart] == 3): vecParticleProcName.append('Compton')
				if (temp_proc_id[jpart] == 4): vecParticleProcName.append('Bremsstrahlung')
				if (temp_proc_id[jpart] == 5): vecParticleProcName.append('eIonization')
				if (temp_proc_id[jpart] == 6): vecParticleProcName.append('Annihilation')
				if (temp_proc_id[jpart] == 7): vecParticleProcName.append('Conversion')
				if (temp_proc_id[jpart] == 8): vecParticleProcName.append('hIonization')
			

frac_part = 100.*(float(countPart)/float(countDet))


part_e_min = np.min(vecParticleEnergy)
part_e_max = np.max(vecParticleEnergy)
n_bins = int((part_e_max - part_e_min)/E_bin_part)

# histograms for the slab
N_array_part, bin_array_part = np.histogram(vecParticleEnergy, bins=n_bins)
flux_array_part = np.zeros(len(N_array_part))
err_flux_array_part = np.zeros(len(N_array_part))
err_ene_part = np.zeros(len(N_array_part))
ene_array_part = np.zeros(len(N_array_part))
left_ene_array_part = np.zeros(len(N_array_part))

for jn in xrange(len(N_array_part)):
    err_ene_part[jn] = (bin_array_part[jn+1] - bin_array_part[jn])/2.
    ene_array_part[jn] = bin_array_part[jn] + err_ene_part[jn]
    left_ene_array_part[jn] = bin_array_part[jn]
    energy_bin_part = bin_array_part[jn+1] - bin_array_part[jn]
    flux_array_part[jn] = float(N_array_part[jn])/(Time*energy_bin_part*det_area)
    err_flux_array_part[jn] = (np.sqrt(float(N_array_part[jn])))/(Time*energy_bin_part*det_area)

det_e_min = np.min(vecDetEnergy)
det_e_max = np.max(vecDetEnergy)
n_bins = int((det_e_max - det_e_min)/E_bin)

# histograms for the slab
N_array_det, bin_array_det = np.histogram(vecDetEnergy, bins=n_bins)
flux_array_det = np.zeros(len(N_array_det))
err_flux_array_det = np.zeros(len(N_array_det))
err_ene_det = np.zeros(len(N_array_det))
ene_array_det = np.zeros(len(N_array_det))
left_ene_array_det = np.zeros(len(N_array_det))

for jn in xrange(len(N_array_det)):
    err_ene_det[jn] = (bin_array_det[jn+1] - bin_array_det[jn])/2.
    ene_array_det[jn] = bin_array_det[jn] + err_ene_det[jn]
    left_ene_array_det[jn] = bin_array_det[jn]
    energy_bin_det = bin_array_det[jn+1] - bin_array_det[jn]
    flux_array_det[jn] = float(N_array_det[jn])/(Time*energy_bin_det*det_area)
    err_flux_array_det[jn] = (np.sqrt(float(N_array_det[jn])))/(Time*energy_bin_det*det_area)



count_ProcName = Counter(vecParticleProcName)
ProcNames = count_ProcName.keys()
ProcValues = count_ProcName.values()
ProcExplode = []
ProcColor = []
colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'orange', 'thistle', 'lightgrey', 'tomato']
for jproc in xrange(len(ProcValues)):
	ProcExplode.append(0.05)
	ProcColor.append(colors[jproc])


# Plot the results
fig1 = plt.figure(1, figsize=(10, 6))
ax1 = fig1.add_subplot(111)

# part

if type == 0:
	ax1.bar(left_ene_array_part, flux_array_part, width=2.*err_ene_part, edgecolor="yellow", facecolor='yellow', lw = 2, label=str(len(vecParticleEnergy))+' '+type_name+' ('+str(round(frac_part,1))+'% of bkg counts)', alpha=0.3)
	ax1.errorbar(ene_array_part, flux_array_part, xerr=err_ene_part, yerr=err_flux_array_part, capsize=0, fmt='o', color='yellow', lw = 2, ecolor='black')
if type == 1:
	ax1.bar(left_ene_array_part, flux_array_part, width=2.*err_ene_part, edgecolor="blue", facecolor='blue', lw = 2, label=str(len(vecParticleEnergy))+' '+type_name+' ('+str(round(frac_part,1))+'% of bkg counts)', alpha=0.3)
	ax1.errorbar(ene_array_part, flux_array_part, xerr=err_ene_part, yerr=err_flux_array_part, capsize=0, fmt='o', color='blue', lw = 2, ecolor='black')


ax1.set_title(source_name+' induced '+type_name+' flux')
ax1.set_xlabel("Energy [keV]")
ax1.set_ylabel("Part. cm$^{-2}$ s$^{-1}$ keV$^{-1}$")
ax1.set_ylim(0.0000001, 0.01)
#ax.set_xscale('log')
ax1.set_yscale('log')

ax1.legend(numpoints=1)


# Plot the results
fig2 = plt.figure(2, figsize=(10, 6))
ax2 = fig2.add_subplot(111)

# part

ax2.bar(left_ene_array_det, flux_array_det, width=2.*err_ene_det, edgecolor="lawngreen", facecolor='lawngreen', lw = 2,label='T = '+str(round(Time/1000., 1))+' ks', alpha=1)
ax2.errorbar(ene_array_det, flux_array_det, xerr=err_ene_det, yerr=err_flux_array_det, capsize=0, fmt='.', color='lawngreen', lw = 2, ecolor='black')


ax2.set_title(source_name+' induced background rate')
ax2.set_xlabel("Energy [keV]")
ax2.set_ylabel("cts cm$^{-2}$ s$^{-1}$ keV$^{-1}$")
ax2.set_ylim(0.0000001, 0.01)
#ax.set_xscale('log')
ax2.set_yscale('log')

ax2.legend(numpoints=1)

# Plot the results
fig3 = plt.figure(3, figsize=(6, 6))
ax3 = fig3.add_subplot(111)

ax3.pie(ProcValues, explode=ProcExplode, labels=ProcNames, colors=ProcColor, autopct='%1.1f%%', shadow=True, startangle=90)
ax3.set_title(type_name+' origin process')


hdulist_part.close()
hdulist_det.close()
plt.show()
