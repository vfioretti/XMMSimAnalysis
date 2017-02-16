"""
 filterParticleFlux.py  -  description
 ---------------------------------------------------------------------------------
 filtering the particle flux entering the detector after 
 interaction with XMM proton shield
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python filterParticleFlux.py filedir N_file N_in source
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - N_in = number of simulated particles (in milions)
 - particle source = [0 = CXB]
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2016/10/31: creation date
"""

from astropy.io import fits
import numpy as np
import math
import sys, os


# Import the input parameters
arg_list = sys.argv
filedir = arg_list[1]
N_fits = int(arg_list[2])
N_in = int(arg_list[3])
source = int(arg_list[4])

# Photons
vecPhotEnergy = []
vecPhotEventID = []
vecPhotCountEnergy = []
vecPhotPosX = []
vecPhotPosY = []
vecPhotProcID = []
vecPhotPrimaryFlag = []

# Electrons
vecElecEnergy = []
vecElecEventID = []
vecElecCountEnergy = []
vecElecPosX = []
vecElecPosY = []
vecElecProcID = []
vecElecPrimaryFlag = []

# Positrons
vecPosEnergy = []
vecPosEventID = []
vecPosCountEnergy = []
vecPosPosX = []
vecPosPosY = []
vecPosProcID = []
vecPosPrimaryFlag = []

# Total detector counts
vecCountEnergy = []
vecEventID = []

# Volume set-up
shield_id = 10
det_id = 20
box_id = 1000
hat_id = 2000
sphere_id = 3000

if (source == 0): source_name = 'CXB'

for jfits in xrange(N_fits):
    
    print '%%%%%%%%%%%%%% READING THELSim FILE: ', jfits+1
    hdulist = fits.open(filedir+'/xyz.'+str(jfits)+'.fits.gz')
    
    tbdata = hdulist[1].data
    
    evt_id = tbdata.field('EVT_ID')
    vol_id = tbdata.field('VOLUME_ID')
    track_id = tbdata.field('TRK_ID')
    parent_id = tbdata.field('PARENT_TRK_ID')
    gtime = tbdata.field('GTIME_ENT')
    e_part = tbdata.field('E_KIN_ENT')
    part_id = tbdata.field('PARTICLE_ID')
    e_dep = tbdata.field('E_DEP')
    pos_dep_x = tbdata.field('X_ENT')
    pos_dep_y = tbdata.field('Y_ENT')
    proc_id = tbdata.field('PROCESS_ID')

    # ---------------------------------

    index = 0
    where_sameevent = [index]
    temp_index = index
    while (where_sameevent[-1] < (len(evt_id)-1)):
			while evt_id[temp_index] == evt_id[temp_index+1]:
				where_sameevent.append(temp_index+1)
				temp_index += 1
				print temp_index
			else:
				sameev_event_id = evt_id[where_sameevent]
				sameev_vol_id = vol_id[where_sameevent]
				sameev_track_id = track_id[where_sameevent]
				sameev_parent_id = parent_id[where_sameevent]
				sameev_gtime = gtime[where_sameevent]
				sameev_e_part = e_part[where_sameevent]
				sameev_part_id = part_id[where_sameevent]
				sameev_e_dep = e_dep[where_sameevent]
				sameev_pos_dep_x = pos_dep_x[where_sameevent]
				sameev_pos_dep_y = pos_dep_y[where_sameevent]
				sameev_proc_id = proc_id[where_sameevent]
		
				sameev_together = zip(sameev_gtime, sameev_vol_id, sameev_e_part)
				sort_time = sorted(sameev_together)
		
				gtime_sorted = [x[0] for x in sort_time]
				vol_id_sorted = [x[1] for x in sort_time]
				e_part_sorted = [x[2] for x in sort_time]
				
				if  (vol_id_sorted[0] == sphere_id):
					first_energy = e_part_sorted[0]
					
					gtime_sorted = np.array(gtime_sorted)
					vol_id_sorted = np.array(vol_id_sorted)
					e_part_sorted = np.array(e_part_sorted)

					gtime_sorted = gtime_sorted[vol_id_sorted != sphere_id]
					e_part_sorted = e_part_sorted[vol_id_sorted != sphere_id]
					vol_id_sorted = vol_id_sorted[vol_id_sorted != sphere_id]				
				else:
					print "No touching the sphere!!!!!!"
				
				if (vol_id_sorted.size):
				 if (vol_id_sorted[0] == shield_id):
					where_det = np.where(sameev_vol_id == det_id)
					if (where_det[0].size):
				 		track_id_shield = sameev_track_id[np.where(sameev_vol_id == shield_id)]
						e_dep_det = sameev_e_dep[where_det]
						e_part_det = sameev_e_part[where_det]
						part_id_det = sameev_part_id[where_det]
						track_id_det = sameev_track_id[where_det]
						pos_dep_x_det = sameev_pos_dep_x[where_det]
						pos_dep_y_det = sameev_pos_dep_y[where_det]
						proc_id_det = sameev_proc_id[where_det]
						track_id_det = sameev_track_id[where_det]
						parent_id_det = sameev_parent_id[where_det]
						sum_edep = 0.
						for jdet in xrange(len(where_det[0])):
							if e_dep_det[jdet] != 0.:
								sum_edep = sum_edep + e_dep_det[jdet]
						if (sum_edep != 0):
							vecCountEnergy.append(sum_edep)
							vecEventID.append(evt_id[index])
						
							for jdet in xrange(len(where_det[0])):								
								where_from_shield = np.where(track_id_shield == track_id_det[jdet])
								if (where_from_shield[0].size):	
									flag_hit = 0
									if e_part_det[jdet] == first_energy:
										leakage = 1
									else:
										leakage = 0
										
									if (e_dep_det[jdet] == 0):
										track_id_primary = track_id_det[jdet]
										ene_shower = e_dep_det[np.where(parent_id_det == track_id_primary)]
										hit = [x for x in ene_shower if x > 0.0]
										if hit: flag_hit = 1
									else:
										flag_hit = 1
									
									if flag_hit:
										if part_id_det[jdet] == 22:
											vecPhotEnergy.append(e_part_det[jdet])
											vecPhotEventID.append(evt_id[index])
											vecPhotCountEnergy.append(e_dep_det[jdet])
											vecPhotPosX.append(pos_dep_x_det[jdet])
											vecPhotPosY.append(pos_dep_y_det[jdet])
											vecPhotProcID.append(proc_id_det[jdet])
											vecPhotPrimaryFlag.append(leakage)
										else:
											if part_id_det[jdet] == 11:
												vecElecEnergy.append(e_part_det[jdet])
												vecElecEventID.append(evt_id[index])
												vecElecCountEnergy.append(e_dep_det[jdet])
												vecElecPosX.append(pos_dep_x_det[jdet])
												vecElecPosY.append(pos_dep_y_det[jdet])
												vecElecProcID.append(proc_id_det[jdet])
												vecElecPrimaryFlag.append(leakage)

											else:
												if part_id_det[jdet] == -11:
													vecPosEnergy.append(e_part_det[jdet])
													vecPosEventID.append(evt_id[index])
													vecPosCountEnergy.append(e_dep_det[jdet])
													vecPosPosX.append(pos_dep_x_det[jdet])
													vecPosPosY.append(pos_dep_y_det[jdet])
													vecPosProcID.append(proc_id_det[jdet])
													vecPosPrimaryFlag.append(leakage)

												else:
													print "UNKNOWN PARTICLE!!!!!!!!!!!!!!!!!!!!!!"

				N_event_eq = len(where_sameevent)
				if (evt_id[where_sameevent[-1]+1] != evt_id[-1]):
					index = where_sameevent[N_event_eq-1] + 1
					where_sameevent = [index]
					temp_index = index

				else:
					first_sameevent = where_sameevent[-1]+1
					where_sameevent = np.arange(first_sameevent, len(evt_id), 1)

				 	sameev_event_id = evt_id[where_sameevent]
					sameev_vol_id = vol_id[where_sameevent]
					sameev_track_id = track_id[where_sameevent]
					sameev_parent_id = parent_id[where_sameevent]
					sameev_gtime = gtime[where_sameevent]
					sameev_e_part = e_part[where_sameevent]
					sameev_part_id = part_id[where_sameevent]
					sameev_e_dep = e_dep[where_sameevent]
					sameev_pos_dep_x = pos_dep_x[where_sameevent]
					sameev_pos_dep_y = pos_dep_y[where_sameevent]
					sameev_proc_id = proc_id[where_sameevent]
					
					sameev_together = zip(sameev_gtime, sameev_vol_id, sameev_e_part)
					sort_time = sorted(sameev_together)
					
					gtime_sorted = [x[0] for x in sort_time]
					vol_id_sorted = [x[1] for x in sort_time]
					e_part_sorted = [x[2] for x in sort_time]
					
					if  (vol_id_sorted[0] == sphere_id):
					   first_energy = e_part_sorted[0]
					   	
					   gtime_sorted = np.array(gtime_sorted)
					   vol_id_sorted = np.array(vol_id_sorted)
					   e_part_sorted = np.array(e_part_sorted)

					   gtime_sorted = gtime_sorted[vol_id_sorted != sphere_id]
					   e_part_sorted = e_part_sorted[vol_id_sorted != sphere_id]
					   vol_id_sorted = vol_id_sorted[vol_id_sorted != sphere_id]				
					else:
					   print "No touching the sphere!!!!!!"
					
					if (vol_id_sorted.size):
					 if (vol_id_sorted[0] == shield_id):
						where_det = np.where(sameev_vol_id == det_id)
						if (where_det[0].size):
							track_id_shield = sameev_track_id[np.where(sameev_vol_id == shield_id)]
							e_dep_det = sameev_e_dep[where_det]
							e_part_det = sameev_e_part[where_det]
							part_id_det = sameev_part_id[where_det]
							track_id_det = sameev_track_id[where_det]
							pos_dep_x_det = sameev_pos_dep_x[where_det]
							pos_dep_y_det = sameev_pos_dep_y[where_det]
							proc_id_det = sameev_proc_id[where_det]
							track_id_det = sameev_track_id[where_det]
							parent_id_det = sameev_parent_id[where_det]
							sum_edep = 0.
							for jdet in xrange(len(where_det[0])):
								if e_dep_det[jdet] != 0.:
									sum_edep = sum_edep + e_dep_det[jdet]
							if (sum_edep != 0):
								vecCountEnergy.append(sum_edep)
								vecEventID.append(evt_id[index])
						
								for jdet in xrange(len(where_det[0])):								
									where_from_shield = np.where(track_id_shield == track_id_det[jdet])
									if (where_from_shield[0].size):	
										flag_hit = 0
										if e_part_det[jdet] == first_energy:
										   leakage = 1
										else:
										   leakage = 0
										   
										if (e_dep_det[jdet] == 0):
											track_id_primary = track_id_det[jdet]
											ene_shower = e_dep_det[np.where(parent_id_det == track_id_primary)]
											hit = [x for x in ene_shower if x > 0.0]
											if hit: flag_hit = 1
										else:
											flag_hit = 1
										
										if flag_hit:
											if part_id_det[jdet] == 22:
												vecPhotEnergy.append(e_part_det[jdet])
												vecPhotEventID.append(evt_id[index])
												vecPhotCountEnergy.append(e_dep_det[jdet])
												vecPhotPosX.append(pos_dep_x_det[jdet])
												vecPhotPosY.append(pos_dep_y_det[jdet])
												vecPhotProcID.append(proc_id_det[jdet])
												vecPhotPrimaryFlag.append(leakage)
											else:
												if part_id_det[jdet] == 11:
													vecElecEnergy.append(e_part_det[jdet])
													vecElecEventID.append(evt_id[index])
													vecElecCountEnergy.append(e_dep_det[jdet])
													vecElecPosX.append(pos_dep_x_det[jdet])
													vecElecPosY.append(pos_dep_y_det[jdet])
													vecElecProcID.append(proc_id_det[jdet])
													vecElecPrimaryFlag.append(leakage)
												else:
													if part_id_det[jdet] == -11:
														vecPosEnergy.append(e_part_det[jdet])
														vecPosEventID.append(evt_id[index])
														vecPosCountEnergy.append(e_dep_det[jdet])
														vecPosPosX.append(pos_dep_x_det[jdet])
														vecPosPosY.append(pos_dep_y_det[jdet])
														vecPosProcID.append(proc_id_det[jdet])
														vecPosPrimaryFlag.append(leakage)
													else:
													    print "UNKNOWN PARTICLE!!!!!!!!!!!!!!!!!!!!!!"

					hdulist.close()
					break

# write FITS for each product

path_output = './'+source_name+'/'+str(N_in)+'mil/'
if not os.path.exists(path_output):
	os.makedirs(path_output)


# Photons
Phot_col1 = fits.Column(name='EVENT_ID', format='1J', array=np.array(vecPhotEventID))
Phot_col2 = fits.Column(name='ENERGY', format='1D', array=np.array(vecPhotEnergy))
Phot_col3 = fits.Column(name='E_DEP', format='1D', array=np.array(vecPhotCountEnergy))
Phot_col4 = fits.Column(name='POS_X', format='1D', array=np.array(vecPhotPosX))
Phot_col5 = fits.Column(name='POS_Y', format='1D', array=np.array(vecPhotPosY))
Phot_col6 = fits.Column(name='PROC_ID', format='1J', array=np.array(vecPhotProcID))
Phot_col7 = fits.Column(name='PRIMARY_FLAG', format='1J', array=np.array(vecPhotPrimaryFlag))

Phot_cols = fits.ColDefs([Phot_col1, Phot_col2, Phot_col3, Phot_col4, Phot_col5, Phot_col6, Phot_col7])
Phot_tbhdu = fits.BinTableHDU.from_columns(Phot_cols)
Phot_tbhdu.writeto(path_output+source_name+'_'+str(N_in)+'mil_XMMSim_Photons.fits.gz', clobber=1)

# Electrons
Elec_col1 = fits.Column(name='EVENT_ID', format='1J', array=np.array(vecElecEventID))
Elec_col2 = fits.Column(name='ENERGY', format='1D', array=np.array(vecElecEnergy))
Elec_col3 = fits.Column(name='E_DEP', format='1D', array=np.array(vecElecCountEnergy))
Elec_col4 = fits.Column(name='POS_X', format='1D', array=np.array(vecElecPosX))
Elec_col5 = fits.Column(name='POS_Y', format='1D', array=np.array(vecElecPosY))
Elec_col6 = fits.Column(name='PROC_ID', format='1J', array=np.array(vecElecProcID))
Elec_col7 = fits.Column(name='PRIMARY_FLAG', format='1J', array=np.array(vecElecPrimaryFlag))

Elec_cols = fits.ColDefs([Elec_col1, Elec_col2, Elec_col3, Elec_col4, Elec_col5, Elec_col6, Elec_col7])
Elec_tbhdu = fits.BinTableHDU.from_columns(Elec_cols)
Elec_tbhdu.writeto(path_output+source_name+'_'+str(N_in)+'mil_XMMSim_Electrons.fits.gz', clobber=1)

# Positrons
Pos_col1 = fits.Column(name='EVENT_ID', format='1J', array=np.array(vecPosEventID))
Pos_col2 = fits.Column(name='ENERGY', format='1D', array=np.array(vecPosEnergy))
Pos_col3 = fits.Column(name='E_DEP', format='1D', array=np.array(vecPosCountEnergy))
Pos_col4 = fits.Column(name='POS_X', format='1D', array=np.array(vecPosPosX))
Pos_col5 = fits.Column(name='POS_Y', format='1D', array=np.array(vecPosPosY))
Pos_col6 = fits.Column(name='PROC_ID', format='1J', array=np.array(vecPosProcID))
Pos_col7 = fits.Column(name='PRIMARY_FLAG', format='1J', array=np.array(vecPosPrimaryFlag))

Pos_cols = fits.ColDefs([Pos_col1, Pos_col2, Pos_col3, Pos_col4, Pos_col5, Pos_col6, Pos_col7])
Pos_tbhdu = fits.BinTableHDU.from_columns(Pos_cols)
Pos_tbhdu.writeto(path_output+source_name+'_'+str(N_in)+'mil_XMMSim_Positrons.fits.gz', clobber=1)

# Detector
Det_col1 = fits.Column(name='EVENT_ID', format='1J', array=np.array(vecEventID))
Det_col2 = fits.Column(name='E_DEP', format='1D', array=np.array(vecCountEnergy))

Det_cols = fits.ColDefs([Det_col1, Det_col2])
Det_tbhdu = fits.BinTableHDU.from_columns(Det_cols)
Det_tbhdu.writeto(path_output+source_name+'_'+str(N_in)+'mil_XMMSim_Detector.fits.gz', clobber=1)

