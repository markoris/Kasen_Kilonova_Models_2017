import h5py, glob
import numpy as np

def get_spectrum(spec_file):
	nu = np.array(spec_file['nu'],dtype='d')
	
	# specific luminosity (ergs/s/Hz) 
	# this is a 2D array, Lnu[times][nu]
	Lnu_all   = np.array(spec_file['Lnu'],dtype='d')

	# if you want thing in Flambda (ergs/s/Angstrom)
	c = 2.99e10
	lam  = c/nu*1e8
	Llam = Lnu_all*nu**2.0/c/1e8
	
	# divide by 10 pc to get absolute spectrum

	pc = 3.086e18 # parsec in cm
	D = 10*pc # 10 parsecs
	Llam / (4 * np.pi * D**2)
	
	return Llam

dyn_models_xlan_hi = glob.glob('systematic_kilonova_model_grid/*Xlan1e-1*')
dyn_models_xlan_lo = [fname.replace('Xlan1e-1.0', 'Xlan1e-2.0') for fname in dyn_models_xlan_hi]
wind_models = [fname.replace('Xlan1e-1.0', 'Xlan1e-9.0') for fname in dyn_models_xlan_hi]


for i in range(len(dyn_models_xlan_hi)):
	# load higher X_La fraction model (X_La = 0.1)
	dyn_xlan_hi = h5py.File(dyn_models_xlan_hi[i], 'r')
	# load lower X_La fraction model (X_La = 0.01)
	dyn_xlan_lo = h5py.File(dyn_models_xlan_lo[i], 'r')
	# get their spectra
	spec_xlan_hi = get_spectrum(dyn_xlan_hi)
	spec_xlan_lo = get_spectrum(dyn_xlan_lo)
	# take average, corresponding to X_La = 0.055, close to ~0.053 for LANL dynamical ejecta
	spec_dyn = (spec_xlan_hi + spec_xlan_lo)/2

	# TODO: add for loop over Xlan 1e-9 components, add spectra together
	# TODO: convolve spectra with filters, save all 3000 sims to an hdf5 file with individual keys for each band (grizyJHK)

	for j in range(len(wind_models)):
		wind = h5py.File(wind_models[j], 'r')
		spec_wind = get_spectrum(wind)

		# ADD SPECTRA FLUXES TOGETHER
		# note: this ignores all potential photon + chemical reprocessing
		spec_2c = spec_dyn + spec_wind	
	
		params_dyn = dyn_models_xlan_hi[i].split('_')[3:5]
		params_wind = wind_models[j].split('_')[3:5]
		
		params = 'Run_SS_Xlan0.055_Xlan1e-9_md%f_vd%f_mw%f_vw%f.dat 

	# times confirmed matching	
	times_xlan_hi = np.array(dyn_xlan_hi['time'])
	times_xlan_lo = np.array(dyn_xlan_lo['time'])
	times_xlan_hi /= 3600.0/24.0
	times_xlan_lo /= 3600.0/24.0

