import h5py, glob
import numpy as np

def get_spectrum(spec_file):
	nu = np.array(spec_file['nu'],dtype='d')
	
	# specific luminosity (ergs/s/Hz) 
	# this is a 2D array, Lnu[times][nu]
	Lnu_all   = np.array(spec_file['Lnu'],dtype='d')

	# if you want thing in Flambda (ergs/s/Angstrom)
	c = 2.99e10 # cm/s
	lam  = c/nu*1e8 # Angstroms
	Llam = Lnu_all*nu**2.0/c/1e8
	
	# divide by 10 pc to get absolute spectrum

	pc = 3.086e18 # parsec in cm
	D = 10*pc # 10 parsecs
	Llam / (4 * np.pi * D**2)
	
	return Llam, lam

# NOTE: one of the high-lanthanide fraction models is missing from the grid: m0.1_vk0.30_fd1.0_Xlan1e-1 MISSING!
dyn_models_xlan_hi = glob.glob('systematic_kilonova_model_grid/*Xlan1e-1*')
dyn_models_xlan_lo = [fname.replace('Xlan1e-1.0', 'Xlan1e-2.0') for fname in dyn_models_xlan_hi]
wind_models = [fname.replace('Xlan1e-1.0', 'Xlan1e-9.0') for fname in dyn_models_xlan_hi]

n_sims = len(dyn_models_xlan_hi)*len(wind_models)

for i in range(len(dyn_models_xlan_hi)//10):
	# load higher X_La fraction model (X_La = 0.1)
	dyn_xlan_hi = h5py.File(dyn_models_xlan_hi[i], 'r')
	# load lower X_La fraction model (X_La = 0.01)
	dyn_xlan_lo = h5py.File(dyn_models_xlan_lo[i], 'r')
	# get their spectra
	spec_xlan_hi, lam = get_spectrum(dyn_xlan_hi)
	spec_xlan_lo, _ = get_spectrum(dyn_xlan_lo)
	# take average, corresponding to X_La = 0.055, close to ~0.053 for LANL dynamical ejecta
	spec_dyn = (spec_xlan_hi + spec_xlan_lo)/2

	times = np.array(dyn_xlan_lo['time'])
	times = np.round(times/3600.0/24.0, 2) # days

	# TODO: convolve spectra with filters, save all 3000 sims to an hdf5 file with individual keys for each band (grizyJHK)

	for j in range(len(wind_models)//10):
		wind = h5py.File(wind_models[j], 'r')
		spec_wind, _ = get_spectrum(wind)

		# ADD SPECTRA FLUXES TOGETHER
		# note: this ignores all potential photon + chemical reprocessing
		spec_2c = spec_dyn + spec_wind

		# take only spectral data where time is positive
		time_positive = np.where(times > 0)[0]
		spec_2c = spec_2c[time_positive, :]
	
		params_dyn = dyn_models_xlan_hi[i].split('/')[1].split('_')[3:5]
		params_wind = wind_models[j].split('/')[1].split('_')[3:5]
		
		# Naming convention to match LANL simulation filenames
		#params = 'Run_SS_Xlan0.055_Xlan1e-9_md%s_vd%s_mw%s_vw%s.dat' % (params_dyn[0][1:], params_dyn[1][2:], params_wind[0][1:], params_wind[1][2:])
		params = np.array([params_dyn[0][1:], params_dyn[1][2:], params_wind[0][1:], params_wind[1][2:]], dtype=float)

		try:
			params_array = np.concatenate((params_array, params[None, :]), axis=0)
			spectra_array = np.concatenate((spectra_array, spec_2c[None, :]), axis=0)
		except NameError:
			params_array = params[None, :]
			spectra_array = spec_2c[None, :]	

	print('done with %d / %d dyn models' % (i*54+j, 54*55))

# save params + spectra arrays for future use
print('writing spectra to hdf5...')
h5f = h5py.File('hdf5_data/spectra.h5', 'a')
header = h5f.create_dataset('header', (1))
header.attrs['units_wavelength'] = 'Angstrom'
header.attrs['units_spectral_energy_density'] = 'erg / s / cm^2 / Angstrom'
header.attrs['spectra_array_structure'] = '[N, 250, 1629]: N = number of simulations, 250 = observation times, 1629 = spectral energy density in relevant wavelength bin'
header.attrs['observation_times'] = np.array2string(times, precision=4, max_line_width=50)
header.attrs['units_observation_times'] = 'days'
h5f.create_dataset('params', data=params_array, compression="gzip", chunks=True, maxshape=(None, 4)) 
h5f.create_dataset('spectra', data=spectra_array, compression="gzip", chunks=True, maxshape=(None, 250, 1629))
h5f.close()

# convolve with observational filters

from conversions import spec_to_mags

bands = 'grizyJHK'

lcs =  {'g': [],
	'r': [],
	'i': [],
	'z': [],
	'y': [],
	'J': [],
	'H': [],
	'K': []}
# for each spectrum over time and wavelength
for i in range(spectra_array.shape[0]):
	wav, spec = np.copy(lam), np.copy(spectra_array[i])
	# for each filter, create a light curve ...
	lc =   {'g': [],
		'r': [],
		'i': [],
		'z': [],
		'y': [],
		'J': [],
		'H': [],
		'K': []}
	for b in bands:
		# ... by converting the spectra at each time to a broadband filter magnitude
		for t in range(spec.shape[0]):
			flux_band = spec_to_mags(wav, spec[t, :], band=b)
			lc[b] = np.concatenate((lc[b], flux_band[0]), axis=0)
	lcs[b] = np.concatenate((lcs[b], lc[b]), axis=0)
print(lcs)
