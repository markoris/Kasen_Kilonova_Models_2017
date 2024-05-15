import numpy, h5py
import matplotlib.pyplot as plt
import numpy as np

band = 'K'

lc = h5py.File('hdf5_data/lcs_%s.h5' % band)

times = np.array(lc['header'].attrs['observation_times'][1:-1].split(), dtype=float)
time_positive = np.where(times > 0)[0]
times = times[time_positive]

lcs = lc['spectra']

for i in range(lcs.shape[0]):
	plt.plot(times, lcs[i], alpha=0.1)
	if i % 100 == 0: print(i)

plt.xscale('log')
plt.ylim(-40, -20)
plt.gca().invert_yaxis()

plt.savefig('lcs_%s.png' % band)
