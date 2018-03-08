import numpy as np
import matplotlib.pyplot as plt
import isodist
import os
import glob
import pickle
from tqdm import tqdm

ISO_PATH = os.getenv('ISODIST_DATA')

#load isochrones
path = ISO_PATH+'/parsec-sdss-2mass/'
Zs = [float(file[79:85]) for file in glob.glob(os.path.join(path, '*.dat.gz'))]

isochrone = isodist.PadovaIsochrone(Z=Zs, parsec=True)
logages = isochrone.logages()
logages = logages[logages >= 8.]
Zs = isochrone.Zs()

#Generate isochrone grid
def isochrone_grid(logages, Zs, imf):
	grid = []
	for lage in tqdm(range(len(logages))):
		for met in Zs:
			#print 'generating grid entry for Age='+str(round((10**lage)/1e9,1))+'Gyr and Fe/H='+str(round(isodist.Z2FEH(met),2))+'...'
			iso = isochrone(logages[lage], Z=met)
			H_iso = iso['H'][1:]
			M_iso = iso['M_ini'][1:]
			logg_iso = iso['logg'][1:]
			J_K_iso = iso['J'][1:]-iso['Ks'][1:]
			lage_iso = np.ones(len(iso['H'][1:]))*logages[lage]-9
			met_iso = np.ones(len(iso['H'][1:]))*met
			m_weights = (imf(iso['M_ini'][1:], int=True)-imf(iso['M_ini'][:-1], int=True))
			entry = np.dstack((lage_iso, met_iso, M_iso, logg_iso, H_iso, J_K_iso, m_weights))[0]
			grid.extend(entry)
	grid = np.array(grid)
	return grid

imf= isodist.imf.chabrier2003
grid = isochrone_grid(logages, Zs, imf)
savfile = open('../savs/Padova_grid_chabrier2003.sav', 'w')
pickle.dump(grid, savfile)
savfile.close()


imf= isodist.imf.lognormalChabrier2001
grid = isochrone_grid(logages, Zs, imf)
savfile = open('../savs/Padova_grid_lognormalchabrier2001.sav', 'w')
pickle.dump(grid, savfile)
savfile.close()

