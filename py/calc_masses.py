import numpy
import numpy as np
import densprofiles
import define_rgbsample
import pickle
from isodist import Z2FEH
from galpy.util import bovy_coords
from fitDens import _setup_densfunc
import os
from scipy.integrate import quad
from scipy import interpolate
import multiprocessing
from galpy.util import multi

def load_isochrones(gridfile):
    iso_file = open(gridfile)
    iso_grid = pickle.load(iso_file)
    return iso_grid

def calc_normalisation(params, nbin, iso_grid,
                       fehbin=[-0.1,0.0], 
                       agebin=[1.,3.], 
                       loggcut=[1.8,3.0],
                       teffcut=[4000,5000], 
                       type='brokenexpflare',
                       verbose=True,
                       fitIndx=None,
                       weights = 'padova',
                       distance_cut = False,
                       lowermass = None):
    #first get the values necessary from the isochrone grid
    #make a mask for giant stars (+ J-K cut)
    if teffcut == None:
        giants = (iso_grid[:,3] >= loggcut[0])&(iso_grid[:,3] < loggcut[1])&(iso_grid[:,5] > 0.5)
    else:
        giants = (iso_grid[:,3] >= loggcut[0])&(iso_grid[:,3] < loggcut[1])&(iso_grid[:,5] > 0.5)&(10**iso_grid[:,7] >= teffcut[0])&(10**iso_grid[:,7] < teffcut[1])
    #make a mask for the age and feh bin
    if agebin == None:
    	bin = (10**iso_grid[:,0] >= 0.)&(10**iso_grid[:,0] < 13.)&\
          (Z2FEH(iso_grid[:,1]) >= fehbin[0])&(Z2FEH(iso_grid[:,1]) < fehbin[1])
    else:
    	bin = (10**iso_grid[:,0] >= agebin[0])&(10**iso_grid[:,0] < agebin[1])&\
          (Z2FEH(iso_grid[:,1]) >= fehbin[0])&(Z2FEH(iso_grid[:,1]) < fehbin[1])
    
    if lowermass != None:
		giants *= iso_grid[:,2] >= lowermass
		bin *= iso_grid[:,2] >= lowermass
    if len(iso_grid[:,0][bin]) < 1:
        fehs = np.unique(Z2FEH(iso_grid[:,1]))
        cfehbin = fehbin[0]+((fehbin[1]-fehbin[0])/2)
        feh_offsets = np.fabs(fehs-cfehbin)
        ind = np.argmin(feh_offsets)
        cfeh = fehs[ind]
        bin = (10**iso_grid[:,0] >= agebin[0])&(10**iso_grid[:,0] < agebin[1])&\
          (Z2FEH(iso_grid[:,1]) == cfeh)
    #find the average giant mass
    mass = iso_grid[:,2]
    if weights == 'padova':
        weight = iso_grid[:,6]*(10**iso_grid[:,0]/iso_grid[:,1])
    if weights == 'basti':
        weight = iso_grid[:,6]
    av_mass = np.sum(mass[giants&bin]*weight[giants&bin])/np.sum(weight[giants&bin])
    #find the ratio between giants and the total stellar pop. for this bin
    mass_total = mass[bin]
    weight_total = weight[bin]
    mass_bin = mass[giants&bin]
    weight_bin = weight[giants&bin]
    m_ratio = np.sum(mass_bin*weight_bin)/np.sum(mass_total*weight_total)
    #now compute and sum the rate for this density function
    #load the raw selection function
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    #load the effective selection function
	if agebin == None:
		with open('../essf/maps/essf_rgb_green15_modelmh_feh'+str(round(fehbin[0],1))+'.sav','rb') as savefile:
			locations= pickle.load(savefile)
			effsel= pickle.load(savefile)
			distmods= pickle.load(savefile)
		with open('../essf/maps/essf_rgb_marshall06_modelmh_feh'+str(round(fehbin[0],1))+'.sav','rb') as savefile:
			mlocations= pickle.load(savefile)
			meffsel= pickle.load(savefile)
			mdistmods= pickle.load(savefile)
    if agebin != None:
        if agebin[0] < 1.:
            with open('../essf/maps/essf_rgb_green15_modelmh_feh'+str(round(fehbin[0],1))+'_age'+str(round(1.0,1))+'.sav','rb') as savefile:
                locations= pickle.load(savefile)
                effsel= pickle.load(savefile)
                distmods= pickle.load(savefile)
            with open('../essf/maps/essf_rgb_marshall06_modelmh_feh'+str(round(fehbin[0],1))+'_age'+str(round(1.0,1))+'.sav','rb') as savefile:
                mlocations= pickle.load(savefile)
                meffsel= pickle.load(savefile)
                mdistmods= pickle.load(savefile)
        if agebin[0] > 0.9:
            with open('../essf/maps/essf_rgb_green15_modelmh_feh'+str(round(fehbin[0],1))+'_age'+str(round(agebin[0],1))+'.sav','rb') as savefile:
                locations= pickle.load(savefile)
                effsel= pickle.load(savefile)
                distmods= pickle.load(savefile)
            with open('../essf/maps/essf_rgb_marshall06_modelmh_feh'+str(round(fehbin[0],1))+'_age'+str(round(agebin[0],1))+'.sav','rb') as savefile:
                mlocations= pickle.load(savefile)
                meffsel= pickle.load(savefile)
                mdistmods= pickle.load(savefile)
			
    # Fill in regions not covered by Marshall map
    meffsel[meffsel < -0.5]= effsel[meffsel < -0.5]
    if fitIndx is None:
            fitIndx= numpy.ones(len(mlocations),dtype='bool') #True-betwDiskIndx 
    locations, effsel, distmods = np.array(mlocations)[fitIndx], np.array(meffsel)[fitIndx], mdistmods
    #get the density function and set it up to find the normalisation (surfdens=True)
    rdensfunc= _setup_densfunc(type)
    densfunc= lambda x: rdensfunc(x,None,None,params=params, surfdens=True)
    #evaluate surface density at R0 for the density normalisation (always 1. if R_b > R0)
    R0 = densprofiles._R0
    Rb = np.exp(params[3])
    dens_norm = densfunc(densprofiles._R0)
    #set up the density function again with surfdens=False for the rate calculation
    rdensfunc= _setup_densfunc(type)
    densfunc= lambda x,y,z: rdensfunc(x,y,z,params=params, surfdens=False)
    ds= 10.**(distmods/5.-2.)
    #imply the distance cut if distance_cut == True
    if distance_cut == True:
    	distmods = distmods[ds <= 3.]
    	ds= ds[ds <= 3.]
    	effsel = effsel[:,:len(ds)]
    #Compute the grid of R, phi and Z for each location
    Rgrid, phigrid, zgrid= [], [], []
    for loc in locations:
        lcen, bcen= apo.glonGlat(loc)
        XYZ= bovy_coords.lbd_to_XYZ(lcen*numpy.ones_like(ds),
                                    bcen*numpy.ones_like(ds),
                                    ds,
                                    degree=True)
        Rphiz= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],XYZ[:,1],XYZ[:,2],
                                            Xsun=define_rgbsample._R0,
                                            Zsun=define_rgbsample._Z0)
        Rgrid.append(Rphiz[:,0])
        phigrid.append(Rphiz[:,1])
        zgrid.append(Rphiz[:,2])
    Rgrid= numpy.array(Rgrid)
    phigrid= numpy.array(phigrid)
    zgrid= numpy.array(zgrid)
    # Now compute rate(R) for each location and combine
    effsel*= numpy.tile(ds**2.*(distmods[1]-distmods[0])*(ds*np.log(10)/5.),(effsel.shape[0],1))
    tdens= densfunc(Rgrid,phigrid,zgrid)/dens_norm
    rate= tdens*effsel
    sumrate = np.sum(rate)
    #calculate normalisation N(R0)
    norm = (nbin/sumrate)
    #convert units (Kpc^2 > pc^2, deg > rad etc)
    norm *= 1e-6*(180/np.pi)**2
    #compute mass in bin using values from isochrones
    bin_mass = (norm*av_mass)/m_ratio
    if verbose==True:
        print bin_mass 
    return bin_mass, norm, m_ratio, (av_mass*1e-6*(180/np.pi)**2)/(sumrate*m_ratio)

def calculate_bin_error(samples, fehbin, agebin, nbin, iso_grid, 
						type='brokenexpflare',
						loggcut=[1.8,3.0],
						teffcut=[4000,5000], 
						n_sampling=1000, 
						progress=True, 
						mp=True, 
						fitIndx=None,
						weights = 'padova',
						distance_cut = False,
						lowermass = None):
    randsamples = np.random.permutation(samples.T)[:n_sampling]
    m_sample = np.zeros(np.shape(randsamples)[0])
    if multi == False:
        for ii,params in enumerate(randsamples):
            if progress==True:
                print ''+str(round(float(ii)/float(n_sampling)*100,2))+'% complete!'
            m = calc_normalisation(params, nbin , iso_grid, fehbin = fehbin, agebin=agebin, loggcut=loggcut, teffcut=teffcut, type=type, verbose=False, fitIndx=fitIndx, gridfile=gridfile, weights=weights, distance_cut = distance_cut, lowermass=lowermass)[0]
            m_sample[ii] = m
    if mp == True:
         m_sample= multi.parallel_map((lambda x: calc_normalisation(randsamples[x], nbin, iso_grid, fehbin=fehbin, agebin=agebin,loggcut=loggcut, teffcut=teffcut, type=type, verbose=False, fitIndx=fitIndx, distance_cut=distance_cut, lowermass=lowermass)[0]),\
         								range(np.shape(randsamples)[0]),numcores=numpy.amin([np.shape(randsamples)[0], multiprocessing.cpu_count()/2]))
    median = np.percentile(m_sample, 50)
    lowerr = np.percentile(m_sample, 16)
    uperr = np.percentile(m_sample, 84)
    return m_sample, median, lowerr, uperr 
    
