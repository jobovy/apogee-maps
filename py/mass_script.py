import pickle
import numpy as np
import numpy
import sys
import multiprocessing
from tqdm import tqdm
import define_rgbsample
from galpy.util import multi
import os
try:
    reload(calc_masses)
except NameError:
    import calc_masses
    
selectFile= '../savs/selfunc-nospdata.sav'
if os.path.exists(selectFile):
    with open(selectFile,'rb') as savefile:
        apo= pickle.load(savefile)
with open('../essf/maps/essf_rgb_green15_modelmh_feh-0.0_age1.0.sav','rb') as savefile:
    locations= pickle.load(savefile)
    effsel= pickle.load(savefile)
    distmods= pickle.load(savefile)
with open('../essf/maps/essf_rgb_marshall06_modelmh_feh-0.0_age1.0.sav','rb') as savefile:
    mlocations= pickle.load(savefile)
    meffsel= pickle.load(savefile)
    mdistmods= pickle.load(savefile)
# Fill in regions not covered by Marshall map
meffsel[meffsel < -0.5]= effsel[meffsel < -0.5]

# Get (lcen,bcen) for each location
lcen= numpy.zeros(len(locations))
bcen= numpy.zeros(len(locations))
hmax= numpy.zeros(len(locations))
for ii,loc in enumerate(locations):
    if loc in apo.list_fields():
        tlcen, tbcen= apo.glonGlat(loc)
        lcen[ii]= tlcen
        bcen[ii]= tbcen
        hmax[ii]= apo.Hmax(loc,cohort='long')
        if numpy.isnan(hmax[ii]):
            hmax[ii]= apo.Hmax(loc,cohort='medium')
            if numpy.isnan(hmax[ii]):
                hmax[ii]= apo.Hmax(loc,cohort='short')
    if loc not in apo.list_fields():
        lcen[ii] = numpy.nan
        bcen[ii] = numpy.nan
        hmax[ii]= numpy.nan
        
# Get the locations of various subsamples
medbIndx = numpy.fabs(bcen) > 6.
highbIndx= numpy.fabs(bcen) > 10.
outDiskIndx= (lcen > 150.)*(lcen < 250.)*(True-highbIndx)
betwDiskIndx= (lcen <= 150.)*(lcen >= 70.)*(True-highbIndx)
inDiskIndx= (lcen < 70.)*(lcen >= 25.)*(True-highbIndx)
bulgeIndx= ((lcen < 25.)+(lcen > 335.))*(True-highbIndx)
brightIndx= (hmax <= 12.21)
mediumIndx= (hmax > 12.21)*(hmax <= 12.81)
faintIndx= (hmax > 12.81)

ldata= None
data_highbIndx= None
data_outDiskIndx= None
data_betwDiskIndx= None
data_inDiskIndx= None
data_bulgeIndx= None
data_brightIndx= None
data_mediumIndx= None
data_faintIndx= None

def load_data(subsample='lowlow', add_ages=False, agebin=[0.,1.], fehbin=[0.,0.2], afebin=None, agetype='Martig', corrections=False):
    global ldata
    global data_highbIndx
    global data_outDiskIndx
    global data_betwDiskIndx
    global data_inDiskIndx
    global data_bulgeIndx
    global data_brightIndx
    global data_mediumIndx
    global data_faintIndx
    if subsample.lower() == 'all':
        ldata= define_rgbsample.get_rgbsample(agetype=agetype)
    elif subsample.lower() == 'alllowalpha':
    	ldata= define_rgbsample.get_rgbsample()
    	ldata= ldata[ldata[define_rgbsample._AFETAG] < 0.1]
    elif subsample.lower() == 'lowlow':
        ldata= define_rgbsample.get_lowlowsample()
    elif subsample.lower() == 'highfeh':
        ldata= define_rgbsample.get_highfehsample()
    elif subsample.lower() == 'highalpha':
        ldata= define_rgbsample.get_highalphasample()
    elif subsample.lower() == 'solar':
        ldata= define_rgbsample.get_solarsample()
    elif subsample.lower() == 'fehage':
        ldata= define_rgbsample.get_fehage(agebin=agebin, fehbin=fehbin, afebin=afebin, agetype=agetype, apply_corrections=corrections)
    # Get the indices of the various subsamples defined above
    data_highbIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[highbIndx]: data_highbIndx[ii]= True
    data_outDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[outDiskIndx]: data_outDiskIndx[ii]= True
    data_betwDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[betwDiskIndx]: data_betwDiskIndx[ii]= True
    data_inDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[inDiskIndx]: data_inDiskIndx[ii]= True
    data_bulgeIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[bulgeIndx]: data_bulgeIndx[ii]= True
    data_brightIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[brightIndx]: data_brightIndx[ii]= True
    data_mediumIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[mediumIndx]: data_mediumIndx[ii]= True
    data_faintIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[faintIndx]: data_faintIndx[ii]= True            


def loadeffsel_maps(sample='rgb', fehbin=-0.6, agebin=1.0):
    global locations, effsel, distmods, mlocations, meffsel, mdistmods
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    with open('../essf/maps/essf_'+sample+'_green15_modelmh_feh'+str(round(fehbin,1))+'_age'+str(round(agebin,1))+'.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel= pickle.load(savefile)
        distmods= pickle.load(savefile)
    with open('../essf/maps/essf_'+sample+'_marshall06_modelmh_feh'+str(round(fehbin,1))+'_age'+str(round(agebin,1))+'.sav','rb') as savefile:
        mlocations= pickle.load(savefile)
        meffsel= pickle.load(savefile)
        mdistmods= pickle.load(savefile)
    # Fill in regions not covered by Marshall map
    meffsel[meffsel < -0.5]= effsel[meffsel < -0.5]

    # Get (lcen,bcen) for each location
    lcen= numpy.zeros(len(locations))
    bcen= numpy.zeros(len(locations))
    hmax= numpy.zeros(len(locations))
    for ii,loc in enumerate(locations):
        if loc in apo.list_fields():
            tlcen, tbcen= apo.glonGlat(loc)
            lcen[ii]= tlcen
            bcen[ii]= tbcen
            hmax[ii]= apo.Hmax(loc,cohort='long')
            if numpy.isnan(hmax[ii]):
                hmax[ii]= apo.Hmax(loc,cohort='medium')
                if numpy.isnan(hmax[ii]):
                    hmax[ii]= apo.Hmax(loc,cohort='short')
        if loc not in apo.list_fields():
            lcen[ii] = numpy.nan
            bcen[ii] = numpy.nan
            hmax[ii]= numpy.nan

    # Get the locations of various subsamples
    highbIndx= numpy.fabs(bcen) > 10.
    outDiskIndx= (lcen > 150.)*(lcen < 250.)*(True-highbIndx)
    betwDiskIndx= (lcen <= 150.)*(lcen >= 70.)*(True-highbIndx)
    inDiskIndx= (lcen < 70.)*(lcen >= 25.)*(True-highbIndx)
    bulgeIndx= ((lcen < 25.)+(lcen > 335.))*(True-highbIndx)
    brightIndx= (hmax <= 12.21)
    mediumIndx= (hmax > 12.21)*(hmax <= 12.81)
    faintIndx= (hmax > 12.81)
    
#Load from saved file
savfile = open('../out/paramsRGB_brokenexpflare_01dex2gyrbins_low.dat', 'rb')obj = pickle.load(savfile)
agebins, fehbins, numbins, paramt, samples = obj
samples = np.array(samples)
    
#make a grid of the best fit values (median of the mcmc sampling)
agebincent = (agebins[:-1]+agebins[1:])/2.
fehbincent = (fehbins[:-1]+fehbins[1:])/2.

agebins = [1.,3.,5.,7.,9.,11.,13.0]
fehbins = np.arange(-0.6,0.3,0.1)
numbins = []

load_data(subsample='fehage', add_ages=True, agebin=[0., 13.5], fehbin=[-0.6,0.2], afebin='low', agetype='Martig', corrections=True)
dat = ldata#[data_medbIndx]
for j in range(0,len(fehbins)-1):
    numbinfeh = []
    for i in range(0,len(agebins)-1):
        mask = (dat['Age'] >= agebins[i])&(dat['Age'] < agebins[i+1])&(dat['FE_H'] >= fehbins[j])&(dat['FE_H'] < fehbins[j+1])
        ldata = dat[mask]
        num_bin = len(ldata)+(int(0.25*len(ldata)))    ##### ADDED MISSSING STARS PLEASE NOTE
        numbinfeh.append(num_bin)
    numbins.append(numbinfeh)

grid = []
dgrid = []
for j in range(0,len(fehbincent)):
    vals = []
    dvals = []
    for i in range(0,len(agebincent)):
        hRin = np.median(samples[j,i,0,::10])
        dhRin = np.std(samples[j,i,0,::10])
        hz = np.median(samples[j,i,1,::10])
        dhz = np.std(samples[j,i,1,::10])
        hRout = np.median(samples[j,i,2,::10])
        dhRout = np.std(samples[j,i,2,::10])
        rb = np.median(samples[j,i,3,::10])
        drb = np.std(samples[j,i,3,::10])
        rf = np.median(samples[j,i,4,::10])
        drf = np.std(samples[j,i,4,::10])
        nsamples = 1000
        low = np.sort(1/samples[j,i,1,::10])[int(round(0.025*nsamples))]
        high = np.sort(1/samples[j,i,1,::10])[int(round((1.-0.025)*nsamples))]
        val = [hRin, hz, hRout, rb, rf]
        dval = [dhRin, dhz, dhRout, drb, drf]
        vals.append(val)
        dvals.append(dval)
    grid.append(vals)
    dgrid.append(dvals)
grid= np.array(grid)
dgrid = np.array(dgrid)

type = 'brokenexpflare'

tagebins = [1.,3.,5.,7.,9.,11.,13.0]
agebins = [0.,3.,5.,7.,9.,11.,13.0]
fehbins = np.arange(-0.6,0.3,0.1)
    
massgrid = np.zeros((len(fehbincent),len(agebincent), 5))
m_samplegrid = np.zeros((len(fehbincent),len(agebincent), 1000))
print 'Calculating surface-mass density for low alpha populations...'
for j in tqdm(range(0,len(fehbincent))):
    for i in range(0,len(agebincent)):
        #print 'calculating mass for feh '+str(fehbins[j])+' age '+str(agebins[i])
        iso_grid = calc_masses.load_isochrones('../savs/Padova_grid_lognormalchabrier2001.sav')
        sys.stdout.flush()
        loadeffsel_maps(sample='rgb', fehbin=fehbins[j], agebin=tagebins[i])
        m = calc_masses.calc_normalisation(grid[j,i], np.array(numbins)[j,i], iso_grid, fehbin = [fehbins[j], fehbins[j+1]], agebin=[agebins[i], agebins[i+1]],loggcut=[1.8,3.0], teffcut=[4100,5100], type=type, fitIndx=None, weights='padova', distance_cut=False, lowermass=None)
        #print 'sampling for error...'
        sys.stdout.flush()
        m_sample, m_med, m_low, m_up = calc_masses.calculate_bin_error(samples[j,i], [fehbins[j], fehbins[j+1]], [agebins[i], agebins[i+1]], np.array(numbins)[j,i], iso_grid, loggcut=[1.8,3.0], teffcut=[4100,5100], progress=False, n_sampling=1000, fitIndx=None, weights='padova', distance_cut=False, lowermass=None, type=type)
        nbin = np.array(numbins)[j,i]
        m_low_poisson = m[3]*(nbin-np.sqrt(nbin))-(m[0]-m_low)
        m_up_poisson = m[3]*(nbin+np.sqrt(nbin))+(m_up-m[0])
        print m_med, m_med - m_low, m_up - m_med , m_low_poisson, m_up_poisson
        massgrid[j,i] = [m[0], m_low, m_up, m_low_poisson, m_up_poisson]
        m_samplegrid[j,i] = m_sample
        sys.stdout.flush()
        
samples = np.array(samples)
obj = [agebins, fehbins, numbins, paramt, samples, massgrid, m_samplegrid]
savfile = open('../out/paramsRGB_brokenexpflare_01dex2gyrbins_low_mass.dat', 'w')
pickle.dump(obj, savfile)
savfile.close()

    
#Load from saved file
savfile = open('../out/paramsRGB_brokenexpflare_01dex2gyrbins_high.dat', 'rb')obj = pickle.load(savfile)
agebins, fehbins, numbins, paramt, samples = obj
samples = np.array(samples)
    
#make a grid of the best fit values (median of the mcmc sampling)
agebincent = (agebins[:-1]+agebins[1:])/2.
fehbincent = (fehbins[:-1]+fehbins[1:])/2.


agebins = [1.,3.,5.,7.,9.,11.,13.0]
fehbins = np.arange(-0.6,0.3,0.1)
numbins = []
load_data(subsample='fehage', add_ages=True, agebin=[0., 13.5], fehbin=[-0.6,0.2], afebin='high', agetype='Martig', corrections=True)
dat = ldata#[data_medbIndx]
for j in range(0,len(fehbins)-1):
    numbinfeh = []
    for i in range(0,len(agebins)-1):
        mask = (dat['Age'] >= agebins[i])&(dat['Age'] < agebins[i+1])&(dat['FE_H'] >= fehbins[j])&(dat['FE_H'] < fehbins[j+1])
        ldata = dat[mask]
        num_bin = len(ldata)+(int(0.25*len(ldata)))    ##### ADDED MISSSING STARS PLEASE NOTE
        numbinfeh.append(num_bin)
    numbins.append(numbinfeh)

grid = []
dgrid = []
for j in range(0,len(fehbincent)):
    vals = []
    dvals = []
    for i in range(0,len(agebincent)):
        hRin = np.median(samples[j,i,0,::10])
        dhRin = np.std(samples[j,i,0,::10])
        hz = np.median(samples[j,i,1,::10])
        dhz = np.std(samples[j,i,1,::10])
        hRout = np.median(samples[j,i,2,::10])
        dhRout = np.std(samples[j,i,2,::10])
        rb = np.median(samples[j,i,3,::10])
        drb = np.std(samples[j,i,3,::10])
        rf = np.median(samples[j,i,4,::10])
        drf = np.std(samples[j,i,4,::10])
        nsamples = 1000
        low = np.sort(1/samples[j,i,1,::10])[int(round(0.025*nsamples))]
        high = np.sort(1/samples[j,i,1,::10])[int(round((1.-0.025)*nsamples))]
        val = [hRin, hz, hRout, rb, rf]
        dval = [dhRin, dhz, dhRout, drb, drf]
        vals.append(val)
        dvals.append(dval)
    grid.append(vals)
    dgrid.append(dvals)
grid= np.array(grid)
dgrid = np.array(dgrid)

type = 'brokenexpflare'
tagebins = [1.,3.,5.,7.,9.,11.,13.0]
agebins = [0.,3.,5.,7.,9.,11.,13.0]
fehbins = np.arange(-0.6,0.3,0.1)

print 'calculating surface-mass densities for high alpha populations...'
massgrid = np.zeros((len(fehbincent),len(agebincent), 5))
m_samplegrid = np.zeros((len(fehbincent),len(agebincent), 1000))
for j in tqdm(range(0,len(fehbincent))):
    for i in range(0,len(agebincent)):
        #print 'calculating mass for feh '+str(fehbins[j])+' age '+str(agebins[i])
        iso_grid = calc_masses.load_isochrones('../savs/Padova_grid_lognormalchabrier2001.sav')
        sys.stdout.flush()
        loadeffsel_maps(sample='rgb', fehbin=fehbins[j], agebin=tagebins[i])
        m = calc_masses.calc_normalisation(grid[j,i], np.array(numbins)[j,i], iso_grid, fehbin = [fehbins[j], fehbins[j+1]], agebin=[agebins[i], agebins[i+1]],loggcut=[1.8,3.0], teffcut=[4100,5100], type=type, fitIndx=None, weights='padova', distance_cut=False, lowermass=None)
        #print 'sampling for error...'
        sys.stdout.flush()
        m_sample, m_med, m_low, m_up = calc_masses.calculate_bin_error(samples[j,i], [fehbins[j], fehbins[j+1]], [agebins[i], agebins[i+1]], np.array(numbins)[j,i], iso_grid, loggcut=[1.8,3.0], teffcut=[4100,5100], progress=False, n_sampling=1000, fitIndx=None, weights='padova', distance_cut=False, lowermass=None, type=type)
        nbin = np.array(numbins)[j,i]
        m_low_poisson = m[3]*(nbin-np.sqrt(nbin))-(m[0]-m_low)
        m_up_poisson = m[3]*(nbin+np.sqrt(nbin))+(m_up-m[0])
        print m_med, m_med - m_low, m_up - m_med , m_low_poisson, m_up_poisson
        massgrid[j,i] = [m[0], m_low, m_up, m_low_poisson, m_up_poisson]
        m_samplegrid[j,i] = m_sample
        sys.stdout.flush()
        
samples = np.array(samples)
obj = [agebins, fehbins, numbins, paramt, samples, massgrid, m_samplegrid]
savfile = open('../out/paramsRGB_brokenexpflare_01dex2gyrbins_high_mass.dat', 'w')
pickle.dump(obj, savfile)
savfile.close()
