
import densprofiles
import define_rgbsample
import fitDens
import compareDataModel
import numpy
import os, os.path
import pickle
import copy
from galpy.util import bovy_plot, bovy_coords
import mwdust
import fitsio

selectFile= '../savs/selfunc-nospdata.sav'
if os.path.exists(selectFile):
    with open(selectFile,'rb') as savefile:
        apo= pickle.load(savefile)
with open('../essf/essf_green15_rgb.sav','rb') as savefile:
    locations= pickle.load(savefile)
    effsel= pickle.load(savefile)
    distmods= pickle.load(savefile)
with open('../essf/essf_marshall06_rgb.sav','rb') as savefile:
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
            
            
def sample(type='brokenexp',fitIndx=None,data_fitIndx=None,nsamples=3000,dmap='green15', ret=True):
    if fitIndx is None:
        fitIndx= numpy.ones(len(locations),dtype='bool') #True-betwDiskIndx
        data_fitIndx= numpy.ones(len(ldata),dtype='bool') #True-data_betwDiskIndx
    if dmap == 'green15':
        tlocations= copy.deepcopy(locations)
        teffsel= copy.deepcopy(effsel)
        tdistmods= copy.deepcopy(distmods)
    elif dmap.lower() == 'marshall06':        
        tlocations= copy.deepcopy(mlocations)
        teffsel= copy.deepcopy(meffsel)
        tdistmods= copy.deepcopy(mdistmods)
    elif dmap.lower() == 'zero':        
        tlocations= copy.deepcopy(zlocations)
        teffsel= copy.deepcopy(zeffsel)
        tdistmods= copy.deepcopy(zdistmods)
    elif dmap.lower() == 'drimmel03':        
        tlocations= copy.deepcopy(dlocations)
        teffsel= copy.deepcopy(deffsel)
        tdistmods= copy.deepcopy(ddistmods)
    bf, samples= fitDens.fitDens(ldata[data_fitIndx],numpy.array(tlocations)[fitIndx],copy.deepcopy(teffsel)[fitIndx],
                                 tdistmods,type=type,
                                 nsamples=nsamples,mcmc=True)
    labels= []
    for ii in range(len(bf)): labels.append(r"$\mathrm{param}\ %i$" % ii)
    corner.corner(samples.T,quantiles=[0.16, 0.5, 0.84],labels=labels,
                         show_titles=True, title_args={"fontsize": 12})
    if ret == True:
        return bf, samples.T
    else:
        return None

def ret_sample(type='brokenexp',fitIndx=None,data_fitIndx=None,nsamples=3000,dmap='green15'):
    if fitIndx is None:
        fitIndx= numpy.ones(len(locations),dtype='bool') #True-betwDiskIndx
        data_fitIndx= numpy.ones(len(ldata),dtype='bool') #True-data_betwDiskIndx
    if dmap == 'green15':
        tlocations= copy.deepcopy(locations)
        teffsel= copy.deepcopy(effsel)
        tdistmods= copy.deepcopy(distmods)
    elif dmap.lower() == 'marshall06':        
        tlocations= copy.deepcopy(mlocations)
        teffsel= copy.deepcopy(meffsel)
        tdistmods= copy.deepcopy(mdistmods)
    elif dmap.lower() == 'zero':        
        tlocations= copy.deepcopy(zlocations)
        teffsel= copy.deepcopy(zeffsel)
        tdistmods= copy.deepcopy(zdistmods)
    elif dmap.lower() == 'drimmel03':        
        tlocations= copy.deepcopy(dlocations)
        teffsel= copy.deepcopy(deffsel)
        tdistmods= copy.deepcopy(ddistmods)
    bf, samples= fitDens.fitDens(ldata[data_fitIndx],numpy.array(tlocations)[fitIndx],copy.deepcopy(teffsel)[fitIndx],
                                 tdistmods,type=type,
                                 nsamples=nsamples,mcmc=True)
    
    return bf, samples

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
    
agebins = numpy.arange(1.,14.,2.)
fehbins = numpy.arange(-0.6,0.3,0.1)
paramt = []
samples = []
numbins = []
load_data(subsample='fehage', add_ages=True, agebin=[0., 13.], fehbin=[-0.6,0.2], afebin='highnew', agetype='Martig', corrections=True)
dat = ldata
for j in range(0,len(fehbins)-1):
    paramtfeh = []
    samplesfeh = []
    numbinfeh = []
    
    for i in range(0,len(agebins)-1):
        loadeffsel_maps(fehbin=fehbins[j], agebin=agebins[i])
        mask = (dat['Age'] >= agebins[i])&(dat['Age'] < agebins[i+1])&(dat['FE_H'] >= fehbins[j])&(dat['FE_H'] < fehbins[j+1])
        
        ldata = dat[mask]
        print len(ldata)
        fitIndx= None #True-highbIndx
        data_fitIndx= None #True-data_highbIndx
        bf, sample = ret_sample(type='brokenexpflare',fitIndx=fitIndx,data_fitIndx=data_fitIndx,nsamples=10000,dmap='marshall06')
        num_bin = len(dat[mask])
        paramtfeh.append(bf)
        samplesfeh.append(sample)
        numbinfeh.append(num_bin)
    paramt.append(paramtfeh)
    samples.append(samplesfeh)
    numbins.append(numbinfeh)

samples = np.array(samples)
obj = [agebins, fehbins, numbins, paramt, samples]
savfile = open('../savs/paramsRGB_brokenexpflare_01dex2gyrbins_high.dat', 'w')
pickle.dump(obj, savfile)

agebins = numpy.arange(1.,14.,2.)
fehbins = numpy.arange(-0.6,0.3,0.1)
paramt = []
samples = []
numbins = []
load_data(subsample='fehage', add_ages=True, agebin=[0., 13.], fehbin=[-0.6,0.2], afebin='lownew', agetype='Martig', corrections=True)
dat = ldata
for j in range(0,len(fehbins)-1):
    paramtfeh = []
    samplesfeh = []
    numbinfeh = []
    
    for i in range(0,len(agebins)-1):
        loadeffsel_maps(fehbin=fehbins[j], agebin=agebins[i])
        mask = (dat['Age'] >= agebins[i])&(dat['Age'] < agebins[i+1])&(dat['FE_H'] >= fehbins[j])&(dat['FE_H'] < fehbins[j+1])
        
        ldata = dat[mask]
        print len(ldata)
        fitIndx= None #True-highbIndx
        data_fitIndx= None #True-data_highbIndx
        bf, sample = ret_sample(type='brokenexpflare',fitIndx=fitIndx,data_fitIndx=data_fitIndx,nsamples=10000,dmap='marshall06')
        num_bin = len(dat[mask])
        paramtfeh.append(bf)
        samplesfeh.append(sample)
        numbinfeh.append(num_bin)
    paramt.append(paramtfeh)
    samples.append(samplesfeh)
    numbins.append(numbinfeh)

samples = np.array(samples)
obj = [agebins, fehbins, numbins, paramt, samples]
savfile = open('../savs/paramsRGB_brokenexpflare_01dex2gyrbins_low.dat', 'w')
pickle.dump(obj, savfile)