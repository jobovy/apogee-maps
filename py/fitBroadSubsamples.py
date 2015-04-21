###############################################################################
# fitBroadSubsamples.py: do and report fits for the broad sub-samples
###############################################################################
import sys
import os, os.path
import pickle
import copy
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot, bovy_coords, save_pickles
import triangle
import mwdust
import densprofiles
import define_rcsample
import fitDens
import compareDataModel
_NSAMPLES= 3000
_LW= 1.5
# Globals
locations= None
distmods= None
effsel= None
effsel_mar= None
effsel_drim= None
effsel_sale= None
effsel_zero= None
lcen= None
bcen= None
highbIndx= None
outDiskIndx= None
betwDiskIndx= None 
inDiskIndx= None
ldata= None
data_highbIndx= None
data_outDiskIndx= None
data_betwDiskIndx= None
data_inDiskIndx= None
data_bulgeIndx= None
data_brightIndx= None
data_mediumIndx= None
data_faintIndx= None
def fitBroadSubsamples(sample,savename):
    # Setup the selection function
    setup_selection_function()
    # Load the data
    load_data(sample)
    # Do the fits
    global bf_exp, bf_brexp, bf_twoexp
    global ml_exp, ml_brexp, ml_twoexp
    global samples_exp, samples_brexp, samples_twoexp
    if os.path.exists(savename):
        with open(savename,'rb') as savefile:
            bf_exp= pickle.load(savefile)
            bf_brexp= pickle.load(savefile)
            bf_twoexp= pickle.load(savefile)
            ml_exp= pickle.load(savefile)
            ml_brexp= pickle.load(savefile)
            ml_twoexp= pickle.load(savefile)
            samples_exp= pickle.load(savefile)
            samples_brexp= pickle.load(savefile)
            samples_twoexp= pickle.load(savefile)
    else:
        # Perform fits of
        # a) expplusconst
        bf_exp, samples_exp, ml_exp= fit(type='expplusconst')
        # b) brokenexpflare
        bf_brexp, samples_brexp, ml_brexp= fit(type='brokenexpflare')
        # c) brokentwoexp
        bf_twoexp, samples_twoexp, ml_twoexp= fit(type='brokentwoexp')
        save_pickles(savename,
                     bf_exp,bf_brexp,bf_twoexp,
                     ml_exp,ml_brexp,ml_twoexp,
                     samples_exp,samples_brexp,samples_twoexp)
    # Do the rest of the fits as justfit
    global bf_brexp_g15, ml_brexp_g15, bf_brexp_drim, ml_brexp_drim
    global bf_brexp_sale, ml_brexp_sale, bf_brexp_zero, ml_brexp_zero
    bf_brexp_g15, ml_brexp_g15= fit(type='brokenexpflare',dmap='green15',
                                    justfit=True,init=bf_brexp) 
    bf_brexp_drim, ml_brexp_drim= fit(type='brokenexpflare',dmap='drimmel03',
                                      justfit=True,init=bf_brexp)
    bf_brexp_sale, ml_brexp_sale= fit(type='brokenexpflare',dmap='sale14',
                                      justfit=True,init=bf_brexp)
    bf_brexp_zero, ml_brexp_zero= fit(type='brokenexpflare',dmap='zero',
                                      justfit=True,init=bf_brexp)
    return None

def fit(type='brokenexpflare',dmap='marshall06',init=None,justfit=False):
    # Get the correct selection function
    tlocations= copy.deepcopy(locations)
    tdistmods= copy.deepcopy(distmods)
    if dmap == 'green15':
        teffsel= copy.deepcopy(effsel)
    elif dmap.lower() == 'marshall06':        
        teffsel= copy.deepcopy(effsel_mar)
    elif dmap.lower() == 'sale14':        
        teffsel= copy.deepcopy(effsel_sale)
    elif dmap.lower() == 'drimmel03':        
        teffsel= copy.deepcopy(effsel_drim)
    elif dmap.lower() == 'zero':        
        teffsel= copy.deepcopy(effsel_zero)
    # Now fit, sample, fit
    fitOut= fitDens.fitDens(ldata,
                            numpy.array(tlocations),
                            copy.deepcopy(teffsel),
                            tdistmods,
                            type=type,
                            init=init,
                            retMaxL=True,
                            nsamples=_NSAMPLES,mcmc=not justfit)
    if justfit: return fitOut
    bf, samples, maxl= fitOut
    bf, maxl= fitDens.fitDens(ldata,
                              numpy.array(tlocations),
                              copy.deepcopy(teffsel),
                              tdistmods,
                              type=type,
                              init=numpy.median(samples,axis=1),
                              retMaxL=True)
    return (bf,samples,maxl)

def setup_selection_function():
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    # Green et al. 
    global locations, distmods, effsel, effsel_mar, effsel_drim, effsel_sale
    global effsel_zero
    with open('../essf/essf_green15.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel= pickle.load(savefile)
        distmods= pickle.load(savefile)
    # Marshall et al. (2006)
    with open('../essf/essf_marshall06.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_mar= pickle.load(savefile)
    # Fill in regions not covered by Sale map
    effsel_mar[effsel_mar < -0.5]= effsel[effsel_mar < -0.5]
    # Sale et al. (2014)
    with open('../essf/essf_sale14.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_sale= pickle.load(savefile)
    # Fill in regions not covered by Marshall map
    effsel_sale[effsel_sale < -0.5]= effsel[effsel_sale < -0.5]
    # Drimmel et al (2003)
    with open('../essf/essf_drimmel03.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_drim= pickle.load(savefile)
    # Zero
    with open('../essf/essf_zero.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_zero= pickle.load(savefile)
    # Get (lcen,bcen) for each location
    global lcen, bcen
    lcen= numpy.zeros(len(locations))
    bcen= numpy.zeros(len(locations))
    for ii,loc in enumerate(locations):
        tlcen, tbcen= apo.glonGlat(loc)
        lcen[ii]= tlcen
        bcen[ii]= tbcen
    # Get the locations of various subsamples
    global highbIndx, outDiskIndx, betwDiskIndx, inDiskIndx
    highbIndx= numpy.fabs(bcen) > 10.
    outDiskIndx= (lcen > 140.)*(lcen < 250.)*(True-highbIndx)
    betwDiskIndx= (lcen <= 140.)*(lcen >= 70.)*(True-highbIndx)
    inDiskIndx= (lcen < 70.)*(True-highbIndx)
    return None

def load_data(sample):
    global ldata
    global data_highbIndx
    global data_outDiskIndx
    global data_betwDiskIndx
    global data_inDiskIndx    
    global data_bulgeIndx    
    global data_brightIndx
    global data_mediumIndx
    global data_faintIndx
    if sample.lower() == 'all':
        ldata= define_rcsample.get_rcsample()
    elif sample.lower() == 'alllowalpha':
        ldata= define_rcsample.get_rcsample()
        ldata= ldata[ldata[define_rcsample._AFETAG] < 0.1]
    elif sample.lower() == 'lowlow':
        ldata= define_rcsample.get_lowlowsample()
    elif sample.lower() == 'highfeh':
        ldata= define_rcsample.get_highfehsample()
    elif sample.lower() == 'highalpha':
        ldata= define_rcsample.get_highalphasample()
    elif sample.lower() == 'solar':
        ldata= define_rcsample.get_solarsample()
    # Get the indices of the various subsamples defined above
    data_highbIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[highbIndx]:
            data_highbIndx[ii]= True
    data_outDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[outDiskIndx]:
            data_outDiskIndx[ii]= True
    data_betwDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[betwDiskIndx]:
            data_betwDiskIndx[ii]= True
    data_inDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[inDiskIndx]:
            data_inDiskIndx[ii]= True
    return None

def plotCompareData(sample,savename,plotname):
    locs= ['highb','outdisk','meddisk','indisk']
    indices= [highbIndx,outDiskIndx,betwDiskIndx,inDiskIndx]
    data_indices= [data_highbIndx,data_outDiskIndx,
                   data_betwDiskIndx,data_inDiskIndx]
    # Full prediction for numbers
    Xs,pd= compareDataModel.predict_spacedist(bf_brexp,
                                              numpy.array(locations),
                                              copy.deepcopy(effsel_mar),
                                              distmods,
                                              type='brokenexpflare',coord='dm')
    for loc, index, data_index in zip(locs,indices,data_indices):
        bovy_plot.bovy_print()
        # High |b|
        # Marshall is fiducial
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_mar)[index],
                                                   distmods,type='brokenexpflare',
                                                   coord='dm')
        bovy_plot.bovy_hist(ldata['RC_DM_H'][data_index],
                            histtype='step',normed=True,
                            lw=2.*_LW,
                            range=[7.,15.5],
                            bins=round(numpy.sqrt(numpy.sum(data_index))),
                            yrange=[0.,
                                    1.2*numpy.amax(pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]))],
                            color='k',
                            xlabel=r'$\mu$')
        bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),color='r',
                            lw=2.*_LW,overplot=True,zorder=12)
        bovy_plot.bovy_text(r'$%i / %i\ \mathrm{stars}$' \
                                % (numpy.sum(data_index),
                                   len(ldata))
                            +'\n'+r'$%i / %i\ \mathrm{predicted}$' \
                                % (numpy.sum(pdt)/numpy.sum(pd)*len(ldata),len(ldata)),
                            top_left=True,size=14.)
        # Green
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_g15,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel)[index],
                                                   distmods,type='brokenexpflare',
                                                   coord='dm')
        bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),color='b',
                            lw=_LW,overplot=True,zorder=13)
        # Drimmel
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_drim,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_drim)[index],
                                                   distmods,type='brokenexpflare',
                                                   coord='dm')
        bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),color='gold',
                            lw=_LW,overplot=True,zorder=12)
        # Sale
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_sale,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_sale)[index],
                                                   distmods,type='brokenexpflare',
                                                   coord='dm')
        bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),color='c',
                            lw=_LW,overplot=True,zorder=12)
        # Zero
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_zero,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_zero)[index],
                                                   distmods,type='brokenexpflare',
                                                   coord='dm')
        bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),color='k',
                            ls='--',lw=_LW,overplot=True,zorder=10)
        if sample.lower() == 'lowlow' or sample.lower() == 'highalpha':
            if loc.lower() == 'highb':
                bovy_plot.bovy_text(r'$|b| > 10^\circ$',title=True,size=14.)
            elif loc.lower() == 'indisk':
                bovy_plot.bovy_text(r'$l < 70^\circ, |b| \leq 10^\circ$',
                                    title=True,size=14.)
            elif loc.lower() == 'meddisk':
                bovy_plot.bovy_text(r'$70^\circ \leq l \leq 140^\circ, |b| \leq 10^\circ$',
                                    title=True,size=14.)
            elif loc.lower() == 'outdisk':
                bovy_plot.bovy_text(r'$140^\circ < l < 250^\circ, |b| \leq 10^\circ$',
                                    title=True,size=14.)
        bovy_plot.bovy_end_print(plotname.replace('LOC',loc))
    return None

if __name__ == '__main__':
    # Input:
    #  - sample: 'lowlow', 'solar', 'highfeh', 'highalpha'
    #  - savename: name of the file for the pickle
    #  - tablename: name of the file for the table
    #  - plotname: name of the file for the plot
    sample= sys.argv[1]
    savename= sys.argv[2]
    tablename= sys.argv[3]
    plotname= sys.argv[4]
    # First, do the fits
    fitBroadSubsamples(sample,savename)
    # Then write the table
    #writeTable(sample,savename,tablename)
    # And make the plot comparing data and model
    plotCompareData(sample,savename,plotname)
