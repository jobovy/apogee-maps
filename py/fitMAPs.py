###############################################################################
# fitMAPs: fit the density distribution of MAPs
###############################################################################
import sys
import os, os.path
import pickle
import copy
import numpy
from galpy.util import save_pickles
import define_rcsample
import fitDens
_NSAMPLES= 50000
# Globals
locations= None
distmods= None
effsel= None
effsel_mar= None
effsel_drim= None
effsel_sale= None
effsel_zero= None
def fitMAPs(type,savefilename):
    setup_selection_function()
    if os.path.exists(savefilename):
        with open(savefilename,'rb') as savefile:
            bf= pickle.load(savefile)
            samples= pickle.load(savefile)
            bf_g15= pickle.load(savefile)
            samples_g15= pickle.load(savefile)
            bf_zero= pickle.load(savefile)
            samples_zero= pickle.load(savefile)
            bii= pickle.load(savefile)
    else:
        bf= []
        samples= []
        bf_g15= []
        samples_g15= []
        bf_zero= []
        samples_zero= []
        bii= 0
    maps= define_rcsample.MAPs()
    for ii, map in enumerate(maps.map()):
        if ii < bii: continue
        tbf, ts= fitmap(map,type=type,dmap='marshall06')
        bf.append(tbf)
        samples.append(ts)
        tbf, ts= fitmap(map,type=type,dmap='green15')
        bf_g15.append(tbf)
        samples_g15.append(ts)
        tbf, ts= fitmap(map,type=type,dmap='zero')
        bf_zero.append(tbf)
        samples_zero.append(ts)
        print ii, numpy.median(samples[-1],axis=1)
        save_pickles(savefilename,bf,samples,
                     bf_g15,samples_g15,
                     bf_zero,samples_zero,
                     ii+1)
    return None

def fitmap(tdata,type='brokenexp',dmap='marshall06'):
    tlocations= copy.deepcopy(locations)
    tdistmods= copy.deepcopy(distmods)
    if dmap == 'green15':
        teffsel= copy.deepcopy(effsel)
    elif dmap.lower() == 'marshall06':        
        teffsel= copy.deepcopy(effsel_mar)
    elif dmap.lower() == 'zero':        
        teffsel= copy.deepcopy(effsel_zero)
    bf, samples= fitDens.fitDens(tdata,
                                 numpy.array(tlocations),
                                 copy.deepcopy(teffsel),
                                 tdistmods,type=type,verbose=False,
                                 mcmc=True,nsamples=_NSAMPLES)
    bf= fitDens.fitDens(tdata,
                        numpy.array(tlocations),
                        copy.deepcopy(teffsel),
                        tdistmods,type=type,verbose=False,
                        init=numpy.median(samples,axis=1))
    return (bf, samples)

def setup_selection_function():
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
    return None

if __name__ == '__main__':
    fitMAPs(sys.argv[1],sys.argv[2])
