###############################################################################
# fitDens.py: module with routines to fit the density profile of a set of stars
###############################################################################
import os, os.path
import pickle
import numpy
from scipy import optimize
from galpy.util import bovy_coords
import bovy_mcmc
import acor
import densprofiles
import define_rcsample
def fitDens(data,
            locations,effsel,distmods,
            type='exp',
            mcmc=False,nsamples=10000,
            verbose=True,
            init=None,
            retMaxL=False,
            pos_keys = ['RC_GALR_H', 'RC_GALPHI_H', 'RC_GALZ_H']):
    """
    NAME:
       fitDens
    PURPOSE:
       fit the density profile of a set of stars
    INPUT:
       data - recarray with the data
       locations - locations of the APOGEE effective selection function
       effsel - array (nloc,nD) of the effective selection function, includes area of the field
       distmods - grid of distance moduli on which the effective selection function is pre-computed
       type= ('exp') type of density profile to fit
       mcmc= (False) run MCMC or not
       nsamples= (10000) number of MCMC samples to obtain
       verbose= (True) set this to False for no optimize convergence messages
       init= (None) if set, these are the initial conditions
       retMaxL= (False) if True, return the maximum likelihood
       pos_keys = the keys in the data table to R, phi and Z
    OUTPUT:
       (best-fit, samples, maximum-likelihood) based on options
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
       2017-03-04 - added pos_keys kwarg - Mackereth (ARI)
    """
    # Setup the density function and its initial parameters
    densfunc= _setup_densfunc(type)
    if init is None:
        init= _setup_initparams_densfunc(type,data, pos_keys = pos_keys)
    # Setup the integration of the effective volume
    effsel, Rgrid, phigrid, zgrid= _setup_effvol(locations,effsel,distmods)
    # Get the data's R,phi,z
    dataR= data[pos_keys[0]]
    dataphi= data[pos_keys[1]]
    dataz= data[pos_keys[2]]
    # Optimize
    out= optimize.fmin(lambda x: _mloglike(x,densfunc,type,
                                           dataR,dataphi,dataz,
                                           effsel,Rgrid,phigrid,zgrid),
                       init,disp=verbose)
    if 'explinflare' in type:
        step= [0.2,0.2,0.2,0.2,0.02] # need small step for Rfinv
    else:
        step= 0.2
    if mcmc:
        samples= bovy_mcmc.markovpy(out,
                                    step,
                                    lambda x: loglike(x,densfunc,type,
                                                      dataR,dataphi,dataz,
                                                      effsel,Rgrid,
                                                      phigrid,zgrid),
                                    (),
                                    isDomainFinite=[[False,False] for ii in range(len(out))],
                                    domain= [[0.,0.] for ii in range(len(out))],
                                    nsamples=nsamples,
                                    nwalkers=2*len(out))
        if verbose: print_samples_qa(samples)
        out= (out,numpy.array(samples).T,)
    else:
        out= (out,)
    if retMaxL:
        out= out+(loglike(out[0],densfunc,type,dataR,dataphi,dataz,
                          effsel,Rgrid,
                          phigrid,zgrid),)
    return out

def _mloglike(*args,**kwargs):
    """Minus the log likelihood"""
    return -loglike(*args,**kwargs)

def loglike(params,densfunc,type,
            dataR,dataphi,dataz,
            effsel,Rgrid,phigrid,zgrid):
    """
    NAME:
       loglike
    PURPOSE:
       compute the log likelihood of the data given a density profile
    INPUT:
       params - parameters of the density
       densfunc -  a function that evaluates the density as densfunc(R,phi,z,params=params)
       dataR, dataphi, dataz - the cylindrical Galactocentric coordinates of the data
       effsel - array (nloc,nD), includes D^2*Delta D factor
       Rgrid, phigrid, zgrid - array (nloc,nD) of the cylindrical Galactocentric coordinates corresponding to the (l_loc,b_loc,D) of the effective selection function
    OUTPUT:
       log likelihood
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    # Check priors
    if not _check_range_params_densfunc(params,type):
        return -numpy.finfo(numpy.dtype(numpy.float64)).max
    #raise NotImplementedError("Need to implement priors")
    # Setup the density function
    tdensfunc= lambda x,y,z: densfunc(x,y,z,params=params)
    # Evaluate the log density at the data
    datadens= numpy.log(tdensfunc(dataR,dataphi,dataz))
    # Evaluate the effective volume
    teffvol= effvol(tdensfunc,effsel,Rgrid,phigrid,zgrid)
    # Put it all together
    return numpy.sum(datadens)-len(dataR)*numpy.log(teffvol)

def effvol(densfunc,effsel,Rgrid,phigrid,zgrid):
    """
    NAME:
       effvol
    PURPOSE:
       calculate the effective volume
    INPUT:
       densfunc - function that returns the density when called with R,phi,z
       effsel - array (nloc,nD), includes D^2*Delta D factor
       Rgrid, phigrid, zgrid - array (nloc,nD) of the cylindrical Galactocentric coordinates corresponding to the (l_loc,b_loc,D) of the effective selection function
    OUTPUT:
       effective volume
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    # Evaluate the density
    tdens= densfunc(Rgrid,phigrid,zgrid)
    return numpy.sum(effsel*tdens)

############################## DENSITY FUNCTIONS #############################
def _setup_densfunc(type):
    """Return the density function for this type"""
    if type.lower() == 'exp':
        return densprofiles.expdisk
    elif type.lower() == 'expplusconst':
        return densprofiles.expdiskplusconst
    elif type.lower() == 'twoexp':
        return densprofiles.twoexpdisk
    elif type.lower() == 'brokenexp':
        return densprofiles.brokenexpdisk
    elif type.lower() == 'tribrokenexp':
        return densprofiles.tribrokenexpdisk
    elif type.lower() == 'symbrokenexp':
        return densprofiles.symbrokenexpdisk
    elif type.lower() == 'brokenexpflare':
        return densprofiles.brokenexpflaredisk
    elif type.lower() == 'tribrokenexpflare':
        return densprofiles.tribrokenexpflaredisk
    elif type.lower() == 'tribrokenexpfixedflare':
        return densprofiles.tribrokenexpfixedflaredisk
    elif type.lower() == 'brokentwoexp':
        return densprofiles.brokentwoexpdisk
    elif type.lower() == 'brokentwoexpflare':
        return densprofiles.brokentwoexpflaredisk
    elif type.lower() == 'tribrokentwoexp':
        return densprofiles.tribrokentwoexpdisk
    elif type.lower() == 'gaussexp':
        return densprofiles.gaussexpdisk
    elif type.lower() == 'brokenquadexp':
        return densprofiles.brokenquadexpdisk
    elif type.lower() == 'symbrokenquadexp':
        return densprofiles.symbrokenquadexpdisk
    elif type.lower() == 'brokenexpfixedspiral':
        return densprofiles.brokenexpdiskfixedspiral
    elif type.lower() == 'tribrokenexplinflare':
        return densprofiles.tribrokenexplinflaredisk
    elif type.lower() == 'tribrokenexpinvlinflare':
        return densprofiles.tribrokenexpinvlinflaredisk
def _setup_initparams_densfunc(type,data, pos_keys=['RC_GALR_H', 'RC_GALPHI_H', 'RC_GALZ_H']):
    """Return the initial parameters of the density for this type, might depend on the data"""
    if type.lower() == 'exp':
        return [1./3.,1./0.3]
    elif type.lower() == 'expplusconst':
        return [1./3.,1./0.3,numpy.log(0.1)]
    elif type.lower() == 'twoexp':
        return [1./3.,1./0.3,1./4.,1./0.5,densprofiles.logit(0.5)]
    elif type.lower() == 'brokenexp':
        return [-1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]]))]
    elif type.lower() == 'tribrokenexp':
        return [1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]]))]
    elif type.lower() == 'symbrokenexp':
        return [0.4,1./0.3,numpy.log(10.)]
    elif type.lower() == 'brokenexpflare':
        return [-1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                 -1./5.]
    elif type.lower() == 'tribrokenexpflare':
        return [1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                 -1./5.]
    elif type.lower() == 'tribrokenexpfixedflare':
        return [1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]]))]
    elif type.lower() == 'brokentwoexp':
        return [-1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                 densprofiles.logit(0.5),1./0.8]
    elif type.lower() == 'brokentwoexpflare':
        return [-1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                 densprofiles.logit(0.5),1./0.8,-0.2]
    elif type.lower() == 'tribrokentwoexp':
        return [1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                 densprofiles.logit(0.5),1./0.8]
    elif type.lower() == 'gaussexp':
        return [1./3.,1./0.3,numpy.log(10.)]
    elif type.lower() == 'brokenquadexp':
        return [1./3.,1./0.3,1./3.,numpy.log(10.)]
    elif type.lower() == 'symbrokenquadexp':
        return [1./3.,1./0.3,numpy.log(10.)]
    elif type.lower() == 'brokenexpfixedspiral':
        return [1./6.,1./0.3,1./2.,numpy.log(14.),numpy.log(1.)]
    elif type.lower() == 'tribrokenexplinflare':
        return [1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                0.]
    elif type.lower() == 'tribrokenexpinvlinflare':
        return [1./3.,1./0.3,1./3.,numpy.log(numpy.median(data[pos_keys[0]])),
                -1./5.]
def _check_range_params_densfunc(params,type):
    """Check that the current parameters are in a reasonable range (prior)"""
    # Second parameter is always a scale height, which we don't allow neg.
    if params[1] < 0.: return False
    if type.lower() == 'exp':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.exp(params[2]) > 1.: return False
    elif type.lower() == 'expplusconst':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.exp(params[2]) > 1.: return False
        if params[2] < -20.: return False
    elif type.lower() == 'twoexp':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.fabs(params[2]) > 2.: return False
        if numpy.fabs(params[3]) > 20.: return False
    elif type.lower() == 'brokenexp':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.fabs(params[2]) > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'tribrokenexp':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'symbrokenexp':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.exp(params[2]) > 16.: return False
        if numpy.exp(params[2]) < 1.: return False
    elif type.lower() == 'brokenexpflare':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.fabs(params[2]) > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'tribrokenexpflare':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'tribrokenexpfixedflare':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'brokentwoexp':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.fabs(params[2]) > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
        if params[4] < -7.: return False
        if params[4] > 0.: return False #make 2nd less dominant
        if params[5] < 0.: return False
        if numpy.fabs(params[5]) > 20.: return False
    elif type.lower() == 'brokentwoexpflare':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.fabs(params[2]) > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
        if params[4] < -7.: return False
        if params[4] > 0.: return False #make 2nd less dominant
        if params[5] < 0.: return False
        if numpy.fabs(params[5]) > 20.: return False
    elif type.lower() == 'tribrokentwoexp':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
        if params[4] < -7.: return False
        if params[4] > 0.: return False #make 2nd less dominant
        if params[5] < 0.: return False
        if params[5] > 20.: return False
    elif type.lower() == 'gaussexp':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.exp(params[2]) > 16.: return False
        if numpy.exp(params[2]) < 1.: return False
    elif type.lower() == 'brokenquadexp':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'symbrokenquadexp':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.exp(params[2]) > 16.: return False
        if numpy.exp(params[2]) < 1.: return False
    elif type.lower() == 'brokenexpfixedspiral':
        if numpy.fabs(params[0]) > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if numpy.fabs(params[2]) > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'tribrokenexplinflare':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    elif type.lower() == 'tribrokenexpinvlinflare':
        if params[0] < 0.: return False
        if params[0] > 2.: return False
        if numpy.fabs(params[1]) > 20.: return False
        if params[2] < 0.: return False
        if params[2] > 2.: return False
        if numpy.exp(params[3]) > 16.: return False
        if numpy.exp(params[3]) < 1.: return False
    return True

########################### EFFECTIVE VOLUME SETUP ############################
def _setup_effvol(locations,effsel,distmods):
    # First restore the APOGEE selection function (assumed pre-computed)
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    # Now compute the necessary coordinate transformations
    ds= 10.**(distmods/5-2.)
    Rgrid, phigrid, zgrid= [], [], []
    for loc in locations:
        lcen, bcen= apo.glonGlat(loc)
        XYZ= bovy_coords.lbd_to_XYZ(lcen*numpy.ones_like(ds),
                                    bcen*numpy.ones_like(ds),
                                    ds,
                                    degree=True)
        Rphiz= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],XYZ[:,1],XYZ[:,2],
                                            Xsun=define_rcsample._R0,
                                            Zsun=define_rcsample._Z0)
        Rgrid.append(Rphiz[:,0])
        phigrid.append(Rphiz[:,1])
        zgrid.append(Rphiz[:,2])
    Rgrid= numpy.array(Rgrid)
    phigrid= numpy.array(phigrid)
    zgrid= numpy.array(zgrid)
    # Also need to multiply in distance factors
    effsel*= numpy.tile(ds**3.*(distmods[1]-distmods[0]),(effsel.shape[0],1))
    return (effsel,Rgrid,phigrid,zgrid)

# From vclos before...
def print_samples_qa(samples):
    print "Mean, standard devs, acor tau, acor mean, acor s ..."
    for kk in range(len(samples[0])):
        xs= numpy.array([s[kk] for s in samples])
        #Auto-correlation time
        tau, m, s= acor.acor(xs)
        print numpy.mean(xs), numpy.std(xs), tau, m, s
