###############################################################################
# fitDens.py: module with routines to fit the density profile of a set of stars
###############################################################################
import os, os.path
import pickle
import numpy
from scipy import optimize
from galpy.util import bovy_coords
import densprofiles
import define_rcsample
def fitDens(data,
            locations,effsel,distmods,
            type='exp',
            mcmc=False,
            verbose=True):
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
    OUTPUT:
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    # Setup the density function and its initial parameters
    densfunc= _setup_densfunc(type)
    init= _setup_initparams_densfunc(type,data)
    # Setup the integration of the effective volume
    effsel, Rgrid, phigrid, zgrid= _setup_effvol(locations,effsel,distmods)
    # Get the data's R,phi,z
    dataR= data['RC_GALR_H']
    dataphi= data['RC_GALPHI_H']
    dataz= data['RC_GALZ_H']
    # Optimize
    out= optimize.fmin_powell(lambda x: _mloglike(x,densfunc,
                                                  dataR,dataphi,dataz,
                                                  effsel,Rgrid,phigrid,zgrid),
                              init,disp=verbose)
    print out
    return out

def _mloglike(*args,**kwargs):
    """Minus the log likelihood"""
    return -loglike(*args,**kwargs)

def loglike(params,densfunc,
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
def _setup_initparams_densfunc(type,data):
    """Return the initial parameters of the density for this type, might depend on the data"""
    if type.lower() == 'exp':
        return [1./3.,1./0.3]
    elif type.lower() == 'expplusconst':
        return [1./3.,1./0.3,numpy.log(0.1)]

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
                                            Ysun=0.,
                                            Zsun=define_rcsample._Z0)
        Rgrid.append(Rphiz[0])
        phigrid.append(Rphiz[1])
        zgrid.append(Rphiz[2])
    Rgrid= numpy.array(Rgrid)
    phigrid= numpy.array(phigrid)
    zgrid= numpy.array(zgrid)
    # Also need to multiply in distance factors
    effsel*= numpy.tile(ds**3.*(distmods[1]-distmods[0]),(effsel.shape[0],1))
    return (effsel,Rgrid,phigrid,zgrid)
