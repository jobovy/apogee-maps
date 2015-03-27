###############################################################################
# compareDataModel.py: module with data/model comparisons
###############################################################################
import os, os.path
import pickle
import numpy
from scipy import interpolate
from galpy.util import bovy_coords
from fitDens import _setup_densfunc
import define_rcsample
def predict_spacedist(params,
                      locations,effsel,distmods,
                      type='exp',
                      coord='Z'):
    """
    NAME:
       predict_spacedist
    PURPOSE:
       predict the spatial distribution
    INPUT:
       params - parameters of the density profile
       locations - locations of the APOGEE effective selection function to consider
       effsel - array (nloc,nD) of the effective selection function, includes area of the field
       distmods - grid of distance moduli on which the effective selection function is pre-computed
       type= ('exp') type of density profile to fit      
       coord= ('dm', 'X', or 'Z')
    OUTPUT:
       (R,model(R))
    HISTORY:
       2015-03-26 - Written - Bovy (IAS)
    """
    if coord.lower() == 'x':
        # Grid in X
        Xs= numpy.linspace(0.,20.,301)
    elif coord.lower() == 'z':
        # Grid in X
        Xs= numpy.linspace(0.,20.,301)
    elif coord.lower() == 'dm':
        # Grid in X
        Xs= numpy.linspace(7.,15.5,301)
    # Setup the density function
    rdensfunc= _setup_densfunc(type)
    densfunc= lambda x,y,z: rdensfunc(x,y,z,params=params)
    # Restore the APOGEE selection function (assumed pre-computed)
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    # Now compute the necessary coordinate transformations
    ds= 10.**(distmods/5-2.)
    Rgrid, phigrid, zgrid, Xgrid= [], [], [], []
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
        Xgrid.append(Rphiz[0]*numpy.cos(Rphiz[1]))
    Rgrid= numpy.array(Rgrid)
    phigrid= numpy.array(phigrid)
    zgrid= numpy.array(zgrid)
    Xgrid= numpy.array(Xgrid)
    # Now compute rate(R) for each location and combine
    effsel*= numpy.tile(ds**3.*(distmods[1]-distmods[0]),(effsel.shape[0],1))
    tdens= densfunc(Rgrid,phigrid,zgrid)
    rate= tdens*effsel
    out= numpy.zeros((len(locations),len(Xs)))
    for ii in range(len(locations)):
        if coord.lower() == 'x':
            # Jacobian
            tjac= numpy.fabs((numpy.roll(distmods,-1)-distmods)/\
                                 (numpy.roll(Xgrid[ii],-1)-Xgrid[ii]))
            tjac[-1]= tjac[-2]
            tXs= Xgrid[ii,rate[ii] > 0.]
        elif coord.lower() == 'z':
            # Jacobian
            tjac= numpy.fabs((numpy.roll(distmods,-1)-distmods)/\
                                 (numpy.roll(zgrid[ii],-1)-zgrid[ii]))
            tjac[-1]= tjac[-2]
            tXs= zgrid[ii,rate[ii] > 0.]
        elif coord.lower() == 'dm':
            # Jacobian
            tjac= numpy.ones_like(Xs)
            tXs= distmods[rate[ii] > 0.]
        sindx= numpy.argsort(tXs)
        tXs= tXs[sindx]
        trate= rate[ii,rate[ii] > 0.][sindx]
        tjac= tjac[rate[ii] > 0.][sindx]
        ipthis= numpy.log(trate*tjac+10.**-8.)
        baseline= numpy.polynomial.Polynomial.fit(tXs,ipthis,4)
        ipthis= ipthis/baseline(tXs)
        sp= interpolate.InterpolatedUnivariateSpline(tXs,ipthis,k=3)
        tindx= (Xs >= numpy.amin(tXs))\
            *(Xs <= numpy.amax(tXs))
        out[ii,tindx]= (numpy.exp(sp(Xs[tindx])*baseline(Xs[tindx]))-10.**-8.)
    out[numpy.isinf(out)]= 0.
    return (Xs,out.sum(axis=0))


