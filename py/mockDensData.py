###############################################################################
# mockDensData.py: generate mock data following a given density
###############################################################################
import os, os.path
import pickle
import multiprocessing
from optparse import OptionParser
import numpy
from scipy import ndimage
import fitsio
from galpy.util import bovy_coords, multi
import mwdust
import define_rcsample
import fitDens
import densprofiles
dmap= None
dmapg15= None
apo= None
def generate(locations,
             type='exp',
             sample='lowlow',
             extmap='green15',
             nls=101,
             nmock=1000,
             H0=-1.49,
             _dmapg15=None,             
             ncpu=1):
    """
    NAME:
       generate
    PURPOSE:
       generate mock data following a given density
    INPUT:
       locations - locations to be included in the sample
       type= ('exp') type of density profile to sample from
       sample= ('lowlow') for selecting mock parameters
       extmap= ('green15') extinction map to use ('marshall06' and others use Green15 to fill in unobserved regions)
       nls= (101) number of longitude bins to use for each field
       nmock= (1000) number of mock data points to generate
       H0= (-1.49) absolute magnitude (can be array w/ sampling spread)
       ncpu= (1) number of cpus to use to compute the probability
    OUTPUT:
       mockdata recarray with tags 'RC_GALR_H', 'RC_GALPHI_H', 'RC_GALZ_H'
    HISTORY:
       2015-04-03 - Written - Bovy (IAS)
    """
    if isinstance(H0,float): H0= [H0]
    # Setup the density function and its initial parameters
    rdensfunc= fitDens._setup_densfunc(type)
    mockparams= _setup_mockparams_densfunc(type,sample)
    densfunc= lambda x,y,z: rdensfunc(x,y,z,params=mockparams)   
    # Setup the extinction map
    global dmap
    global dmapg15
    if _dmapg15 is None: dmapg15= mwdust.Green15(filter='2MASS H')
    else: dmapg15= _dmapg15
    if isinstance(extmap,mwdust.DustMap3D.DustMap3D):
        dmap= extmap
    elif extmap.lower() == 'green15':
        dmap= dmapg15
    elif extmap.lower() == 'marshall06':
        dmap= mwdust.Marshall06(filter='2MASS H')
    elif extmap.lower() == 'sale14':
        dmap= mwdust.Sale14(filter='2MASS H')
    elif extmap.lower() == 'drimmel03':
        dmap= mwdust.Drimmel03(filter='2MASS H')
    # Use brute-force rejection sampling to make no approximations
    # First need to estimate the max probability to use in rejection;
    # Loop through all locations and compute sampling probability on grid in 
    # (l,b,D)
    # First restore the APOGEE selection function (assumed pre-computed)
    global apo
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    # Now compute the necessary coordinate transformations and evaluate the 
    # maximum probability
    distmods= numpy.linspace(7.,15.5,301)
    ds= 10.**(distmods/5-2.)
    nbs= nls
    lnprobs= numpy.empty((len(locations),len(distmods),nbs,nls))
    radii= []
    lcens, bcens= [], []
    lnprobs= multi.parallel_map(lambda x: _calc_lnprob(locations[x],nls,nbs,
                                                       ds,distmods,
                                                       H0,
                                                       densfunc),
                                range(len(locations)),
                                numcores=numpy.amin([len(locations),
                                                     multiprocessing.cpu_count(),ncpu]))
    lnprobs= numpy.array(lnprobs)
    for ll, loc in enumerate(locations):
        lcen, bcen= apo.glonGlat(loc)
        rad= apo.radius(loc)
        radii.append(rad) # save for later
        lcens.append(lcen[0])
        bcens.append(bcen[0])
    maxp= (numpy.exp(numpy.nanmax(lnprobs))-10.**-8.)*1.1 # Just to be sure
    # Now generate mock data using rejection sampling
    nout= 0
    arlocations= numpy.array(locations)
    arradii= numpy.array(radii)
    arlcens= numpy.array(lcens)
    arbcens= numpy.array(bcens)
    out= numpy.recarray((nmock,),
                        dtype=[('RC_DIST_H','f8'),
                               ('RC_DM_H','f8'),
                               ('RC_GALR_H','f8'),
                               ('RC_GALPHI_H','f8'),
                               ('RC_GALZ_H','f8')])
    while nout < nmock:
        nnew= 2*(nmock-nout)
        # nnew new locations
        locIndx= numpy.floor(numpy.random.uniform(size=nnew)*len(locations)).astype('int')
        newlocations= arlocations[locIndx]
        # Point within these locations
        newds_coord= numpy.random.uniform(size=nnew)
        newds= 10.**((newds_coord*(numpy.amax(distmods)-numpy.amin(distmods))\
            +numpy.amin(distmods))/5.-2.)
        newdls_coord= numpy.random.uniform(size=nnew)
        newdls= newdls_coord*2.*arradii[locIndx]\
            -arradii[locIndx]
        newdbs_coord= numpy.random.uniform(size=nnew)
        newdbs= newdbs_coord*2.*arradii[locIndx]\
            -arradii[locIndx]
        newr2s= newdls**2.+newdbs**2.
        keepIndx= newr2s < arradii[locIndx]**2.
        newlocations= newlocations[keepIndx]
        newds_coord= newds_coord[keepIndx]
        newdls_coord= newdls_coord[keepIndx]
        newdbs_coord= newdbs_coord[keepIndx]
        newds= newds[keepIndx]
        newdls= newdls[keepIndx]
        newdbs= newdbs[keepIndx]
        newls= newdls+arlcens[locIndx][keepIndx]
        newbs= newdbs+arbcens[locIndx][keepIndx]
        # Reject?
        tps= numpy.zeros_like(newds)
        for nloc in list(set(newlocations)):
            lindx= newlocations == nloc
            pindx= arlocations == nloc
            coord= numpy.array([newds_coord[lindx]*(len(distmods)-1.),
                                newdbs_coord[lindx]*(nbs-1.),
                                newdls_coord[lindx]*(nls-1.)])
            tps[lindx]= \
                numpy.exp(ndimage.interpolation.map_coordinates(\
                    lnprobs[pindx][0],
                    coord,cval=-10.,
                    order=1))-10.**-8.
        XYZ= bovy_coords.lbd_to_XYZ(newls,newbs,newds,degree=True)
        Rphiz= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],XYZ[:,1],XYZ[:,2],
                                            Xsun=define_rcsample._R0,
                                            Ysun=0.,
                                            Zsun=define_rcsample._Z0)
        testp= numpy.random.uniform(size=len(newds))*maxp
        keepIndx= tps > testp 
        if numpy.sum(keepIndx) > nmock-nout:
            rangeIndx= numpy.zeros(len(keepIndx),dtype='int')
            rangeIndx[keepIndx]= numpy.arange(numpy.sum(keepIndx))
            keepIndx*= (rangeIndx < nmock-nout)
        out['RC_DIST_H'][nout:nout+numpy.sum(keepIndx)]= newds[keepIndx]
        out['RC_DM_H'][nout:nout+numpy.sum(keepIndx)]= newds_coord[keepIndx]*(numpy.amax(distmods)-numpy.amin(distmods))\
            +numpy.amin(distmods)
        out['RC_GALR_H'][nout:nout+numpy.sum(keepIndx)]= Rphiz[0][keepIndx]
        out['RC_GALPHI_H'][nout:nout+numpy.sum(keepIndx)]= Rphiz[1][keepIndx]
        out['RC_GALZ_H'][nout:nout+numpy.sum(keepIndx)]= Rphiz[2][keepIndx]
        nout= nout+numpy.sum(keepIndx)
    return (out,lnprobs)

def _setup_mockparams_densfunc(type,sample):
    """Return the parameters of the mock density for this type"""
    if type.lower() == 'exp':
        if sample.lower() == 'lowlow':
            return [0.,1./0.3]
        elif sample.lower() == 'solar':
            return [1./3.,1./0.3]
        else:
            return [1./3.,1./0.3]
    elif type.lower() == 'expplusconst':
        if sample.lower() == 'lowlow':
            return [0.,1./0.3,numpy.log(0.1)]
        else:
            return [1./3.,1./0.3,numpy.log(0.1)]
    elif type.lower() == 'twoexp':
        return [1./3.,1./0.3,1./4.,1./0.5,densprofiles.logit(0.5)]
    elif type.lower() == 'brokenexp':
        if sample.lower() == 'lowlow':
            return [-0.2,1./.3,0.2,numpy.log(11.)]
        elif sample.lower() == 'solar':
            return [-1./6.,1./0.3,1./2.,numpy.log(14.)]
        else:
            return [-1./6.,1./0.3,1./2.,numpy.log(14.)]
    elif type.lower() == 'gaussexp':
        if sample.lower() == 'lowlow':
            return [.4,1./0.3,numpy.log(11.)]
        else:
            return [1./3.,1./0.3,numpy.log(10.)]

def _calc_lnprob(loc,nls,nbs,ds,distmods,H0,densfunc):
    lcen, bcen= apo.glonGlat(loc)
    rad= apo.radius(loc)
    ls= numpy.linspace(lcen-rad,lcen+rad,nls)
    bs= numpy.linspace(bcen-rad,bcen+rad,nbs)
    # Tile these
    tls= numpy.tile(ls,(len(ds),len(bs),1))
    tbs= numpy.swapaxes(numpy.tile(bs,(len(ds),len(ls),1)),1,2)
    tds= numpy.tile(ds,(len(ls),len(bs),1)).T
    XYZ= bovy_coords.lbd_to_XYZ(tls.flatten(),
                                tbs.flatten(),
                                tds.flatten(),
                                degree=True)
    Rphiz= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],XYZ[:,1],XYZ[:,2],
                                        Xsun=define_rcsample._R0,
                                        Ysun=0.,
                                        Zsun=define_rcsample._Z0)
    # Evaluate probability density
    tH= numpy.tile(distmods.T,(1,len(ls),len(bs),1))[0].T
    for ii in range(tH.shape[1]):
        for jj in range(tH.shape[2]):
            try:
                tH[:,ii,jj]+= dmap(ls[jj],bs[ii],ds)
            except (IndexError, TypeError,ValueError):
                try:
                    tH[:,ii,jj]+= dmapg15(ls[jj],bs[ii],ds)
                except IndexError: # assume zero outside
                    pass
    tH= tH.flatten()+H0[0]
    ps= densfunc(Rphiz[0],Rphiz[1],Rphiz[2])*apo(loc,tH)\
        *numpy.fabs(numpy.cos(tbs.flatten()/180.*numpy.pi))\
        *tds.flatten()**3.
    return numpy.log(numpy.reshape(ps,(len(distmods),nbs,nls))\
                         +10.**-8.)

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the mock datawill be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("--type",dest='type',default='exp',
                      help="Type of density profile")
    parser.add_option("--sample",dest='sample',default='lowlow',
                      help="Sample parameter for mock parameters")
    parser.add_option("--H0",dest='H0',default=-1.49,type='float',
                      help="RC absolute magnitude")
    parser.add_option("--nls",dest='nls',default=101,type='int',
                      help="Number of longitudes to bin each field in")
    parser.add_option("--nmock",dest='nmock',default=20000,type='int',
                      help="Number of mock samples to generate")
    # Dust map to use
    parser.add_option("--extmap",dest='extmap',default='green15',
                      help="Dust map to use ('Green15', 'Marshall03', 'Drimmel03', 'Sale14', or 'zero'")
    # Multiprocessing?
    parser.add_option("-m","--multi",dest='multi',default=1,type='int',
                      help="number of cpus to use")
    return parser

if __name__ == '__main__':
    parser= get_options()
    options, args= parser.parse_args()
    data= define_rcsample.get_rcsample()
    locations= list(set(list(data['LOCATION_ID'])))
    #locations= [4240,4242]
    out= generate(locations,
                  type=options.type,
                  sample=options.sample,
                  extmap=options.extmap,
                  nls=options.nls,
                  nmock=options.nmock,
                  H0=options.H0,
                  ncpu=options.multi)
    fitsio.write(args[0],out[0],clobber=True)
    
