###############################################################################
# dust: tools to deal with the dust maps
###############################################################################
import os, os.path
import numpy
import h5py
import healpy
def dist2distmod(dist):
    """dist in kpc to distance modulus"""
    return 5.*numpy.log10(dist)+10.
def distmod2dist(distmod):
    """distance modulus to distance in kpc"""
    return 10.**(distmod/5.-2.)

_greendir= os.path.join(os.getenv('DUST_DIR'),'green15')
_GREEN15DISTMODS= numpy.linspace(4.,19.,31)
_GREEN15DISTS= distmod2dist(_GREEN15DISTMODS)
# we load these in the load_green15 function, but then re-use them
_PRELOADGREEN15= False
_PRELOADGREEN15SAMPLES= False
pix_info= None
dsamples= None
nsamples= None
best_fit= None

def load_green15(dist,nside_out=None,samples=False):
    """
    NAME:
       load_green15
    PURPOSE:
       load the Green et al. (2015) dust map at a given distance
    INPUT:
       dist - distance in kpc (nearest bin in the map will be returned)
       nside_out=  desired output nside (default is max in dust map)
       samples= (False) if True, return nsamples maps from the samples, otherwise the bestfit is returned
    OUTPUT:
       map(s) in HEALPIX nested scheme
    HISTORY:
       2015-03-04 - Written - Bovy (IAS)
    """
    # Load the data
    if not _PRELOADGREEN15:
        global pix_info
        global best_fit
        with h5py.File(os.path.join(_greendir,'dust-map-3d.h5'),'r') as greendata:
            pix_info= greendata['/pixel_info'][:]
            best_fit= greendata['/best_fit'][:]
    if not _PRELOADGREEN15SAMPLES and samples:
        global dsamples
        global nsamples
        with h5py.File(os.path.join(_greendir,'dust-map-3d.h5'),'r') as greendata:
            dsamples= greendata['/samples'][:]
            nsamples= dsamples.shape[0]
    # Distance pixel
    tpix= numpy.argmin(numpy.fabs(dist-_GREEN15DISTS))
    # Construct an empty map at the highest HEALPix resolution present in the map; code snippets adapted from http://argonaut.skymaps.info/usage
    nside_max= numpy.max(pix_info['nside'])
    npix= healpy.pixelfunc.nside2npix(nside_max)
    if samples:
        pix_val= numpy.empty((nsamples,npix),dtype='f8')
    else:
        pix_val= numpy.empty(npix,dtype='f8')
    pix_val[:] = numpy.nan
    # Fill the upsampled map
    for nside in numpy.unique(pix_info['nside']):
        # Get indices of all pixels at current nside level
        indx= pix_info['nside'] == nside
        # Extract A_H of each selected pixel
        if samples:
            pix_val_n= 0.460*dsamples[:,indx,tpix]
        else:
            pix_val_n= 0.460*best_fit[indx,tpix]
        # Determine nested index of each selected pixel in upsampled map
        mult_factor = (nside_max/nside)**2
        pix_idx_n = pix_info['healpix_index'][indx]*mult_factor
        # Write the selected pixels into the upsampled map
        for offset in range(mult_factor):
            if samples:               
                pix_val[:,pix_idx_n+offset] = pix_val_n
            else:
                pix_val[pix_idx_n+offset] = pix_val_n[:]
    # If the desired nside is less than the maximum nside in the map, degrade
    if not nside_out is None and nside_out < nside_max:
        if samples:
            for ii in range(nsamples):
                pix_val[ii]= healpy.pixelfunc.ud_grade(pix_val[ii],
                                               nside_out,pess=False,
                                               order_in='NEST', 
                                               order_out='NEST')
        else:
            pix_val= healpy.pixelfunc.ud_grade(pix_val,
                                               nside_out,pess=False,
                                               order_in='NEST', 
                                               order_out='NEST')
    return pix_val
