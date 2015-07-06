###############################################################################
# plot_distanceintegral_smallpatch.py: make a useful plot of the distance 
#                                      integration in the effective volume
#                                      for a small patch of the sky with better
#                                      distance resolution
###############################################################################
import sys
import os, os.path
import pickle
import numpy
import h5py
from scipy import signal, integrate, interpolate
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import save_pickles, bovy_plot, multi
import multiprocessing
import densprofiles
import dust
import gaia_rc
# magnitude limits for the survey
_GMIN= 3.
_GMAX= 20.
_NGSAMPLES= 100
_HIRESGREEN15DISTMODS= numpy.linspace(4.,19.,61)
_HIRESGREEN15DISTS= dust.distmod2dist(_HIRESGREEN15DISTMODS)
_NSIDE_HIRES= 1024
def read_patch():
    with h5py.File(os.path.join(os.getenv('DATADIR'),'bovy','apogee-maps',
                                'l60_b03.h5'),'r') as greendata:
            pix_info= greendata['/pixel_info'][:]
            best_fit= greendata['/best_fit'][:]
    # Upgrade the Nside=512 pixels to Nside=1024
    healpix_index_out= numpy.zeros(len(pix_info)+3*numpy.sum(pix_info['nside'] == 512),
                                   dtype='int64')
    best_fit_out= numpy.zeros((len(pix_info)+3*numpy.sum(pix_info['nside'] == 512),61))
    healpix_index_out[:len(pix_info)-numpy.sum(pix_info['nside'] == 512)]=\
        pix_info['healpix_index'][pix_info['nside'] == 1024]
    best_fit_out[:len(pix_info)-numpy.sum(pix_info['nside'] == 512)]=\
        best_fit[pix_info['nside'] == 1024]
    for ii in range(numpy.sum(pix_info['nside'] == 512)):
        healpix_index_out[len(pix_info)-numpy.sum(pix_info['nside'] == 512)+4*ii:len(pix_info)-numpy.sum(pix_info['nside'] == 512)+4*ii+4]=\
            pix_info['healpix_index'][pix_info['nside'] == 512][ii]*4+numpy.arange(4)
        best_fit_out[len(pix_info)-numpy.sum(pix_info['nside'] == 512)+4*ii:len(pix_info)-numpy.sum(pix_info['nside'] == 512)+4*ii+4]=\
            best_fit[pix_info['nside'] == 512][ii]
    return (healpix_index_out,best_fit_out*0.460)

def plot_distanceintegral_smallpatch(savename,plotname):
    if os.path.exists(savename):
        with open(savename,'rb') as savefile:
            area= pickle.load(savefile)
    else:
        # Load the patch
        hpIndx, dmap= read_patch()
        # For samping over the absolute magnitude distribution
        iso= gaia_rc.load_iso()
        Gsamples= gaia_rc.sample_Gdist(iso,n=_NGSAMPLES)
        # l and b of the pixels
        theta, phi= healpy.pixelfunc.pix2ang(_NSIDE_HIRES,hpIndx,nest=True)
        cosb= numpy.sin(theta)
        area= multi.parallel_map(lambda x: distanceIntegrandHires(\
                _HIRESGREEN15DISTS[x],theta,phi,cosb,Gsamples),
                                 range(len(_HIRESGREEN15DISTS)),
                                 numcores=numpy.amin([16,
                                                      len(_HIRESGREEN15DISTS),
                                                      multiprocessing.cpu_count()]))

        save_pickles(savename,area)
    # Plot the power spectrum
    if True:
        psdx, psd= signal.periodogram(area*_HIRESGREEN15DISTS**3./numpy.sum(area*_HIRESGREEN15DISTS**3.),
                                      fs=(_HIRESGREEN15DISTMODS[1]-_HIRESGREEN15DISTMODS[0]),
                                      detrend=lambda x: x,scaling='spectrum')
        bovy_plot.bovy_print(fig_height=3.)
        matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
        bovy_plot.bovy_plot(psdx[1:]*2.*numpy.pi,psd[1:],
                            'k-',loglog=True,
                            xlabel=r'$2\pi\,k_\mu\,(\mathrm{mag}^{-1})$',
                            ylabel=r'$P_k$',
                            xrange=[0.04,4.])
        bovy_plot.bovy_text(r'$\mathrm{normalized}\ D^3\,\nu_*(\mu|\theta)\,\textswab{S}(\mu)$',
                            bottom_left=True,size=16.)
        bovy_plot.bovy_end_print(plotname)               
    else:
        bovy_plot.bovy_print(fig_height=3.)
        matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
        bovy_plot.bovy_plot(_HIRESGREEN15DISTMODS,
                            area*_HIRESGREEN15DISTS**3.,
                            'k-',
                            xlabel=r'$\mu\,(\mathrm{mag}^{-1})$',
                            ylabel=r'$D^3\,\nu_*(\mu|\theta)\,\textswab{S}(\mu)$')
        bovy_plot.bovy_end_print(plotname)
    spl= interpolate.InterpolatedUnivariateSpline(_HIRESGREEN15DISTMODS,
                                                  area*_HIRESGREEN15DISTS**3.,
                                                  k=5)
    fthder= [spl.derivatives(dm)[4] for dm in _HIRESGREEN15DISTMODS]
    print "Simpson error= ", 0.5**4./180.*numpy.mean(numpy.fabs(fthder))/integrate.simps(area*_HIRESGREEN15DISTS**3.,dx=0.5)
    return None

def distanceIntegrandHires(dist,theta,phi,cosb,Gsamples,hpIndx,dmap):
    # Calculate the density
    densmap= densprofiles.expdist(phi,numpy.pi/2.-theta,
                                  dist*numpy.ones(len(theta)),glon=True,
                                  params=[1./3.,1./0.3])
    # Sample over the distribution of MG
    combinedmask= numpy.zeros_like(dmap)
    G0= 0.68+dust.dist2distmod(dist)
    #if dust.dist2distmod(dist) == 9.5 or dust.dist2distmod(dist) == 10.5:
    #    print numpy.sum(numpy.isnan(combinedmap))
    for jj in range(_NGSAMPLES):
        combinedmask+= ((dmap > (_GMIN-G0-Gsamples[jj]+0.68))\
                            *(dmap < (_GMAX-G0-Gsamples[jj]+0.68))).astype('float')
    combinedmask/= _NGSAMPLES
    #if dust.dist2distmod(dist) == 9.5 or dust.dist2distmod(dist) == 10.5:
    #    print numpy.sum(combinedmask)
    #    print dust.dist2distmod(dist), numpy.sum(cosb*densmap*combinedmask)
    # Compute cross correlation
    return float(numpy.sum(cosb*densmap*combinedmask))

if __name__ == '__main__':
    plot_distanceintegral_smallpatch(sys.argv[1], # savefilename
                                     sys.argv[2]) # plotfilename
    
