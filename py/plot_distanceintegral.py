###############################################################################
# plot_distanceintegral.py: make a useful plot of the distance integration in
#                           the effective volume
###############################################################################
import sys
import os, os.path
import pickle
import numpy
from scipy import signal, integrate, interpolate
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import save_pickles, bovy_plot, multi
import multiprocessing
#from matplotlib import pyplot
#from matplotlib.ticker import NullFormatter
#from scipy import interpolate
import densprofiles
import dust
import gaia_rc
# nside to work at, 2048 is the max
_NSIDE= 2048
# magnitude limits for the survey
_GMIN= 3.
_GMAX= 20.
_NGSAMPLES= 10000
_DEGTORAD= numpy.pi/180.
def plot_distanceintegral(savename,plotname,rmcenter=False):
    if os.path.exists(savename):
        with open(savename,'rb') as savefile:
            area= pickle.load(savefile)
    else:
        # For samping over the absolute magnitude distribution
        iso= gaia_rc.load_iso()
        Gsamples= gaia_rc.sample_Gdist(iso,n=_NGSAMPLES)
        # l and b of the pixels
        theta, phi= healpy.pixelfunc.pix2ang(_NSIDE,
                                             numpy.arange(healpy.pixelfunc.nside2npix(_NSIDE)),
                                             nest=True)
        cosb= numpy.sin(theta)
        area= multi.parallel_map(lambda x: distanceIntegrand(\
                dust._GREEN15DISTS[x],cosb,Gsamples,rmcenter),
                                 range(len(dust._GREEN15DISTS)),
                                 numcores=numpy.amin([16,
                                                      len(dust._GREEN15DISTS),
                                                      multiprocessing.cpu_count()]))

        save_pickles(savename,area)
    # Plot the power spectrum
    if True:
        psdx, psd= signal.periodogram(area*dust._GREEN15DISTS**3./numpy.sum(area*dust._GREEN15DISTS**3.),
                                      fs=1./(dust._GREEN15DISTMODS[1]-dust._GREEN15DISTMODS[0]),
                                      detrend=lambda x: x,scaling='spectrum')
        bovy_plot.bovy_print(fig_height=3.)
        matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
        bovy_plot.bovy_plot(psdx[1:],psd[1:],
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
        bovy_plot.bovy_plot(dust._GREEN15DISTMODS,
                            area*dust._GREEN15DISTS**3.,
                            'k-',
                            xlabel=r'$\mu\,(\mathrm{mag}^{-1})$',
                            ylabel=r'$D^3\,\nu_*(\mu|\theta)\,\textswab{S}(\mu)$')
        bovy_plot.bovy_end_print(plotname)
    spl= interpolate.InterpolatedUnivariateSpline(dust._GREEN15DISTMODS,
                                                  area*dust._GREEN15DISTS**3.,
                                                  k=5)
    fthder= [spl.derivatives(dm)[4] for dm in dust._GREEN15DISTMODS]
    print "Simpson error= ", 0.5**4./180.*numpy.mean(numpy.fabs(fthder))/integrate.simps(area*dust._GREEN15DISTS**3.,dx=0.5)
    return None

def distanceIntegrand(dist,cosb,Gsamples,rmcenter):
    # Calculate the density
    densmap= densprofiles.healpixelate(dist,densprofiles.expdisk,
                                       [1./3.,1./0.3],nside=_NSIDE,
                                       nest=False)
    # Load the dust map
    combinedmap= dust.load_combined(dist,nest=False,nside_out=_NSIDE)
    # Sample over the distribution of MG
    combinedmask= numpy.zeros_like(combinedmap)
    G0= 0.68+dust.dist2distmod(dist)
    if dust.dist2distmod(dist) == 9.5 or dust.dist2distmod(dist) == 10.5:
        print numpy.sum(numpy.isnan(combinedmap))
    for jj in range(_NGSAMPLES):
        combinedmask+= ((combinedmap > (_GMIN-G0-Gsamples[jj]+0.68))\
                            *(combinedmap < (_GMAX-G0-Gsamples[jj]+0.68))).astype('float')
    combinedmask/= _NGSAMPLES
    if dust.dist2distmod(dist) == 9.5 or dust.dist2distmod(dist) == 10.5:
        print numpy.sum(combinedmask)
        print dust.dist2distmod(dist), numpy.sum(cosb*densmap*combinedmask)
    # If rmcenter, rm the center of the MW
    if rmcenter:
        theta, phi= healpy.pixelfunc.pix2ang(_NSIDE,
                                             numpy.arange(healpy.pixelfunc.nside2npix(_NSIDE)),
                                             nest=True)
        combinedmask[((phi < 20.*_DEGTORAD)+(phi > (360.-20.)*_DEGTORAD))\
                         *(numpy.fabs(numpy.pi/2.-theta) < 20.*_DEGTORAD)]= 0.
    # Compute cross correlation
    return float(numpy.sum(cosb*densmap*combinedmask))

if __name__ == '__main__':
    plot_distanceintegral(sys.argv[1], # savefilename
                          sys.argv[2], # plotfilename
                          rmcenter=len(sys.argv) > 3) # remove the central part
    
