###############################################################################
# plot_distanceintegral.py: make a useful plot of the distance integration in
#                           the effective volume
###############################################################################
import sys
import os, os.path
import pickle
import numpy
from scipy import signal
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import save_pickles, bovy_plot
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
def plot_distanceintegral(savename,plotname):
    if os.path.exists(savename):
        with open(savename,'rb') as savefile:
            area= pickle.load(savefile)
    else:
        area= numpy.zeros_like(dust._GREEN15DISTS)
        # For samping over the absolute magnitude distribution
        iso= gaia_rc.load_iso()
        Gsamples= gaia_rc.sample_Gdist(iso,n=_NGSAMPLES)
        # l and b of the pixels
        theta, phi= healpy.pixelfunc.pix2ang(_NSIDE,
                                             numpy.arange(healpy.pixelfunc.nside2npix(_NSIDE)),
                                             nest=True)
        cosb= numpy.sin(theta)
        for ii, dist in enumerate(dust._GREEN15DISTS):
            print "Working on distance %i: %.1f kpc" % (ii,dist)
            # Calculate the density
            densmap= densprofiles.healpixelate(dist,densprofiles.expdisk,
                                               [1./3.,1./0.3],nside=_NSIDE,
                                               nest=False)
            # Load the dust map
            combinedmap= dust.load_combined(dist,nest=False,nside_out=_NSIDE)
            # Sample over the distribution of MG
            combinedmask= numpy.zeros_like(combinedmap)
            G0= 0.68+dust.dist2distmod(dist)
            for jj in range(_NGSAMPLES):
                combinedmask+= ((combinedmap > (_GMIN-G0-Gsamples[jj]+0.68))\
                                    *(combinedmap < (_GMAX-G0-Gsamples[jj]+0.68))).astype('float')
                combinedmask/= _NGSAMPLES
            # Compute cross correlation
            area[ii]= numpy.sum(cosb*densmap*combinedmask)
        save_pickles(savename,area)
    # Plot the power spectrum
    psdx, psd= signal.periodogram(area*dust._GREEN15DISTS**3./numpy.sum(area*dust._GREEN15DISTS**3.),
                                  fs=(dust._GREEN15DISTMODS[1]-dust._GREEN15DISTMODS[0]),
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
    return None

if __name__ == '__main__':
    plot_distanceintegral(sys.argv[1], # savefilename
                          sys.argv[2]) # plotfilename
    
