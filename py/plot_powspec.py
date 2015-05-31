###############################################################################
# plot_powspec: plot the power-spectrum of a density distribution, dust, the 
#               essf, and the cross correlation of these at a given distance
###############################################################################
import gc
import sys
import os, os.path
import pickle
import numpy
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import save_pickles, bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from scipy import interpolate
import densprofiles
import dust
import gaia_rc
# nside to work at, 2048 is the max
_NSIDE=512 #2048
# magnitude limits for the survey
_GMIN= 3.
_GMAX= 20.
_NGSAMPLES= 10000
def plot_powspec(dist,basename,plotname):
    # Density
    densname= basename+'_D%.1f_denscl.sav' % dist
    if not os.path.exists(densname):
        densmap= densprofiles.healpixelate(dist,densprofiles.expdisk,
                                           [1./3.,1./0.3],nside=_NSIDE,
                                           nest=False)
        denscl= healpy.sphtfunc.anafast(densmap,pol=False)
        ell= numpy.arange(len(denscl))
        save_pickles(densname,ell,denscl,densmap)
    else:
        with open(densname,'rb') as savefile:
            ell= pickle.load(savefile)
            denscl= pickle.load(savefile)
#            densmap= pickle.load(savefile)
    # dust map Cl and cross-power with dens
    combinedname= basename+'_D%.1f_combinedcl.sav' % dist
    bestfitloaded= False
    if os.path.exists(combinedname):
        with open(combinedname,'rb') as savefile:
            ell= pickle.load(savefile)
            combinedcl= pickle.load(savefile)
            combinedcr= pickle.load(savefile)
            combinedmcl= pickle.load(savefile)
            combinedmcr= pickle.load(savefile)
            bestfitloaded= True
    if not bestfitloaded:
        # do the best-fit
        combinedmap= dust.load_combined(dist,nest=False,nside_out=_NSIDE)
        # Sample over the distribution of MG
        combinedmask= numpy.zeros_like(combinedmap)
        iso= gaia_rc.load_iso()
        Gsamples= gaia_rc.sample_Gdist(iso,n=_NGSAMPLES)
        for jj in range(_NGSAMPLES):
            combinedmask+= ((combinedmap > (_GMIN-Gsamples[jj]))\
                                *(combinedmap < (_GMAX-Gsamples[jj]))).astype('float')
        combinedmask/= _NGSAMPLES
        combinedcl= healpy.sphtfunc.anafast(combinedmap,pol=False)
        combinedcr= healpy.sphtfunc.anafast(combinedmap,map2=densmap,pol=False)
        combinedmcl= healpy.sphtfunc.anafast(combinedmask,pol=False)
        combinedmcr= healpy.sphtfunc.anafast(combinedmask,map2=densmap,pol=False)
        # Save
        save_pickles(combinedname,ell,combinedcl,combinedcr,
                     combinedmcl,combinedmcr)
        gc.collect()
    # Plot (2l+1)Cl!!
    # Can smooth the masked power spectrum, perhaps underplot the non-smoothed in gray
    # sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],numpy.log(combinedmcl)[1:],k=3,s=300.)
    # sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],numpy.log(numpy.fabs(combinedmcr))[1:],k=3,s=10000.)
    # Plot min, median, and max of samples
    # First plot the power-spectrum, then the cross-correlation, then the
    # cumulative sum
    bovy_plot.bovy_print(fig_height=3.)
    yrange=[10.**-12.,10.],
    bovy_plot.bovy_plot(ell[2::2],
                        (2.*ell[2::2]+1.)*denscl[2::2],
                        'b-',loglog=True,
                        ylabel=r'$(2l+1)\,C_l$',
                        xrange=[0.5,20000],
                        yrange=yrange)
    nullfmt   = NullFormatter()         # no labels
    pyplot.gca().xaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotname)
    # Cross-correlation
    bovy_plot.bovy_print(fig_height=3.)
    sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                     numpy.log(numpy.fabs(combinedmcr))[1:],
                                     k=3,s=100000.)
    bovy_plot.bovy_plot(ell[1:],
                        10.*(2.*ell[1:]+1.)*numpy.exp(sp(numpy.log(ell[1:]))),
                        'r-',overplot=True,zorder=2)
    pyplot.gca().xaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotname.replace('powspec','crosspowspec'))
    effvol= numpy.sum((2.*ell+1.)*combinedmcr)
    bovy_plot.bovy_plot(ell[1:],
                        numpy.fabs((effvol
                         -numpy.cumsum((2.*ell+1.)*combinedmcr))/effvol)[1:],
                        'k-',loglog=True,
                        xlabel=r'$l$',
                        ylabel=r'$\delta\sum_{l}\sum_{m}\nu_{*,lm}\,m^*_{lm}$',
                        xrange=[0.5,20000],
                        yrange=[2.*10.**-13.,20.],
                        zorder=3)
    bovy_plot.bovy_end_print(plotname.replace('powspec','cumulcrosspowspec')) 
    return None

if __name__ == '__main__':
    plot_powspec(float(sys.argv[1]), # distance
                 sys.argv[2], # basename of pickles
                 sys.argv[3]) # plotfilename
