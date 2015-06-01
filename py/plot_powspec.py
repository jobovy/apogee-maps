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
_NSIDE= 2048
# magnitude limits for the survey
_GMIN= 3.
_GMAX= 20.
_NGSAMPLES= 10000
def plot_powspec(dist,basename,plotname):
    # Density
    G0= 0.68+dust.dist2distmod(dist)
    densname= basename+'_D%.1f_denscl.sav' % dist
    if not os.path.exists(densname):
        densmap= densprofiles.healpixelate(dist,densprofiles.expdisk,
                                           [1./3.,1./0.3],nside=_NSIDE,
                                           nest=False)
        denscl= healpy.sphtfunc.anafast(densmap,pol=False)
        densmap2= densprofiles.healpixelate(dist,densprofiles.expdisk,
                                            [1./2.,1./0.9],nside=_NSIDE,
                                            nest=False)
        denscl2= healpy.sphtfunc.anafast(densmap2,pol=False)
        ell= numpy.arange(len(denscl))
        save_pickles(densname,ell,denscl,densmap,denscl2,densmap2)
    else:
        with open(densname,'rb') as savefile:
            ell= pickle.load(savefile)
            denscl= pickle.load(savefile)
            densmap= pickle.load(savefile)
            denscl2= pickle.load(savefile)
            densmap2= pickle.load(savefile)
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
            combinedmcr2= pickle.load(savefile)
            bestfitloaded= True
    if not bestfitloaded:
        # do the best-fit
        combinedmap= dust.load_combined(dist,nest=False,nside_out=_NSIDE)
        print numpy.sum(numpy.isnan(combinedmap))
        combinedmap[numpy.isnan(combinedmap)]= 0.
        # Sample over the distribution of MG
        combinedmask= numpy.zeros_like(combinedmap)
        iso= gaia_rc.load_iso()
        Gsamples= gaia_rc.sample_Gdist(iso,n=_NGSAMPLES)
        print "Computing effective selection function"
        for jj in range(_NGSAMPLES):
            combinedmask+= ((combinedmap > (_GMIN-G0-Gsamples[jj]+0.68))\
                                *(combinedmap < (_GMAX-G0-Gsamples[jj]+0.68))).astype('float')
        combinedmask/= _NGSAMPLES
        print "Computing Cl of extinction map"
        combinedcl= healpy.sphtfunc.anafast(combinedmap,pol=False)
        print "Computing cross of extinction map w/ densmap"
        combinedcr= healpy.sphtfunc.anafast(combinedmap,map2=densmap,pol=False)
        print "Computing Cl of effective selection function"
        combinedmcl= healpy.sphtfunc.anafast(combinedmask,pol=False)
        print "Computing cross of effective selection function w/ densmap"
        combinedmcr= healpy.sphtfunc.anafast(combinedmask,map2=densmap,pol=False)
        print "Computing cross of effective selection function w/ densmap2"
        combinedmcr2= healpy.sphtfunc.anafast(combinedmask,map2=densmap2,pol=False)
        # Save
        save_pickles(combinedname,ell,combinedcl,combinedcr,
                     combinedmcl,combinedmcr,combinedmcr2)
        gc.collect()
    # Plot (2l+1)Cl!!
    # Can smooth the masked power spectrum, perhaps underplot the non-smoothed in gray
    # sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],numpy.log(combinedmcl)[1:],k=3,s=300.)
    # sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],numpy.log(numpy.fabs(combinedmcr))[1:],k=3,s=10000.)
    # First plot the power-spectrum, then the cross-correlation, then the
    # cumulative sum
    bovy_plot.bovy_print(fig_height=3.)
    yrange=[10.**-12.,20.],
    line1= bovy_plot.bovy_plot(ell[1:],
                               (2.*ell[1:]+1.)*combinedcl[1:],
                               'k-',loglog=True,
                               ylabel=r'$(2l+1)\,C_l$',
                               xrange=[0.5,20000],
                               yrange=yrange,
                               zorder=3)
    line2= bovy_plot.bovy_plot(ell[2::2],
                               (2.*ell[2::2]+1.)*denscl[2::2],
                               'b-',overplot=True)
    line3= bovy_plot.bovy_plot(ell[1:],
                               (2.*ell[1:]+1.)*combinedmcl[1:],
                               'r-',overplot=True)
    # Add legend
    if dist == 5.:
        pyplot.legend((line2[0],line1[0],line3[0]),
                      (r'$\mathrm{exp.\ disk\ w/}$'+'\n'+r'$h_R = 3\,\mathrm{kpc},$'+'\n'+r'$h_Z = 0.3\,\mathrm{kpc}$',
                       r'$\mathrm{extinction\ map}$',
                       r'$\mathrm{effective\ selection\ function}$'),
                      loc='lower left',bbox_to_anchor=(.02,.02),
                      numpoints=8,
                      prop={'size':14},
                      frameon=False) 
        bovy_plot.bovy_text(r'$\mathrm{power\ spectrum}$',top_right=True,size=16.)
    nullfmt   = NullFormatter()         # no labels
    pyplot.gca().xaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotname)
    # Cross-correlation
    bovy_plot.bovy_print(fig_height=3.)
    line1= bovy_plot.bovy_plot(ell[1:],
                               (2.*ell[1:]+1.)*numpy.fabs(combinedcr[1:]),
                               'k-',loglog=True,
                               ylabel=r'$(2l+1)\,C^{\mathrm{cross}}_l$',
                               xrange=[0.5,20000],
                               yrange=yrange,
                               zorder=1)
    line2= bovy_plot.bovy_plot(ell[1:],
                               (2.*ell[1:]+1.)*numpy.fabs(combinedmcr[1:]),
                               'r-',overplot=True,zorder=2)
    # Add legend
    if dist == 5.:
        pyplot.legend((line1[0],line2[0]),
                      (r'$\mathrm{extinction\ map}$',
                       r'$\mathrm{effective\ selection}$'+'\n'+r'$\mathrm{function}$'),
                      loc='lower left',bbox_to_anchor=(.02,.02),
                      numpoints=8,
                      prop={'size':14},
                      frameon=False) 
        bovy_plot.bovy_text(r'$\mathrm{cross\ power\ spectrum\ w/\ density}$',top_right=True,size=16.)
    #sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],
    #                                 numpy.log(numpy.fabs(combinedmcr))[1:],
    #                                 k=3,s=100000.)
    #bovy_plot.bovy_plot(ell[1:],
    #                    10.*(2.*ell[1:]+1.)*numpy.exp(sp(numpy.log(ell[1:]))),
    #                    'r-',overplot=True,zorder=2)
    pyplot.gca().xaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotname.replace('powspec','crosspowspec'))
    effvol= numpy.sum((2.*ell+1.)*combinedmcr)
    #effvol2= numpy.sum((2.*ell+1.)*combinedmcr2)
    matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
    line1= bovy_plot.bovy_plot(ell[1:],
                               numpy.fabs(numpy.log(effvol)
                                          -numpy.log(numpy.cumsum((2.*ell+1.)*combinedmcr)))[1:],
                               'k-',loglog=True,
                               xlabel=r'$l$',
                               ylabel=r'$\left|\Delta\ln\sum_{l}\sum_{m}\nu_{*,lm}\,\textswab{S}_{lm}\right|$',
                               xrange=[0.5,20000],
                               yrange=[2.*10.**-13.,20.],
                               zorder=3)
    """
    line2= bovy_plot.bovy_plot(ell[1:],
                               numpy.fabs(numpy.log(effvol)
                                          -numpy.log(numpy.cumsum((2.*ell+1.)*combinedmcr))
                                          -numpy.log(effvol2)
                                          +numpy.log(numpy.cumsum((2.*ell+1.)*combinedmcr2)))[1:],
                               'k--',loglog=True,
                               overplot=True,zorder=2)
    # Add legend
    pyplot.legend((line1[0],line2[0]),
                  (r'$\mathrm{exp.\ disk\ top\ panel}$',
                   r'$\mathrm{relative\ wrt\ exp.\ disk\ w/}$'+'\n'+
                   r'$h_R = 2\,\mathrm{kpc}, h_Z = 0.9\,\mathrm{kpc}$'),
                   loc='lower right',#bbox_to_anchor=(.91,.375),
                   numpoints=8,
                   prop={'size':14},
                   frameon=False) 
    """
    bovy_plot.bovy_end_print(plotname.replace('powspec','cumulcrosspowspec')) 
    return None

if __name__ == '__main__':
    plot_powspec(float(sys.argv[1]), # distance
                 sys.argv[2], # basename of pickles
                 sys.argv[3]) # plotfilename
