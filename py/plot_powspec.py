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
# nside to work at, 2048 is the max
_NSIDE=2048
# magnitude limits for the survey
_HMIN= 7.
_HMAX= 12.8
def plot_powspec(dist,basename,plotname,plane=False):
    """plane is True limits to |b| <= 5"""
    # First calculate everything
    H0= -1.49+dust.dist2distmod(dist)
    if plane:
        npix= healpy.pixelfunc.nside2npix(_NSIDE)
        theta,phi= healpy.pixelfunc.pix2ang(_NSIDE,numpy.arange(npix),
                                            nest=False)
        b= (numpy.pi/2.-theta)*180./numpy.pi
    # Density
    densname= basename+'_D%.1f_denscl.sav' % dist
    if not os.path.exists(densname):
        densmap= densprofiles.healpixelate(dist,densprofiles.expdisk,
                                           [1./3.,1./0.3],nside=_NSIDE,
                                           nest=False)
        if plane:
            densmap[numpy.fabs(b) > 5.]= healpy.UNSEEN
        denscl= healpy.sphtfunc.anafast(densmap,pol=False)
        ell= numpy.arange(len(denscl))
        save_pickles(densname,ell,denscl,densmap)
    else:
        with open(densname,'rb') as savefile:
            ell= pickle.load(savefile)
            denscl= pickle.load(savefile)
#            densmap= pickle.load(savefile)
    # Pan-STARRS dust map Cl and cross-power with dens
    green15name= basename+'_D%.1f_green15cl.sav' % dist
    bestfitloaded= False
    samplesloaded= False
    if os.path.exists(green15name):
        with open(green15name,'rb') as savefile:
            ell= pickle.load(savefile)
            green15cl= pickle.load(savefile)
            green15cr= pickle.load(savefile)
            green15mcl= pickle.load(savefile)
            green15mcr= pickle.load(savefile)
            bestfitloaded= True
            try:
                samplescl= pickle.load(savefile)
                samplescr= pickle.load(savefile)
                samplesmcl= pickle.load(savefile)
                samplesmcr= pickle.load(savefile)
            except EOFError:
                samplesstart= 0
            else:
                samplesloaded= True
            try:
                samplesstart= pickle.load(savefile)
            except EOFError:
                samplesstart= 20
            else:
                samplesloaded= True
    if not bestfitloaded:
        # First do the best-fit
        green15map= dust.load_green15(dist,nest=False,nside_out=_NSIDE)
        if plane:
            green15map[numpy.fabs(b) > 5.]= healpy.UNSEEN
        green15mask= ((green15map > (_HMIN-H0))\
                          *(green15map < (_HMAX-H0))).astype('float')
        green15cl= healpy.sphtfunc.anafast(green15map,pol=False)
        green15cr= healpy.sphtfunc.anafast(green15map,map2=densmap,pol=False)
        green15mcl= healpy.sphtfunc.anafast(green15mask,pol=False)
        green15mcr= healpy.sphtfunc.anafast(green15mask,map2=densmap,pol=False)
        # Save
        save_pickles(green15name,ell,green15cl,green15cr,
                     green15mcl,green15mcr)
        gc.collect()
    # Now work on the samples
    nsamples= 20
    if not samplesloaded:
        samplescl= numpy.empty((nsamples,len(green15cl)))
        samplescr= numpy.empty((nsamples,len(green15cl)))
        samplesmcl= numpy.empty((nsamples,len(green15cl)))
        samplesmcr= numpy.empty((nsamples,len(green15cl)))
    if not samplesloaded or samplesstart < nsamples:
        for samplenum in range(samplesstart,nsamples):
            print "Working on sample %i / %i ..." % (samplenum+1,nsamples)
            green15maps= dust.load_green15(dist,nest=False,nside_out=_NSIDE,
                                           samples=True,samplenum=samplenum)
            if plane:
                green15maps[numpy.fabs(b) > 5.]= healpy.UNSEEN
            green15masks= ((green15maps > (_HMIN-H0))\
                              *(green15maps < (_HMAX-H0))).astype('float')
            samplescl[samplenum]= healpy.sphtfunc.anafast(green15maps,
                                                          pol=False)
            gc.collect()
            samplescr[samplenum]= healpy.sphtfunc.anafast(green15maps,
                                                          map2=densmap,
                                                          pol=False)
            gc.collect()
            samplesmcl[samplenum]= healpy.sphtfunc.anafast(green15masks,
                                                           pol=False)
            gc.collect()
            samplesmcr[samplenum]= healpy.sphtfunc.anafast(green15masks,
                                                           map2=densmap,
                                                           pol=False)
            gc.collect()
            # Save, here in case it crashes after a few samples
            save_pickles(green15name,ell,green15cl,green15cr,
                         green15mcl,green15mcr,
                         samplescl,samplescr,samplesmcl,samplesmcr,
                         samplenum+1)
    # Plot (2l+1)Cl!!
    # Can smooth the masked power spectrum, perhaps underplot the non-smoothed in gray
    # sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],numpy.log(green15mcl)[1:],k=3,s=300.)
    # sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],numpy.log(numpy.fabs(green15mcr))[1:],k=3,s=10000.)
    # Plot min, median, and max of samples
    # First plot the power-spectrum, then the cross-correlation, then the
    # cumulative sum
    bovy_plot.bovy_print(fig_height=3.)
    yrange=[10.**-12.,10.],
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.median(samplescl[:,1:],axis=0),
                        'k-',loglog=True,
                        ylabel=r'$(2l+1)\,C_l$',
                        xrange=[0.5,20000],
                        yrange=yrange,
                        zorder=3)
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.amin(samplescl[:,1:],axis=0),
                        color='0.65',zorder=0,overplot=True)
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.amax(samplescl[:,1:],axis=0),
                        color='0.65',zorder=0,overplot=True)
    bovy_plot.bovy_plot(ell[2::2],
                        (2.*ell[2::2]+1.)*denscl[2::2],
                        'b-',overplot=True)                       
    sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                     numpy.log(numpy.median(samplesmcl,
                                                            axis=0))[1:],
                                     k=3,s=300.)
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.exp(sp(numpy.log(ell[1:]))),
                        'r-',overplot=True,zorder=2)
    spmin= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                        numpy.log(numpy.amin(samplesmcl,
                                                             axis=0))[1:],
                                        k=3,s=300.)
    spmax= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                        numpy.log(numpy.amax(samplesmcl,
                                                             axis=0))[1:],
                                        k=3,s=300.)
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.exp(spmin(numpy.log(ell[1:]))),
                        color='0.65',zorder=0,overplot=True)
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.exp(spmax(numpy.log(ell[1:]))),
                        color='0.65',zorder=0,overplot=True)
    nullfmt   = NullFormatter()         # no labels
    pyplot.gca().xaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotname)
    # Cross-correlation
    bovy_plot.bovy_print(fig_height=3.)
    sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                     numpy.log(numpy.median(numpy.fabs(samplescr[:,1:]),
                                                            axis=0)),
                                     k=3,s=10000.)
    bovy_plot.bovy_plot(ell[1:],
                        (2.*ell[1:]+1.)*numpy.exp(sp(numpy.log(ell[1:]))),
                        'k-',loglog=True,
                        ylabel=r'$(2l+1)\,C_l$',
                        xrange=[0.5,20000],
                        yrange=[10.**-12.,100.],
                        zorder=3)
    spmin= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                        numpy.log(numpy.amin(numpy.fabs(samplescr),
                                                             axis=0))[1:],
                                        k=3,s=100000.)
    spmax= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                        numpy.log(numpy.amax(numpy.fabs(samplescr),
                                                             axis=0))[1:],
                                        k=3,s=100000.)
    pyplot.fill_between(ell[1:],
                        (2.*ell[1:]+1.)*numpy.exp(spmin(numpy.log(ell[1:]))),
                        (2.*ell[1:]+1.)*numpy.exp(spmax(numpy.log(ell[1:]))),
                        color='0.65',zorder=1)
    sp= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                     numpy.log(numpy.median(numpy.fabs(samplesmcr),
                                                            axis=0))[1:],
                                     k=3,s=100000.)
    bovy_plot.bovy_plot(ell[1:],
                        10.*(2.*ell[1:]+1.)*numpy.exp(sp(numpy.log(ell[1:]))),
                        'r-',overplot=True,zorder=2)
    spmin= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                        numpy.log(numpy.amin(numpy.fabs(samplesmcr),
                                                             axis=0))[1:],
                                        k=3,s=100000.)
    spmax= interpolate.UnivariateSpline(numpy.log(ell)[1:],
                                        numpy.log(numpy.amax(numpy.fabs(samplesmcr),
                                                             axis=0))[1:],
                                        k=3,s=100000.)
    pyplot.fill_between(ell[1:],
                        10.*(2.*ell[1:]+1.)*numpy.exp(spmin(numpy.log(ell[1:]))),
                        10.*(2.*ell[1:]+1.)*numpy.exp(spmax(numpy.log(ell[1:]))),
                        color='0.35',zorder=0)
    pyplot.gca().xaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotname.replace('powspec','crosspowspec'))
    effvol= numpy.sum((2.*ell+1.)*numpy.median(samplesmcr,axis=0))
    bovy_plot.bovy_plot(ell[1:],
                        numpy.fabs((effvol
                         -numpy.cumsum((2.*ell+1.)*numpy.median(samplesmcr,axis=0)))/effvol)[1:],
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
                 sys.argv[3],
                 plane=len(sys.argv) > 4) # plotfilename
