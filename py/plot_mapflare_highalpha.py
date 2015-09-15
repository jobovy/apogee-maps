###############################################################################
# plot_mapflare_highalpha.py: make of plot of the flaring of MAPs
###############################################################################
import sys
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
from matplotlib import pyplot, cm
import densprofiles
import define_rcsample
_SKIP= 10
_SIGNIF= 0.025
def plot_mapflare(plotname,linflare=False):
    if not linflare:
        with open('../mapfits/tribrokenexpflare.sav','rb') as savefile:
            bf= numpy.array(pickle.load(savefile))
            samples_brexp= numpy.array(pickle.load(savefile))
    else:
        with open('../mapfits/tribrokenexplinflare.sav','rb') as savefile:
            bf= numpy.array(pickle.load(savefile))
            samples_brexp= numpy.array(pickle.load(savefile))
    plotmaps= [19,26,32,39,45]
    bovy_plot.bovy_print(fig_width=8.,fig_height=9.*3.99/8.98)
    maps= define_rcsample.MAPs()
    cmap= cm.coolwarm
    overplot= False
    for ii, map in enumerate(maps.map()):
        if not ii in plotmaps: continue
        # Create all flaring profiles
        #Rmin= numpy.sort(map['RC_GALR_H'])[int(round(0.005*len(map)))]
        #Rmax= numpy.sort(map['RC_GALR_H'])[numpy.amin([len(map)-1,int(round(0.995*len(map)))])]
        Rs= numpy.linspace(4.,14.,1001)
        samples= samples_brexp[ii,:,::_SKIP]
        nsamples= len(samples[0])
        tRs= numpy.tile(Rs,(nsamples,1)).T
        ldp= numpy.empty((len(Rs),nsamples))
        if linflare:
            ldp= samples[1]\
                *(1.+samples[4]*(tRs-densprofiles._R0)*(numpy.exp(1.)-1.))
        else:
            ldp= samples[1]*numpy.exp(samples[4]*(tRs-densprofiles._R0))
        ldp= 1000./ldp # make it hz instead of its inverse
        # Label and relative normalization
        tfeh= round(numpy.median(map['FE_H'])*20.)/20.
        if tfeh == 0.25: tfeh= 0.3
        if tfeh == -0.0: tfeh= 0.0
        offset= 10.**(-6.*(tfeh+0.1))-0.4
        if tfeh < -0.1: offset/= 1.5
        if tfeh < -0.2: offset/= 1.15
        print ii, tfeh, len(map), offset
        bovy_plot.bovy_plot(Rs,numpy.median(ldp,axis=1)*offset,
                            '-',
                            color=cmap((tfeh+0.5)*0.95/0.5+0.05),
                            lw=2.,overplot=overplot,
                            ylabel=r'$h_Z\,(\mathrm{pc})\times\mathrm{constant}$',
                            xrange=[0.,16.],
                            yrange=[10.**2.,10**5.99],
                            zorder=10+ii,
                            semilogy=True)
        pyplot.fill_between(Rs,
                            numpy.sort(ldp,axis=1)[:,int(round(_SIGNIF*nsamples))]*offset,
                            numpy.sort(ldp,axis=1)[:,int(round((1.-_SIGNIF)*nsamples))]*offset,
                            color=cmap((tfeh+0.5)/0.4),
                            lw=0.,zorder=ii)
        pyplot.plot(Rs,Rs*0.+300.*offset,color=cmap((tfeh+0.5)*0.95/0.5+0.05),
                    ls='--',lw=2.*0.8,zorder=ii+5)
        overplot= True
        if ii == 19:
            bovy_plot.bovy_text(2.,
                                10.**5.6,
                                r'$[\mathrm{Fe/H}]$',size=16.,color='k')
        bovy_plot.bovy_text(2.,numpy.median(ldp,axis=1)[0]*offset,
                            r'$%+.1f$' % tfeh,size=16.,
                            color=cmap((tfeh+0.5)*0.95/0.5+0.05))
    bovy_plot.bovy_text(10.,10.**5.6,
                        r'$\mathrm{high-}[\alpha/\mathrm{Fe}]\ \mathrm{MAPs}$',
                        size=16.)
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_mapflare(sys.argv[1],linflare=len(sys.argv) > 2)
    
