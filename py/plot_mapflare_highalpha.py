###############################################################################
# plot_mapflare_highalpha.py: make of plot of the flaring of MAPs
###############################################################################
import os, os.path
import sys
import csv
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
def plot_mapflare(plotname):
    with open('../mapfits/tribrokenexpflare.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples_brexp= numpy.array(pickle.load(savefile))
    with open('../mapfits/tribrokenexpinvlinflare.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples_brexp_invlin= numpy.array(pickle.load(savefile))
    with open('../mapfits/tribrokenexplinflare.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples_brexp_lin= numpy.array(pickle.load(savefile))
    plotmaps= [19,26,32,39,45]
    bovy_plot.bovy_print(fig_width=8.,fig_height=9.*3.99/8.98)
    maps= define_rcsample.MAPs()
    cmap= cm.coolwarm
    overplot= False
    Rs= numpy.linspace(4.,14.,1001)
    # Setup for saving the profiles
    csvfile= open(os.path.join('..','out','mapflare_highalpha.csv'),'w')
    writer= csv.writer(csvfile,delimiter=',',quoting=csv.QUOTE_NONE)
    writer.writerow(['# Scale height profile for MAPs (Figure 13 in Bovy et al. 2016)'])
    writer.writerow(['# The first line lists the radii at which the scale height profiles'])
    writer.writerow(['# are evaluated'])
    writer.writerow(['# The rest of the file are the scale heights; the 0.025'])
    writer.writerow(['# lower limit and the 0.0975 upper limit (each 1 line)'])
    writer.writerow(['# Different MAPs are separated by a comment line'])
    writer.writerow(['{:.2f}'.format(x) for x in Rs])
    for ii, map in enumerate(maps.map()):
        if not ii in plotmaps: continue
        # Create all flaring profiles
        #Rmin= numpy.sort(map['RC_GALR_H'])[int(round(0.005*len(map)))]
        #Rmax= numpy.sort(map['RC_GALR_H'])[numpy.amin([len(map)-1,int(round(0.995*len(map)))])]
        samples= samples_brexp[ii,:,::_SKIP]
        samples_invlin= samples_brexp_invlin[ii,:,::_SKIP]
        samples_lin= samples_brexp_lin[ii,:,::_SKIP]
        nsamples= len(samples[0])
        tRs= numpy.tile(Rs,(nsamples,1)).T
        ldp= numpy.empty((len(Rs),nsamples))
        ldp= samples[1]*numpy.exp(samples[4]*(tRs-densprofiles._R0))
        ldp_invlin= samples_invlin[1]\
            *(1.+samples_invlin[4]*(tRs-densprofiles._R0)*(numpy.exp(1.)-1.))
        ldp_lin= samples_lin[1]\
            /(1.-(tRs-densprofiles._R0)*samples_lin[4]*(numpy.exp(1.)-1.))
        ldp= 1000./ldp # make it hz instead of its inverse
        ldp_invlin= 1000./ldp_invlin
        ldp_lin= 1000./ldp_lin
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
        #skipped because the median for the other profiles are indistinguishable for that with exponential flaring
        """
            pyplot.fill_between(Rs,
                                numpy.median(ldp_lin,axis=1)*offset,
                                numpy.median(ldp_invlin,axis=1)*offset,
                                color='0.65',
                                lw=0.,zorder=ii-1)
        """
        pyplot.fill_between(Rs,
                            numpy.sort(ldp,axis=1)[:,int(round(_SIGNIF*nsamples))]*offset,
                            numpy.sort(ldp,axis=1)[:,int(round((1.-_SIGNIF)*nsamples))]*offset,
                            color=cmap((tfeh+0.5)/0.4),
                            lw=0.,zorder=ii)
        line,= pyplot.plot(Rs,Rs*0.+300.*offset,color=cmap((tfeh+0.5)*0.95/0.5+0.05),
                           ls='--',lw=2.*0.8,zorder=ii+5)
        line.set_dashes([8,6])        
        overplot= True
        if ii == 19:
            bovy_plot.bovy_text(2.,
                                10.**5.6,
                                r'$[\mathrm{Fe/H}]$',size=16.,color='k')
        bovy_plot.bovy_text(2.,numpy.median(ldp,axis=1)[0]*offset,
                            r'$%+.1f$' % tfeh,size=16.,
                            color=cmap((tfeh+0.5)*0.95/0.5+0.05))
        writer.writerow(['# Low-alpha MAP w/ [Fe/H]=%g' % tfeh])
        writer.writerow(['{:.3f}'.format(x) for x in list(numpy.median(ldp,axis=1))])
        writer.writerow(['{:.3f}'.format(x) for x in list(numpy.sort(ldp,axis=1)[:,int(round(_SIGNIF*nsamples))])])
        writer.writerow(['{:.3f}'.format(x) for x in list(numpy.sort(ldp,axis=1)[:,int(round((1.-_SIGNIF)*nsamples))])])
    csvfile.close()
    bovy_plot.bovy_text(10.,10.**5.6,
                        r'$\mathrm{high-}[\alpha/\mathrm{Fe}]\ \mathrm{MAPs}$',
                        size=16.)
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_mapflare(sys.argv[1])
    
