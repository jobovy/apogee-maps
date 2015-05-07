###############################################################################
# plot_mapsurfdens_highalpha.py: make of plot of the surface-density of the 
#                                high alpha stars
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
def plot_mapsurfdens(plotname):
    with open('../mapfits/tribrokenexpflare.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples_brexp= numpy.array(pickle.load(savefile))
    plotmaps= [19,26,32,39,45]
    bovy_plot.bovy_print(fig_width=8.,fig_height=4.)
    maps= define_rcsample.MAPs()
    cmap= cm.coolwarm
    overplot= False
    for ii, map in enumerate(maps.map()):
        if not ii in plotmaps: continue
        # Create all density profiles
        Rs= numpy.linspace(4.,14.,1001)
        samples= samples_brexp[ii,:,::_SKIP]
        nsamples= len(samples[0])
        tRs= numpy.tile(Rs,(nsamples,1)).T
        ldp= numpy.empty((len(Rs),nsamples))
        Rb= numpy.tile(numpy.exp(samples[3]),(len(Rs),1))
        ihRin= numpy.tile(samples[0],(len(Rs),1))
        ihRout= numpy.tile(samples[2],(len(Rs),1))
        # Rb >= R0
        leRb= (tRs <= Rb)*(Rb >= densprofiles._R0)
        ldp[leRb]= ihRin[leRb]*(tRs[leRb]-densprofiles._R0)
        gtRb= (tRs > Rb)*(Rb >= densprofiles._R0)
        ldp[gtRb]= -ihRout[gtRb]*(tRs[gtRb]-densprofiles._R0)\
            +ihRout[gtRb]*(Rb[gtRb]-densprofiles._R0)\
            +ihRin[gtRb]*(Rb[gtRb]-densprofiles._R0)
        # Rb < R0, normalize outer at R0
        leRb= (tRs <= Rb)*(Rb < densprofiles._R0)
        ldp[leRb]= ihRin[leRb]*(tRs[leRb]-densprofiles._R0)\
            -ihRout[leRb]*(Rb[leRb]-densprofiles._R0)\
            -ihRin[leRb]*(Rb[leRb]-densprofiles._R0)
        gtRb= (tRs > Rb)*(Rb < densprofiles._R0)
        ldp[gtRb]= -ihRout[gtRb]*(tRs[gtRb]-densprofiles._R0)
        # Label and relative normalization
        tfeh= round(numpy.median(map['FE_H'])*20.)/20.
        if tfeh == 0.35: tfeh= 0.4
        if tfeh == -0.0: tfeh= 0.0
        print ii, tfeh, len(map)
        anorm= 10**(-10.*tfeh)
        if tfeh > 0.3: anorm= 10**(-12.*tfeh) 
        if tfeh < -0.4: anorm= 10**(-12.*tfeh)
        anorm= 1./anorm # re-order
        anorm/= 3.
        norm= numpy.exp(numpy.median(ldp,axis=1))[numpy.argmin(numpy.fabs(Rs-densprofiles._R0))]/anorm
        bovy_plot.bovy_plot(Rs,numpy.exp(numpy.median(ldp,axis=1))/norm,
                            '-',
                            color=cmap((tfeh+0.4)*0.95/0.5+0.05),
                            lw=2.,overplot=overplot,
                            ylabel=r'$\Sigma(R)$',
                            xrange=[0.,16.],
                            yrange=[0.000001,90.],
                            semilogy=True)
        pyplot.fill_between(Rs,
                            numpy.exp(numpy.sort(ldp,axis=1)[:,int(round(_SIGNIF*nsamples))])/norm,
                            numpy.exp(numpy.sort(ldp,axis=1)[:,int(round((1.-_SIGNIF)*nsamples))])/norm,
                            color=cmap((tfeh+0.4)/0.4),
                            lw=0.)
        overplot= True
        if ii == 19:
            bovy_plot.bovy_text(2.,
                                10.**1.,
                                r'$[\mathrm{Fe/H}]$',size=16.,color='k')
        bovy_plot.bovy_text(2.,(numpy.exp(numpy.median(ldp,axis=1))/norm)[0],
                            r'$%+.1f$' % tfeh,size=16.,
                            color=cmap((tfeh+0.4)*0.95/0.5+0.05))
    bovy_plot.bovy_text(10.,10.**1.,
                        r'$\mathrm{high-}[\alpha/\mathrm{Fe}]\ \mathrm{MAPs}$',
                        size=16.)
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_mapsurfdens(sys.argv[1])
    
