###############################################################################
# plot_broadsurfdens.py: make of plot of the surface-density of the broad bins
###############################################################################
import sys
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
from matplotlib import pyplot
import densprofiles
_SKIP= 10
_SIGNIF= 0.025
def plot_broadsurfdens(plotname):
    broads= ['lowlow','solar','highfeh','highalpha']
    #anorms= [2.,.25,0.01,50.]
    anorms= [0.004,.20,5.,50.]
    bovy_plot.bovy_print(fig_width=8.,fig_height=3.)
    overplot= False
    for ii, broad in enumerate(broads):
        # Restore the fits
        savename= '../broadfits/%s.sav' % broad
        with open(savename,'rb') as savefile:
            bf_exp= pickle.load(savefile)
            bf_brexp= pickle.load(savefile)
            bf_twoexp= pickle.load(savefile)
            ml_exp= pickle.load(savefile)
            ml_brexp= pickle.load(savefile)
            ml_twoexp= pickle.load(savefile)
            samples_exp= pickle.load(savefile)
            samples_brexp= pickle.load(savefile)
            samples_twoexp= pickle.load(savefile)
        # Create all density profiles
        Rs= numpy.linspace(4.,14.,1001)
        if broad.lower() == 'highalpha':
            samples= samples_exp[:,::_SKIP]
            nsamples= len(samples[0])
            tRs= numpy.tile(Rs,(nsamples,1)).T
            ihRin= numpy.tile(samples[0],(len(Rs),1))
            ldp= -ihRin*(tRs-densprofiles._R0)
        else:
            samples= samples_brexp[:,::_SKIP]
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
            # Rb < R0
            leRb= (tRs <= Rb)*(Rb < densprofiles._R0)
            ldp[leRb]= ihRin[leRb]*(tRs[leRb]-densprofiles._R0)\
                -ihRout[leRb]*(Rb[leRb]-densprofiles._R0)\
                -ihRin[leRb]*(Rb[leRb]-densprofiles._R0)
            gtRb= (tRs > Rb)*(Rb < densprofiles._R0)
            ldp[gtRb]= -ihRout[gtRb]*(tRs[gtRb]-densprofiles._R0)
        norm= numpy.exp(numpy.median(ldp,axis=1))[numpy.argmin(numpy.fabs(Rs-densprofiles._R0))]/anorms[ii]        
        bovy_plot.bovy_plot(Rs,numpy.exp(numpy.median(ldp,axis=1))/norm,
                            'k-',
                            lw=2.,overplot=overplot,
                            xlabel=r'$R\,(\mathrm{kpc})$',
                            ylabel=r'$\Sigma(R)$',
                            xrange=[0.,16.],
                            yrange=[0.0003,900.],
                            semilogy=True)
        pyplot.fill_between(Rs,
                            numpy.exp(numpy.sort(ldp,axis=1)[:,int(round(_SIGNIF*nsamples))])/norm,
                            numpy.exp(numpy.sort(ldp,axis=1)[:,int(round((1.-_SIGNIF)*nsamples))])/norm,
                            color='0.65',
                            lw=0.)
        overplot= True
    # Label
    labelx= 1.
    bovy_plot.bovy_text(labelx,5.,r'$\mathrm{high\ [Fe/H]}$',
                        size=15.,backgroundcolor='w')
    bovy_plot.bovy_text(labelx,0.085,r'$\mathrm{solar}$',
                        size=15.,backgroundcolor='w')
    bovy_plot.bovy_text(labelx,0.001,r'$\mathrm{low\ [Fe/H]}$',
                         size=15.,backgroundcolor='w')
    bovy_plot.bovy_text(labelx,150.,r'$\mathrm{high}\ [\alpha/\mathrm{Fe}]$',
                         size=15.,backgroundcolor='w')
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_broadsurfdens(sys.argv[1])
