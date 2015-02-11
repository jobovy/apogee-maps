import sys
import numpy
from scipy import special
import statsmodels.api as sm
from galpy.util import bovy_plot
import define_rcsample
def plot_rcdistancecomparison(plotfilename):
    # Get the sample
    rcdata= define_rcsample.get_rcsample()
    # Now plot the differece
    bovy_plot.bovy_print()
    levels= special.erf(numpy.arange(1,3)/numpy.sqrt(2.))
    bovy_plot.scatterplot(rcdata['RC_DIST'],
                          (rcdata['RC_DIST_H']-rcdata['RC_DIST'])/rcdata['RC_DIST'],
                          conditional=True,
                          levels=levels,
                          linestyle='none',color='k',marker=',',
                          xrange=[0.,7.49],yrange=[-0.075,0.075],
                          xlabel=r'$M_{K_s}\!-\!\mathrm{based\ distance\,(kpc)}$',
                          ylabel=r'$\mathrm{Fractional\ difference\ of}\ M_H\ \mathrm{vs.}\ M_{K_s}$',
                          onedhistx=True,bins=31)
    bovy_plot.bovy_plot([0.,10.],[0.,0.],'--',lw=2.,color='0.75',overplot=True)
    # Plot lowess
    lowess= sm.nonparametric.lowess
    z= lowess((rcdata['RC_DIST_H']-rcdata['RC_DIST'])/rcdata['RC_DIST'],
              rcdata['RC_DIST'],frac=.3)             
    bovy_plot.bovy_plot(z[:,0],z[:,1],'w--',lw=2.,overplot=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_rcdistancecomparison(sys.argv[1])
