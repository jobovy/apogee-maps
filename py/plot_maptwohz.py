###############################################################################
# plot_maptwohz.py: make a plot of 2nd scale height vs. first scale height
###############################################################################
import sys
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
from matplotlib import pyplot
import densprofiles
import define_rcsample
def plot_maphz(plotname):
    # Load the two fit
    with open('../mapfits/tribrokentwoexp.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples= numpy.array(pickle.load(savefile))
    maps= define_rcsample.MAPs()
    plotthisz1= numpy.zeros(len(bf))+numpy.nan
    plotthisz1e= numpy.zeros(len(bf))+numpy.nan
    plotthisz2= numpy.zeros(len(bf))+numpy.nan
    plotthisz2e= numpy.zeros(len(bf))+numpy.nan
    for ii, map in enumerate(maps.map()):
        hzindx= (True-numpy.isnan(samples[ii,4]))\
            *(True-numpy.isnan(samples[ii,5]))\
            *(densprofiles.ilogit(samples[ii,4]) > 0.15)
        tmed= numpy.median(1./samples[ii,1,hzindx])
        terr= numpy.std(1./samples[ii,1,hzindx])
        plotthisz1[ii]= tmed
        plotthisz1e[ii]= terr
        tmed= numpy.median(1./samples[ii,5,hzindx])
        terr= numpy.std(1./samples[ii,5,hzindx])
        plotthisz2[ii]= tmed
        plotthisz2e[ii]= terr
    plotthisz1[plotthisz1e/plotthisz1 > 0.5]= numpy.nan
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(plotthisz2*1000.,plotthisz1*1000.,'ko',
                        xrange=[0.,1200.],yrange=[0.,1200.],
                        xlabel=r'$2^\mathrm{nd}\ \mathrm{scale\ height\,(pc)}$',
                        ylabel=r'$1^\mathrm{st}\ \mathrm{scale\ height\,(pc)}$',zorder=2)
    bovy_plot.bovy_plot([0,1200],[0,1200],'k--',overplot=True,lw=2.)
    pyplot.errorbar(plotthisz2*1000.,plotthisz1*1000.,
                    yerr=plotthisz1e*1000.,
                    marker='o',color='k',ls='none',zorder=1)
    bovy_plot.bovy_end_print(plotname)
    return None

if __name__ == '__main__':
    plot_maphz(sys.argv[1])
