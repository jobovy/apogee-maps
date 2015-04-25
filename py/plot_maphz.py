###############################################################################
# plot_maphz.py: plot h_Z in the feh,afe plane for MAPs
###############################################################################
import sys
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
from matplotlib import pyplot, cm
import define_rcsample
def plot_maphz(plotname):
    # Load the three fits
    with open('../mapfits/tribrokenexpflare.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples= numpy.array(pickle.load(savefile))
    with open('../mapfits/tribrokenexp.sav','rb') as savefile:
        bfnf= numpy.array(pickle.load(savefile))
        samplesnf= numpy.array(pickle.load(savefile))
    with open('../mapfits/tribrokenexpfixedflare.sav','rb') as savefile:
        bfff= numpy.array(pickle.load(savefile))
        samplesff= numpy.array(pickle.load(savefile))
    maps= define_rcsample.MAPs()
    plotthisz= numpy.zeros(len(bf))+numpy.nan
    plotthisze= numpy.zeros(len(bf))+numpy.nan
    for ii, map in enumerate(maps.map()):
        if numpy.nanmedian(numpy.exp(samples[ii,3])) < 5.:
            tmed= numpy.nanmedian(1./samplesnf[ii,1])
            terr= numpy.nanstd(1./samplesnf[ii,1])
        else:
            tmed= numpy.nanmedian(1./samplesff[ii,1])
            terr= numpy.nanstd(1./samplesff[ii,1])
        plotthisz[ii]= tmed
        plotthisze[ii]= terr
    plotthisz[plotthisze/plotthisz > 0.2]= numpy.nan
    bovy_plot.bovy_print()
    maps.plot(plotthisz*1000.,
              vmin=200.,vmax=1000.,
              minnstar=15,
              zlabel=r'$h_Z\,(\mathrm{pc})$',
              shrink=0.68)
    # Sequences
    haloc= define_rcsample.highalphalocus()
    bovy_plot.bovy_plot(haloc[:,0],haloc[:,1],'-',color='0.75',
                        lw=2.5,overplot=True)
    haloc= define_rcsample.lowalphalocus()
    haloc= haloc[(haloc[:,0] > -0.55)*(haloc[:,0] < 0.225)]
    bovy_plot.bovy_plot(haloc[:,0],haloc[:,1],'-',color='0.75',
                        lw=2.5,overplot=True)
    # Label
    #t= pyplot.text(-0.51,0.235,r'$\mathrm{single}$',
    #                size=16.,color='w')
    #t.set_bbox(dict(alpha=0.5,color=cm.coolwarm(0.),
    #                edgecolor='none'))
    #t= pyplot.text(-0.475,0.195,r'$\mathrm{exponential}$',
    #                size=16.,color='w')
    #t.set_bbox(dict(alpha=0.5,color=cm.coolwarm(0.),
    #                edgecolor='none'))
    pyplot.tight_layout()
    bovy_plot.bovy_end_print(plotname,dpi=300)
    return None

if __name__ == '__main__':
    plot_maphz(sys.argv[1])
