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
def plot_maphz(savefilename,plotname):
    with open(savefilename,'rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples= numpy.array(pickle.load(savefile))
        bf_g15= numpy.array(pickle.load(savefile))
        samples_g15= numpy.array(pickle.load(savefile))
        bf_zero= numpy.array(pickle.load(savefile))
        samples_zero= numpy.array(pickle.load(savefile))
    maps= define_rcsample.MAPs()
    plotthis= numpy.zeros(len(bf))+numpy.nan
    for ii, map in enumerate(maps.map()):
        tmed= numpy.median(1./samples[ii,1,True-numpy.isnan(1./samples[ii,1])])
        plotthis[ii]= tmed
    bovy_plot.bovy_print()
    maps.plot(plotthis*1000.,
              vmin=100.,vmax=1000.,
              minnstar=15,
              zlabel=r'$h_Z\,(\mathrm{pc})$',
              shrink=0.68)
    # Sequences
    #haloc= define_rcsample.highalphalocus()
    #bovy_plot.bovy_plot(haloc[:,0],haloc[:,1],'-',color='0.75',
    #                    lw=2.5,overplot=True)
    #haloc= define_rcsample.lowalphalocus()
    #haloc= haloc[(haloc[:,0] > -0.55)*(haloc[:,0] < 0.225)]
    #bovy_plot.bovy_plot(haloc[:,0],haloc[:,1],'-',color='0.75',
    #                    lw=2.5,overplot=True)
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
    plot_maphz(sys.argv[1],sys.argv[2])
