###############################################################################
# plot_afefeh: the basic [a/Fe] vs. [Fe/H] plot for the data section
###############################################################################
import sys
import matplotlib
import numpy
from scipy import special
matplotlib.use('Agg')
from galpy.util import bovy_plot
from matplotlib import pyplot
import define_rcsample
def plot_afefeh(plotfilename):
    # Load the data
    data= define_rcsample.get_rcsample()
    # Plot the data
    bovy_plot.bovy_print()
    bovy_plot.scatterplot(data[define_rcsample._FEHTAG],
                          data[define_rcsample._AFETAG],
                          'k.',ms=.8,
                          levels=special.erf(numpy.arange(1,2)/numpy.sqrt(2.)),
                          xrange=[-.9,0.5],
                          yrange=[-0.15,0.35],
                          xlabel=r'$[\mathrm{Fe/H}]$',
                          ylabel=define_rcsample._AFELABEL)
    # Overplot sub-samples
    # low alpha, low feh
    lowfeh= define_rcsample._lowlow_lowfeh(0.)
    highfeh= define_rcsample._lowlow_highfeh(0.)
    pyplot.plot([lowfeh,lowfeh],[define_rcsample._lowlow_lowafe(lowfeh),
                                 define_rcsample._lowlow_highafe(lowfeh)],
                'k--',lw=2.)
    pyplot.plot([highfeh,highfeh],[define_rcsample._lowlow_lowafe(highfeh),
                                   define_rcsample._lowlow_highafe(highfeh)],
                'k--',lw=2.)
    pyplot.plot([lowfeh,highfeh],[define_rcsample._lowlow_lowafe(lowfeh),
                                  define_rcsample._lowlow_lowafe(highfeh)],
                'k--',lw=2.)
    pyplot.plot([lowfeh,highfeh],[define_rcsample._lowlow_highafe(lowfeh),
                                  define_rcsample._lowlow_highafe(highfeh)],
                'k--',lw=2.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_afefeh(sys.argv[1])
