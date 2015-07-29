###############################################################################
# plot_distanceintegral_final.py: make a final overview plot of the distance
#                                 integral
###############################################################################
import sys
import pickle
import numpy
from scipy import signal, integrate, interpolate
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from galpy.util import bovy_plot
import dust
def plot_distanceintegral_final(plotname):
    # Reload full area
    with open('../savs/distInt.sav','rb') as savefile:
        area_full= pickle.load(savefile)
    # Reload full area w/o center
    with open('../savs/distIntRmcenter.sav','rb') as savefile:
        area_rmcenter= pickle.load(savefile)
    # Calculate PSD of each
    psdx_full, psd_full= \
        signal.periodogram(area_full*dust._GREEN15DISTS**3./numpy.sum(area_full*dust._GREEN15DISTS**3.),
                                      fs=1./(dust._GREEN15DISTMODS[1]-dust._GREEN15DISTMODS[0]),
                           detrend=lambda x: x,scaling='spectrum')
    psdx_rmcenter, psd_rmcenter= \
        signal.periodogram(area_rmcenter*dust._GREEN15DISTS**3./numpy.sum(area_rmcenter*dust._GREEN15DISTS**3.),
                                      fs=1./(dust._GREEN15DISTMODS[1]-dust._GREEN15DISTMODS[0]),
                           detrend=lambda x: x,scaling='spectrum')
    bovy_plot.bovy_print(fig_height=5.5)
    matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
    line1= bovy_plot.bovy_plot(psdx_full[1:],numpy.sqrt(psd_full[1:]),
                               'k-',loglog=True,
                               xlabel=r'$\mathrm{distance\ resolution}\ k_\mu\,(\mathrm{mag}^{-1})$',
                               ylabel=r'$\mathrm{effective\ volume\ error}\ \sqrt{P_k}$',
                               xrange=[0.04,20.],
                               yrange=[10**-11.,5.])
    line2= bovy_plot.bovy_plot(psdx_rmcenter[1:],numpy.sqrt(psd_rmcenter[1:]),
                               'r-',overplot=True)
    bovy_plot.bovy_plot([1.,10.],[6.*10.**-4.,6.*10.**-7.],'k--',overplot=True)
    bovy_plot.bovy_plot([1.,10.],[2.*10.**-5.,2.*10.**-10.],
                        'r--',overplot=True)
    pyplot.legend((line1[0],line2[0]),
                  (r'$\mathrm{full\ sky}$',
                   r'$\mathrm{excluding}\ |180^\circ-l| > 155^\circ, |b| < 25^\circ$'),
                  loc='upper right',#bbox_to_anchor=(.02,.02),
                  numpoints=8,
                  prop={'size':14},
                  frameon=False)
    bovy_plot.bovy_text(r'$\mathrm{normalized}\ D^3\,\nu_*(\mu|\theta)\,\textswab{S}(\mu)$',
                            bottom_left=True,size=16.)
    bovy_plot.bovy_end_print(plotname)
    return None

if __name__ == '__main__':
    plot_distanceintegral_final(sys.argv[1])
