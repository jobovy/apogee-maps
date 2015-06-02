###############################################################################
# plot_gaia_rcmag.py: plot the distribution of RC abs. mag. for Gaia G
###############################################################################
import sys
import numpy
from galpy.util import bovy_plot
import gaia_rc
def plot_gaia_rcmag(plotname):
    Gs= numpy.linspace(0.1,1.9,201)
    iso= gaia_rc.load_iso()
    ZG= gaia_rc.load_ZG(iso)
    pG= numpy.array([gaia_rc.Gdist(g,ZG) for g in Gs])
    bovy_plot.bovy_print(fig_height=3.)
    bovy_plot.bovy_plot(Gs,pG/numpy.sum(pG)/(Gs[1]-Gs[0]),
                        'k-',lw=1.5,
                        xrange=[1.,0.],
                        yrange=[0.,6.],
                        xlabel=r'$M_G$')
    Gmax= Gs[numpy.argmax(pG)]
    sGs= Gs[Gs < Gmax]
    spG= pG[Gs < Gmax]
    lGs= Gs[Gs > Gmax]
    lpG= pG[Gs > Gmax]
    f= (lGs[numpy.argmin(numpy.fabs(lpG-numpy.nanmax(pG)/2.))]-sGs[numpy.argmin(numpy.fabs(spG-numpy.nanmax(pG)/2.))])/2.355
    m,s= numpy.sum(Gs*pG)/numpy.sum(pG), numpy.sqrt(numpy.sum(Gs**2.*pG)/numpy.sum(pG)-(numpy.sum(Gs*pG)/numpy.sum(pG))**2.)
    print "peak = %f, mean = %f, std. dev. = %f, std.dev from fwhm = %f" % (Gmax,m,s,f)
    bovy_plot.bovy_end_print(plotname)
    return None

if __name__ == '__main__':
    plot_gaia_rcmag(sys.argv[1])
