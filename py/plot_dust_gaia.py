###############################################################################
# plot_dust_gaia: plot the dust-map at X for Gaia
###############################################################################
import sys
import numpy
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot, bovy_coords
import dust
# nside to work at, 2048 is the max
_NSIDE= 2048
def plot_dust_gaia(dist,plotname):
    # Load the dust map
    green15map= dust.load_combined(dist,nest=True,nside_out=_NSIDE)
    dm= dust.dist2distmod(dist)
    green15map[green15map == healpy.UNSEEN]= -1.
    print "%i are NaN" % (numpy.sum(numpy.isnan(green15map)))
    #theta, phi= healpy.pixelfunc.pix2ang(_NSIDE,numpy.arange(healpy.pixelfunc.nside2npix(_NSIDE)),nest=True)
    #print (numpy.pi/2.-theta)[numpy.isnan(green15map)]/numpy.pi*180.
    #print phi[numpy.isnan(green15map)]/numpy.pi*180.
    # plot it
    healpy.visufunc.mollview(green15map,
                             nest=True,
                             xsize=4000,min=0.,
                             max=round(10.*(20.-dm-0.68))/10.,
                             format=r'$%g$',
                             cmap='gist_yarg',
                             title="",
                             unit='$A_G\,(\mathrm{mag})$')
    # Plot outline of Marshall et al.
    ls= numpy.linspace(-100.125,100.125,1001)
    healpy.visufunc.projplot(ls,ls*0.-10.125,'k--',lw=1.2,lonlat=True)
    healpy.visufunc.projplot(ls,ls*0.+10.125,'k--',lw=1.2,lonlat=True)
    bs= numpy.linspace(-10.125,10.125,1001)
    healpy.visufunc.projplot(bs*0.-100.125,bs,'k--',lw=1.2,lonlat=True)
    healpy.visufunc.projplot(bs*0.+100.125,bs,'k--',lw=1.2,lonlat=True)
    # Plot boundary of Green et al. map, dec > -30.
    ras= numpy.linspace(0.,360.,1001)
    decs= ras*0.-30.
    lbs= bovy_coords.radec_to_lb(ras,decs,degree=True)
    ls= lbs[:,0]
    bs= lbs[:,1]
    # rm those within Marshall et al.
    keepIndx= True-((ls > 259.875)*(numpy.fabs(bs) < 10.125)
                    +(ls < 100.125)*(numpy.fabs(bs) < 10.125))
    ls= ls[keepIndx]
    bs= bs[keepIndx]
    healpy.visufunc.projplot(ls,bs,'k-.',lw=1.2,lonlat=True)
    # Labels
    healpy.visufunc.projtext(350.,-8.,r'$\mathrm{Marshall\ et\ al.\ (2006)}$',
                             lonlat=True,size=13.)
    healpy.visufunc.projtext(10.,-60.,r'$\mathrm{Drimmel\ et\ al.\ (2003)}$',
                             lonlat=True,size=13.)
    healpy.visufunc.projtext(160.,40.,r'$\mathrm{Green\ et\ al.\ (2015)}$',
                             lonlat=True,size=13.)
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_dust_gaia(float(sys.argv[1]),sys.argv[2])
