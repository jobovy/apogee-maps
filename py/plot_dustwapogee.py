###############################################################################
# plot_dustwapogee: plot the dust-map at 5 kpc with the APOGEE fields in the 
#                   sample overlayed
###############################################################################
import sys
import numpy
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
import apogee.select.apogeeSelect
import dust
import define_rcsample
# nside to work at, 2048 is the max
_NSIDE=2048
def plot_dustwapogee(plotname):
    # Load the dust map
    green15map= dust.load_green15(5.,nest=True,nside_out=_NSIDE)
    green15map[green15map == healpy.UNSEEN]= -1.
    # plot it
    healpy.visufunc.mollview(green15map,
                             nest=True,
                             xsize=4000,min=0.,max=.8,
                             format=r'$%g$',
                             title='',
                             cmap='gist_yarg',
                             unit='$A_H\,(\mathrm{mag})$')
    # Load the RC data to get the fields
    data= define_rcsample.get_rcsample()
    loc_ids= numpy.array(list(set(data['LOCATION_ID'])))
    # Load the selection function, just to get the field centers
    apo= apogee.select.apogeeSelect(_justprocessobslog=True)
    theta= numpy.empty(len(loc_ids))
    phi= numpy.empty(len(loc_ids))
    for ii,loc_id in enumerate(loc_ids):
        tl, tb= apo.glonGlat(loc_id)
        theta[ii]= (90.-tb)/180.*numpy.pi
        phi[ii]= tl/180.*numpy.pi
    hib= numpy.fabs((numpy.pi/2.-theta)) > (8./180.*numpy.pi)
    healpy.visufunc.projplot(theta[hib],phi[hib],'o',ms=5.,mfc='none',mew=0.8,
                             mec='k')
    lowb= True-hib
    healpy.visufunc.projplot(theta[lowb],phi[lowb],'o',ms=5.,mfc='none',
                             mec='w',mew=0.8)
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_dustwapogee(sys.argv[1])
