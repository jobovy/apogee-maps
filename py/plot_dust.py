###############################################################################
# plot_dust: plot the dust-map at X
###############################################################################
import sys
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
import dust
# nside to work at, 2048 is the max
_NSIDE=2048
def plot_dust(dist,plotname):
    # Load the dust map
    green15map= dust.load_green15(dist,nest=True,nside_out=_NSIDE)
    dm= dust.dist2distmod(dist)
    green15map[green15map == healpy.UNSEEN]= -1.
    # plot it
    healpy.visufunc.mollview(green15map,
                             nest=True,
                             xsize=4000,min=0.,
                             max=round(10.*(12.8-dm+1.49))/10.,
                             format=r'$%g$',
                             title=r'$\mathrm{Dust\ at\ %.1f\,kpc}$' % dist,
                             cmap='gist_yarg',
                             unit='$A_H\,(\mathrm{mag})$')
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_dust(float(sys.argv[1]),sys.argv[2])
