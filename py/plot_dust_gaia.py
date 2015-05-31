###############################################################################
# plot_dust_gaia: plot the dust-map at X for Gaia
###############################################################################
import sys
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot
import dust
# nside to work at, 2048 is the max
_NSIDE= 2048
def plot_dust_gaia(dist,plotname):
    # Load the dust map
    green15map= dust.load_combined(dist,nest=True,nside_out=_NSIDE)
    dm= dust.dist2distmod(dist)
    green15map[green15map == healpy.UNSEEN]= -1.
    # plot it
    healpy.visufunc.mollview(green15map,
                             nest=True,
                             xsize=4000,min=0.,
                             max=round(10.*(20.-dm-0.68))/10.,
                             format=r'$%g$',
                             cmap='gist_yarg',
                             unit='$A_G\,(\mathrm{mag})$')
    bovy_plot.bovy_end_print(plotname)

if __name__ == '__main__':
    plot_dust_gaia(float(sys.argv[1]),sys.argv[2])
