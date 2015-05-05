###############################################################################
# plot_spatial_broad.py: plot the spatial distribution of broad samples
###############################################################################
import sys
from galpy.util import bovy_plot
import define_rcsample
def plot_spatial_broad(plotname):
    # Load subsamples, plot
    load_funcs= [define_rcsample.get_lowlowsample,
                 define_rcsample.get_highalphasample,
                 define_rcsample.get_solarsample,
                 define_rcsample.get_highfehsample]
    names= ['lowlow','highalpha','solar','highfeh']
    for ii in range(4):
        data= load_funcs[ii]()
        if ii == 1:
            ylabel= r'$Z\,(\mathrm{kpc})$'
        else:
            ylabel= ' '
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(data['RC_GALR_H'],
                            data['RC_GALZ_H'],
                            'k.',ms=2.,
                            xrange=[0.,16.],
                            yrange=[-3.,3.],
                            xlabel=r'$R\,(\mathrm{kpc})$',
                            ylabel=ylabel,
                            onedhists=True,
                            bins=31)
        bovy_plot.bovy_end_print(plotname.replace('SAMPLE',names[ii]))
    return None

if __name__ == '__main__':
    plot_spatial_broad(sys.argv[1])
