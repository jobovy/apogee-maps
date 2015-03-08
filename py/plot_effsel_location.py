###############################################################################
# plot_effsel_location: plot the effective selection function for a given 
#                       location
###############################################################################
import os, os.path
import sys
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import save_pickles, bovy_plot
from matplotlib import rc, pyplot
import mwdust
import apogee.select.apogeeSelect
_PLOTDIST= True
def plot_effsel_location(location,plotname):
    # Setup selection function
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    else:
        # Setup selection function
        apo= apogee.select.apogeeSelect()
        # Delete these because they're big and we don't need them
        del apo._specdata
        del apo._photdata
    effselFile= '../savs/effselfunc-%i.sav' % location
    if not os.path.exists(effselFile):
        # Distances at which to calculate the effective selection function
        distmods= numpy.linspace(7.,15.5,101)
        ds= 10.**(distmods/5-2.)
        # Setup default effective selection function
        gd= mwdust.Green15(filter='2MASS H',load_samples=True)
        apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=gd)
        sf_default= apof(location,ds)
        sf_samples= numpy.zeros((20,len(ds)))
        if True:
            for ii in range(20):
                # Swap in a sample for bestfit in the Green et al. (2015) dmap
                gd._intps= numpy.zeros(len(gd._pix_info['healpix_index']),
                                       dtype='object') # need to remove the cache
                gd._best_fit= gd._samples[:,ii,:]
                apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=gd)
                sf_samples[ii]= apof(location,ds)          
        zerodust= mwdust.Zero(filter='2MASS H')
        apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=zerodust)
        sf_zero= apof(location,ds)
        drimmel= mwdust.Drimmel03(filter='2MASS H')
        apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=drimmel)
        sf_drimmel= apof(location,ds)
        marshall= mwdust.Marshall06(filter='2MASS H')
        apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=marshall)
        sf_marshall= apof(location,ds)
        save_pickles(effselFile,distmods,sf_default,sf_samples,sf_zero,
                     sf_drimmel,sf_marshall)
    else:
        with open(effselFile,'rb') as savefile:
            distmods= pickle.load(savefile)
            sf_default= pickle.load(savefile)
            sf_samples= pickle.load(savefile)
            sf_zero= pickle.load(savefile)
            sf_drimmel= pickle.load(savefile)
            sf_marshall= pickle.load(savefile)
    # Now plot
    bovy_plot.bovy_print(fig_height=3.)
    rc('text.latex', preamble=r'\usepackage{amsmath}'+'\n'
       +r'\usepackage{amssymb}'+'\n'+r'\usepackage{yfonts}')
    if _PLOTDIST:
        distmods= 10.**(distmods/5-2.)
        xrange= [0.,12.]
    else:
        xrange=[7.,15.8],
    bovy_plot.bovy_plot(distmods,sf_default,
                        'k-',zorder=10,
                        xrange=xrange,
                        yrange=[0.,1.1*numpy.amax(sf_zero)],
                        xlabel=r'$\mathrm{distance\ modulus}\ \mu$',
                        ylabel=r'$\textswab{S}(\mathrm{location},\mu)$')
    pyplot.fill_between(distmods,numpy.amin(sf_samples,axis=0),
                        numpy.amax(sf_samples,axis=0),color='0.4',
                        zorder=0)
    bovy_plot.bovy_plot(distmods,sf_zero,'b-',overplot=True,zorder=7)
    bovy_plot.bovy_plot(distmods,sf_drimmel,'y-',overplot=True,zorder=8)
    bovy_plot.bovy_plot(distmods,sf_marshall,'r-',overplot=True,zorder=9)
    bovy_plot.bovy_end_print(plotname)
    return None

if __name__ == '__main__':
    #4240 is 30,0
    plot_effsel_location(int(sys.argv[1]),sys.argv[2])
