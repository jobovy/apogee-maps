###############################################################################
# plot_ah_location: plot the range of extinctions effor a given  location
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
from define_rcsample import get_rcsample
_PLOTDIST= True
_LW= 1.5
def plot_ah_location(location,plotname):
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
        save_pickles(selectFile,apo)
    glon, glat= apo.glonGlat(location)
    glon= glon[0]
    glat= glat[0]
    ahFile= '../savs/ah-%i.sav' % location
    if not os.path.exists(ahFile):
        # Distances at which to calculate the extinction
        distmods= numpy.linspace(7.,15.5,301)
        ds= 10.**(distmods/5-2.)
        # Setup Green et al. (2015) dust map
        gd= mwdust.Green15(filter='2MASS H')
        pa, ah= gd.dust_vals_disk(glon,glat,ds,apo.radius(location))
        meanah_default= numpy.sum(numpy.tile(pa,(len(ds),1)).T*ah,axis=0)/numpy.sum(pa)
        stdah_default= numpy.sqrt(numpy.sum(numpy.tile(pa,(len(ds),1)).T\
                                                *ah**2.,axis=0)\
                                      /numpy.sum(pa)-meanah_default**2.)
        # Marshall et al. (2006)
        marshall= mwdust.Marshall06(filter='2MASS H')
        try:
            pa, ah= marshall.dust_vals_disk(glon,glat,ds,apo.radius(location))
        except IndexError:
            meanah_marshall= -numpy.ones_like(ds)
            stdah_marshall= -numpy.ones_like(ds)
        else:
            meanah_marshall= numpy.sum(numpy.tile(pa,(len(ds),1)).T*ah,
                                       axis=0)/numpy.sum(pa)
            stdah_marshall= numpy.sqrt(numpy.sum(numpy.tile(pa,(len(ds),1)).T\
                                                     *ah**2.,axis=0)\
                                       /numpy.sum(pa)-meanah_marshall**2.)
        if True:
            # Drimmel et al. (2003)
            drimmel= mwdust.Drimmel03(filter='2MASS H')
            pa, ah= drimmel.dust_vals_disk(glon,glat,ds,apo.radius(location))
            meanah_drimmel= numpy.sum(numpy.tile(pa,(len(ds),1)).T*ah,axis=0)/numpy.sum(pa)
            stdah_drimmel= numpy.sqrt(numpy.sum(numpy.tile(pa,(len(ds),1)).T\
                                                    *ah**2.,axis=0)\
                                          /numpy.sum(pa)-meanah_drimmel**2.)
        else:
            meanah_drimmel= -numpy.ones_like(ds)
            stdah_drimmel= -numpy.ones_like(ds)
        if False:
            # Sale et al. (2014)
            sale= mwdust.Sale14(filter='2MASS H')
            try:
                pa, ah= sale.dust_vals_disk(glon,glat,ds,apo.radius(location))
            except (TypeError,ValueError):
                meanah_sale= -numpy.ones_like(ds)
                stdah_sale= -numpy.ones_like(ds)
            else:
                meanah_sale= numpy.sum(numpy.tile(pa,(len(ds),1)).T*ah,
                                       axis=0)/numpy.sum(pa)
                stdah_sale= numpy.sqrt(numpy.sum(numpy.tile(pa,(len(ds),1)).T\
                                                     *ah**2.,axis=0)\
                                           /numpy.sum(pa)-meanah_sale**2.)
        else:
                meanah_sale= -numpy.ones_like(ds)
                stdah_sale= -numpy.ones_like(ds)
        save_pickles(ahFile,distmods,meanah_default,stdah_default,
                     meanah_marshall,stdah_marshall,
                     meanah_drimmel,stdah_drimmel,
                     meanah_sale,stdah_sale)
    else:
        with open(ahFile,'rb') as savefile:
            distmods= pickle.load(savefile)
            meanah_default= pickle.load(savefile)
            stdah_default= pickle.load(savefile)
            meanah_marshall= pickle.load(savefile)
            stdah_marshall= pickle.load(savefile)
            meanah_drimmel= pickle.load(savefile)
            stdah_drimmel= pickle.load(savefile)
            meanah_sale= pickle.load(savefile)
            stdah_sale= pickle.load(savefile)
    # Now plot
    bovy_plot.bovy_print(fig_height=3.)
    if _PLOTDIST:
        distmods= 10.**(distmods/5-2.)
        xrange= [0.,12.]
        xlabel=r'$D\,(\mathrm{kpc})$'
    else:
        xrange=[7.,15.8],
        xlabel=r'$\mathrm{distance\ modulus}\ \mu$'
    ylabel=r'$A_H$'
    yrange= [0.,1.2*numpy.amax(numpy.vstack((meanah_default+stdah_default,
                                             meanah_marshall+stdah_marshall,
                                             meanah_drimmel+stdah_drimmel,
                                             meanah_sale+stdah_sale)))]
    line_default= bovy_plot.bovy_plot(distmods,meanah_default,
                                      'b-',lw=_LW,zorder=12,
                                      xrange=xrange,
                                      xlabel=xlabel,
                                      yrange=yrange,
                                      ylabel=ylabel)
    pyplot.fill_between(distmods,
                        meanah_default-stdah_default,
                        meanah_default+stdah_default,
                        hatch='//',facecolor=(0,0,0,0),
                        color='b',lw=0.25,zorder=4)
    line_marshall= bovy_plot.bovy_plot(distmods,meanah_marshall,'r-',lw=_LW,
                                       overplot=True,
                                       zorder=8)
    pyplot.fill_between(distmods,
                        meanah_marshall-stdah_marshall,
                        meanah_marshall+stdah_marshall,
                        hatch='\\\\',facecolor=(0,0,0,0),
                        color='r',lw=0.25,zorder=2)
    line_drimmel= bovy_plot.bovy_plot(distmods,meanah_drimmel,'-',lw=_LW,
                                      color='gold',
                                       overplot=True,
                                       zorder=7)
    pyplot.fill_between(distmods,
                        meanah_drimmel-stdah_drimmel,
                        meanah_drimmel+stdah_drimmel,
                        hatch='///',facecolor=(0,0,0,0),
                        color='gold',lw=0.25,zorder=1)
    line_sale= bovy_plot.bovy_plot(distmods,meanah_sale,'-',lw=_LW,
                                   color='c',
                                   overplot=True,
                                   zorder=9)
    pyplot.fill_between(distmods,
                        meanah_sale-stdah_sale,
                        meanah_sale+stdah_sale,
                        hatch='///',facecolor=(0,0,0,0),
                        color='c',lw=0.25,zorder=3)
    if True:
        data= get_rcsample()
        data= data[data['LOCATION_ID'] == location]
        bovy_plot.bovy_plot(data['RC_DIST'],data['AK_TARG']*1.55,
                            'ko',zorder=20,overplot=True,ms=2.)
    if location == 4378:
        pyplot.legend((line_default[0],line_jkz[0],line_zero[0]),
                      (r'$\mathrm{Green\ et\ al.\ (2015)}$',
                       r'$\mathrm{Green\ et\ al.} + p(M_H)$',
                       r'$\mathrm{zero\ extinction}$'),
                      loc='lower right',#bbox_to_anchor=(.91,.375),
                      numpoints=8,
                      prop={'size':14},
                      frameon=False)
    elif location == 4312:
        pyplot.legend((line_sale[0],line_marshall[0],line_drimmel[0]),
                      (r'$\mathrm{Sale\ et\ al.\ (2014)}$',
                       r'$\mathrm{Marshall\ et\ al.\ (2006)}$',
                       r'$\mathrm{Drimmel\ et\ al.\ (2003)}$'),
                      loc='lower right',#bbox_to_anchor=(.91,.375),
                      numpoints=8,
                      prop={'size':14},
                      frameon=False)                      
    # Label
    lcen, bcen= apo.glonGlat(location)
    if numpy.fabs(bcen) < 0.1: bcen= 0.
    bovy_plot.bovy_text(r'$(l,b) = (%.1f,%.1f)$' % (lcen,bcen),
                        top_right=True,size=16.)
    bovy_plot.bovy_end_print(plotname)
    return None

if __name__ == '__main__':
    #4240 is 30,0
    plot_ah_location(int(sys.argv[1]),sys.argv[2])
