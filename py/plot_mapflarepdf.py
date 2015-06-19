###############################################################################
# plot_mapflarepdf.py: plot and combine the PDFs for the flare for low alpha 
#                      MAPs
###############################################################################
import sys
import os, os.path
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot, save_pickles
from matplotlib import pyplot, cm
import define_rcsample
import extreme_deconvolution as XD
def plot_mapflarepdf(savename,plotname):
    # Load the samples
    with open('../mapfits/tribrokenexpflare.sav','rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples= numpy.array(pickle.load(savefile))
    maps= define_rcsample.MAPs()
    # Loop through the low-alpha MAPs and compute the XD decomposition
    if 'lowalpha' in savename:
        plotmaps= [9,16,23,29,36,43,50,57,64,71]
    else:
        plotmaps= [19,26,32,39,45]
    if not os.path.exists(savename):
        ngauss= 2
        allxamp= numpy.empty((len(plotmaps),ngauss))
        allxmean= numpy.empty((len(plotmaps),ngauss,1))
        allxcovar= numpy.empty((len(plotmaps),ngauss,1,1))
        cnt= 0
        for ii, map in enumerate(maps.map()):
            if not ii in plotmaps: continue
            print ii
            # Fit PDFs with XD
            xamp= numpy.array([0.45,0.5])
            xmean= numpy.array([numpy.mean(samples[ii,4])
                                +numpy.random.normal()*numpy.std(samples[ii,4]),
                                numpy.mean(samples[ii,4])
                                +numpy.random.normal()*numpy.std(samples[ii,4])])[:,numpy.newaxis]
            xcovar= numpy.reshape(numpy.tile(numpy.var(samples[ii,4]),(2,1)),
                                  (2,1,1))
            XD.extreme_deconvolution(samples[ii,4][:,numpy.newaxis],
                                     numpy.zeros((len(samples[ii,4]),1)),
                                     xamp,xmean,xcovar)
            allxamp[cnt]= xamp
            allxmean[cnt]= xmean
            allxcovar[cnt]= xcovar
            cnt+= 1
        save_pickles(savename,allxamp,allxmean,allxcovar)
    else:
        with open(savename,'rb') as savefile:
            allxamp= pickle.load(savefile)
            allxmean= pickle.load(savefile)
            allxcovar= pickle.load(savefile)
    # Now plot
    cmap= cm.coolwarm
    xrange= [-0.37,0.25]
    if 'lowalpha' in savename:
#        xrange= [-0.4,0.2]
        yrange= [0.,30.]
        combDiv= 2.
        colorFunc= lambda x: cmap((x+0.5)*0.95/0.9+0.05)
    else:
#        xrange= [-0.3,0.3]
        yrange= [0.,13.5]
        colorFunc= lambda x: cmap((x+0.4)*0.95/0.5+0.05)
        combDiv= 1.5
    overplot= False
    plotXDFit= True
    cnt= 0
    bovy_plot.bovy_print(axes_labelsize=18,text_fontsize=18,
                         xtick_labelsize=14,ytick_labelsize=14)
    for ii, map in enumerate(maps.map()):
        if not ii in plotmaps: continue
        tfeh= round(numpy.median(map['FE_H'])*20.)/20.
        if tfeh == 0.35: tfeh= 0.4
        if tfeh == -0.0: tfeh= 0.0
        bovy_plot.bovy_hist(samples[ii,4],
                            range=xrange,bins=51,overplot=overplot,
                            yrange=yrange,
                            histtype='step',normed=True,zorder=2,
                            color=colorFunc(tfeh),
                            xlabel=r'$R_{\mathrm{flare}}^{-1}\,(\mathrm{kpc}^{-1})$')
        if plotXDFit:
            txs= numpy.linspace(xrange[0],xrange[1],1001)
            pyplot.plot(txs,1./numpy.sqrt(2.*numpy.pi)*(allxamp[cnt,0]/numpy.sqrt(allxcovar[cnt,0,0,0])*numpy.exp(-0.5*(txs-allxmean[cnt,0,0])**2./allxcovar[cnt,0,0,0])
                                                 +allxamp[cnt,1]/numpy.sqrt(allxcovar[cnt,1,0,0])*numpy.exp(-0.5*(txs-allxmean[cnt,1,0])**2./allxcovar[cnt,1,0,0])),
                        color=colorFunc(tfeh),
                        zorder=1)
        overplot=True
        cnt+= 1
    txs= numpy.linspace(xrange[0],xrange[1],1001)
    comb= numpy.ones_like(txs)
    for ii in range(len(plotmaps)):
        comb*= 1./numpy.sqrt(2.*numpy.pi)*(allxamp[ii,0]/numpy.sqrt(allxcovar[ii,0,0,0])*numpy.exp(-0.5*(txs-allxmean[ii,0,0])**2./allxcovar[ii,0,0,0])
                                           +allxamp[ii,1]/numpy.sqrt(allxcovar[ii,1,0,0])*numpy.exp(-0.5*(txs-allxmean[ii,1,0])**2./allxcovar[ii,1,0,0]))
    comb/= numpy.sum(comb)*(txs[1]-txs[0])
    pyplot.plot(txs,comb/combDiv,'k-',lw=2.,zorder=20)
    pyplot.plot([0.,0.],[0.,50.],'k--',lw=1.5,zorder=0)
    t= pyplot.text(xrange[0]+0.25*(xrange[1]-xrange[0])+0.03*('highalpha' in savename),
                        0.8*yrange[1],
                        r'$R_{\mathrm{flare}}^{-1} = %.2f \pm %.2f\,\mathrm{kpc}^{-1}$' % (numpy.sum(comb*txs)/numpy.sum(comb), numpy.sqrt(numpy.sum(comb*txs**2.)/numpy.sum(comb)-(numpy.sum(comb*txs)/numpy.sum(comb))**2.)),
                        size=18.)
    t.set_bbox(dict(color='w',edgecolor='none'))
    if 'lowalpha' in savename:
        bovy_plot.bovy_text(r'$\mathrm{low-}[\alpha/\mathrm{Fe}]\ \mathrm{MAPs}$',
                            top_left=True,
                            size=16.)
    else:
        bovy_plot.bovy_text(r'$\mathrm{high-}[\alpha/\mathrm{Fe}]\ \mathrm{MAPs}$',
                            top_left=True,
                            size=16.)
    bovy_plot.bovy_end_print(plotname)
    return None

if __name__ == '__main__':
    plot_mapflarepdf(sys.argv[1],sys.argv[2])
