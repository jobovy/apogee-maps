###############################################################################
# plot_afe_spectra.py: plot residual spectra for a [Fe/H] bin to demonstrate 
#                      that we can measure [a/Fe]
###############################################################################
import sys
import os, os.path
import pickle
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
from galpy.util import bovy_plot, save_pickles
import empca
import apogee.tools.read as apread
from apogee.tools import paramIndx
import apogee.spec.plot as splot
import apogee.spec.window as apwindow
import apogee.spec.cannon as apcannon
import apogee.spec.stack as apstack
import define_rcsample
_ERASESTR= "                                                                                "
def plot_afe_spectra(savename,plotname):
    # Load the data
    data= define_rcsample.get_rcsample()
    data= data[data['SNR'] > 200.]
    fehindx= (data['FE_H'] <= -0.25)*(data['FE_H'] > -0.35)
    fehdata= data[fehindx]
    # First compute the residuals and do the EM-PCA smoothing
    if not os.path.exists(savename):
        nspec= len(fehdata)
        spec= numpy.zeros((nspec,7214))
        specerr= numpy.zeros((nspec,7214))
        for ii in range(nspec):
            sys.stdout.write('\r'+"Loading spectrum %i / %i ...\r" % (ii+1,nspec))
            sys.stdout.flush()
            spec[ii]= apread.aspcapStar(fehdata['LOCATION_ID'][ii],
                                        fehdata['APOGEE_ID'][ii],
                                        ext=1,header=False,aspcapWavegrid=True)
            specerr[ii]= apread.aspcapStar(fehdata['LOCATION_ID'][ii],
                                           fehdata['APOGEE_ID'][ii],
                                           ext=2,header=False,
                                           aspcapWavegrid=True)
        teffs= fehdata['FPARAM'][:,paramIndx('teff')]
        loggs= fehdata['FPARAM'][:,paramIndx('logg')]
        metals= fehdata[define_rcsample._FEHTAG]
        cf, s, r= apcannon.quadfit(spec,specerr,
                                   teffs-4800.,loggs-2.85,metals+0.3,
                                   return_residuals=True)
        pr= numpy.zeros_like(r)
        # Deal w/ bad data
        _MAXERR= 0.02
        npca= 8
        pca_input= r
        pca_weights= (1./specerr**2.)
        pca_weights[pca_weights < 1./_MAXERR**2.]= 0. 
        nanIndx= numpy.isnan(pca_input) + numpy.isnan(pca_weights)
        pca_weights[nanIndx]= 0. 
        pca_input[nanIndx]= 0.
        # Run EM-PCA
        m= empca.empca(pca_input,pca_weights,nvec=npca,niter=25)#,silent=False)
        for jj in range(nspec):
            for kk in range(npca):
                pr[jj]+= m.coeff[jj,kk]*m.eigvec[kk]
        save_pickles(savename,pr,r,cf)
    else:
        with open(savename,'rb') as savefile:
            pr= pickle.load(savefile)
    # Now plot the various elements
    colormap= cm.seismic
    colorFunc= lambda afe: afe/0.25
    widths= [3.5,3.]
    yranges= [[-0.05,0.02],[-0.03,0.01]]
    for ee, elem in enumerate(['S','Ca1']):
        for ii in range(5):
            tindx= (fehdata[define_rcsample._AFETAG] > ii*0.05-0.025)\
                *(fehdata[define_rcsample._AFETAG] <= (ii+1)*0.05-0.025)
            args= (apstack.median(pr[tindx][:12]),elem,)
            kwargs= {'markLines':ii==4,
                     'yrange':yranges[ee],
                     'ylabel':'',
                     'cleanZero':False,
                     'zorder':int(numpy.floor(numpy.random.uniform()*5)),
                     'color':colormap(colorFunc(ii*0.05)),
                     'overplot':ii>0,
                     'fig_width':widths[ee]}
            if ii>0: kwargs.pop('fig_width')
            splot.windows(*args,**kwargs)
        bovy_plot.bovy_end_print(plotname.replace('ELEM',
                                                  elem.lower().capitalize()))
    # Also do Mg
    for ii in range(5):
        tindx= (fehdata[define_rcsample._AFETAG] > ii*0.05-0.025)\
            *(fehdata[define_rcsample._AFETAG] <= (ii+1)*0.05-0.025)
        args= (apstack.median(pr[tindx][:12]),)
        kwargs={'startindxs':[3012,3120,3990],
                'endindxs':[3083,3158,4012],
                'yrange':[-0.05,0.02],
                'ylabel':'',
                'cleanZero':False,
                '_markwav':[15745.017,15753.189,15770.055,15958.836],
                'zorder':int(numpy.floor(numpy.random.uniform()*5)),
                'color':colormap(colorFunc(ii*0.05)),
                'overplot':ii>0,
                'fig_width':5.,
                'markLines':True}
        if ii>0: kwargs.pop('fig_width')
        if ii != 4: 
            kwargs.pop('_markwav')
            kwargs.pop('markLines')
        kwargs['_startendskip']= 0
        kwargs['_noxticks']= True
        kwargs['_labelwav']= True
        splot.waveregions(*args,**kwargs)                          
        bovy_plot.bovy_text(r'$\mathrm{Mg}$',
                            top_left=True,fontsize=10,backgroundcolor='w')
    bovy_plot.bovy_end_print(plotname.replace('ELEM','Mg'))
    # Also do Si
    for ii in range(5):
        tindx= (fehdata[define_rcsample._AFETAG] > ii*0.05-0.025)\
            *(fehdata[define_rcsample._AFETAG] <= (ii+1)*0.05-0.025)
        args= (apstack.median(pr[tindx][:12]),)
        kwargs={'startindxs':[4469, 4624,5171, 7205, 7843],
                'endindxs':[4488, 4644,5182, 7243, 7871],
                'yrange':[-0.05,0.02],
                'ylabel':'',
                'cleanZero':False,
                '_markwav':apwindow.lines('Si'),
                'zorder':int(numpy.floor(numpy.random.uniform()*5)),
                'color':colormap(colorFunc(ii*0.05)),
                'overplot':ii>0,
                'fig_width':6.,
                'markLines':True}
        if ii>0: kwargs.pop('fig_width')
        if ii != 4: 
            kwargs.pop('_markwav')
            kwargs.pop('markLines')
        kwargs['_startendskip']= 0
        kwargs['_noxticks']= True
        kwargs['_labelwav']= True
        splot.waveregions(*args,**kwargs)                          
        bovy_plot.bovy_text(r'$\mathrm{Si}$',
                            top_left=True,fontsize=10,backgroundcolor='w')
    bovy_plot.bovy_end_print(plotname.replace('ELEM','Si2'))
    # Also do Oxygen
    for ii in range(5):
        tindx= (fehdata[define_rcsample._AFETAG] > ii*0.05-0.025)\
            *(fehdata[define_rcsample._AFETAG] <= (ii+1)*0.05-0.025)
        args= (apstack.median(pr[tindx][:12]),)
        kwargs={'startlams':[15558,16242,16536,16720],
                'endlams':[15566,16250,16544,16728],
                'yrange':[-0.05,0.02],
                'ylabel':'',
                'cleanZero':False,
                '_markwav':[15562,16246,16539,16723.5],
                'zorder':int(numpy.floor(numpy.random.uniform()*5)),
                'color':colormap(colorFunc(ii*0.05)),
                'overplot':ii>0,
                'fig_width':4.5,
                'markLines':True}
        if ii>0: kwargs.pop('fig_width')
        if ii != 4: 
            kwargs.pop('_markwav')
            kwargs.pop('markLines')
        kwargs['_startendskip']= 0
        kwargs['_noxticks']= True
        kwargs['_labelwav']= True
        splot.waveregions(*args,**kwargs)                          
        bovy_plot.bovy_text(r'$\mathrm{O}$',
                            top_left=True,fontsize=10,backgroundcolor='w')
    bovy_plot.bovy_end_print(plotname.replace('ELEM','O'))
    return None

if __name__ == '__main__':
    numpy.random.seed(1)
    plot_afe_spectra(sys.argv[1],sys.argv[2])
