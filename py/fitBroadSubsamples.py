###############################################################################
# fitBroadSubsamples.py: do and report fits for the broad sub-samples
###############################################################################
import sys
import os, os.path
import pickle
import copy
import numpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot, save_pickles
from matplotlib import pyplot
import define_rcsample
import fitDens
import densprofiles
import compareDataModel
_NSAMPLES= 50000
_LW= 1.5
# Globals
locations= None
distmods= None
effsel= None
effsel_mar= None
effsel_drim= None
effsel_sale= None
effsel_zero= None
lcen= None
bcen= None
highbIndx= None
outDiskIndx= None
betwDiskIndx= None 
inDiskIndx= None
ldata= None
data_highbIndx= None
data_outDiskIndx= None
data_betwDiskIndx= None
data_inDiskIndx= None
data_bulgeIndx= None
data_brightIndx= None
data_mediumIndx= None
data_faintIndx= None
def fitBroadSubsamples(sample,savename):
    # Setup the selection function
    setup_selection_function()
    # Load the data
    load_data(sample)
    # Do the fits
    global bf_exp, bf_brexp, bf_twoexp
    global ml_exp, ml_brexp, ml_twoexp
    global samples_exp, samples_brexp, samples_twoexp
    if os.path.exists(savename):
        with open(savename,'rb') as savefile:
            bf_exp= pickle.load(savefile)
            bf_brexp= pickle.load(savefile)
            bf_twoexp= pickle.load(savefile)
            ml_exp= pickle.load(savefile)
            ml_brexp= pickle.load(savefile)
            ml_twoexp= pickle.load(savefile)
            samples_exp= pickle.load(savefile)
            samples_brexp= pickle.load(savefile)
            samples_twoexp= pickle.load(savefile)
    else:
        # Perform fits of
        # a) expplusconst
        bf_exp, samples_exp, ml_exp= fit(type='expplusconst')
        # b) brokenexpflare
        bf_brexp, samples_brexp, ml_brexp= fit(type='tribrokenexpflare')
        # c) brokentwoexp
        bf_twoexp, samples_twoexp, ml_twoexp= fit(type='tribrokentwoexp')
        save_pickles(savename,
                     bf_exp,bf_brexp,bf_twoexp,
                     ml_exp,ml_brexp,ml_twoexp,
                     samples_exp,samples_brexp,samples_twoexp)
    # Do the rest of the fits as justfit
    global bf_brexp_g15, ml_brexp_g15, bf_brexp_drim, ml_brexp_drim
    global bf_brexp_sale, ml_brexp_sale, bf_brexp_zero, ml_brexp_zero
    bf_brexp_g15, ml_brexp_g15= fit(type='tribrokenexpflare',dmap='green15',
                                    justfit=True,init=bf_brexp) 
    bf_brexp_drim, ml_brexp_drim= fit(type='tribrokenexpflare',dmap='drimmel03',
                                      justfit=True,init=bf_brexp)
    bf_brexp_sale, ml_brexp_sale= fit(type='tribrokenexpflare',dmap='sale14',
                                      justfit=True,init=bf_brexp)
    bf_brexp_zero, ml_brexp_zero= fit(type='tribrokenexpflare',dmap='zero',
                                      justfit=True,init=bf_brexp)
    return None

def fit(type='tribrokenexpflare',dmap='marshall06',init=None,justfit=False):
    # Get the correct selection function
    tlocations= copy.deepcopy(locations)
    tdistmods= copy.deepcopy(distmods)
    if dmap == 'green15':
        teffsel= copy.deepcopy(effsel)
    elif dmap.lower() == 'marshall06':        
        teffsel= copy.deepcopy(effsel_mar)
    elif dmap.lower() == 'sale14':        
        teffsel= copy.deepcopy(effsel_sale)
    elif dmap.lower() == 'drimmel03':        
        teffsel= copy.deepcopy(effsel_drim)
    elif dmap.lower() == 'zero':        
        teffsel= copy.deepcopy(effsel_zero)
    # Now fit, sample, fit
    fitOut= fitDens.fitDens(ldata,
                            numpy.array(tlocations),
                            copy.deepcopy(teffsel),
                            tdistmods,
                            type=type,
                            init=init,
                            retMaxL=True,
                            nsamples=_NSAMPLES,mcmc=not justfit)
    if justfit: return fitOut
    bf, samples, maxl= fitOut
    bf, maxl= fitDens.fitDens(ldata,
                              numpy.array(tlocations),
                              copy.deepcopy(teffsel),
                              tdistmods,
                              type=type,
                              init=numpy.median(samples,axis=1),
                              retMaxL=True)
    return (bf,samples,maxl)

def setup_selection_function():
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    # Green et al. 
    global locations, distmods, effsel, effsel_mar, effsel_drim, effsel_sale
    global effsel_zero
    with open('../essf/essf_green15.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel= pickle.load(savefile)
        distmods= pickle.load(savefile)
    # Marshall et al. (2006)
    with open('../essf/essf_marshall06.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_mar= pickle.load(savefile)
    # Fill in regions not covered by Sale map
    effsel_mar[effsel_mar < -0.5]= effsel[effsel_mar < -0.5]
    # Sale et al. (2014)
    with open('../essf/essf_sale14.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_sale= pickle.load(savefile)
    # Fill in regions not covered by Marshall map
    effsel_sale[effsel_sale < -0.5]= effsel[effsel_sale < -0.5]
    # Drimmel et al (2003)
    with open('../essf/essf_drimmel03.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_drim= pickle.load(savefile)
    # Zero
    with open('../essf/essf_zero.sav','rb') as savefile:
        locations= pickle.load(savefile)
        effsel_zero= pickle.load(savefile)
    # Get (lcen,bcen) for each location
    global lcen, bcen
    lcen= numpy.zeros(len(locations))
    bcen= numpy.zeros(len(locations))
    for ii,loc in enumerate(locations):
        tlcen, tbcen= apo.glonGlat(loc)
        lcen[ii]= tlcen
        bcen[ii]= tbcen
    # Get the locations of various subsamples
    global highbIndx, outDiskIndx, betwDiskIndx, inDiskIndx
    highbIndx= numpy.fabs(bcen) > 10.
    outDiskIndx= (lcen > 140.)*(lcen < 250.)*(True-highbIndx)
    betwDiskIndx= (lcen <= 140.)*(lcen >= 70.)*(True-highbIndx)
    inDiskIndx= (lcen < 70.)*(True-highbIndx)
    return None

def load_data(sample):
    global ldata
    global data_highbIndx
    global data_outDiskIndx
    global data_betwDiskIndx
    global data_inDiskIndx    
    global data_bulgeIndx    
    global data_brightIndx
    global data_mediumIndx
    global data_faintIndx
    if sample.lower() == 'all':
        ldata= define_rcsample.get_rcsample()
    elif sample.lower() == 'alllowalpha':
        ldata= define_rcsample.get_rcsample()
        ldata= ldata[ldata[define_rcsample._AFETAG] < 0.1]
    elif sample.lower() == 'lowlow':
        ldata= define_rcsample.get_lowlowsample()
    elif sample.lower() == 'highfeh':
        ldata= define_rcsample.get_highfehsample()
    elif sample.lower() == 'highalpha':
        ldata= define_rcsample.get_highalphasample()
    elif sample.lower() == 'solar':
        ldata= define_rcsample.get_solarsample()
    # Get the indices of the various subsamples defined above
    data_highbIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[highbIndx]:
            data_highbIndx[ii]= True
    data_outDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[outDiskIndx]:
            data_outDiskIndx[ii]= True
    data_betwDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[betwDiskIndx]:
            data_betwDiskIndx[ii]= True
    data_inDiskIndx= numpy.zeros(len(ldata),dtype='bool')
    for ii in range(len(ldata)):
        if ldata[ii]['LOCATION_ID'] in numpy.array(locations)[inDiskIndx]:
            data_inDiskIndx[ii]= True
    return None

def writeTable(sample,savename,tablename):
    delimiter= ' & '
    types= ['tribrokenexpflare','expplusconst','brokentwoexp']
    densmodels= ['broken exp. w/ flare','single exp.',
                 'broken exp. w/ 2 $h_Z$']
    extmaps= ['\citet{Marshall06a}',
              '\citet{Green15a}',
              '\citet{Sale14a}',
              '\citet{Drimmel03a}',
              'zero']
    with open(tablename,'w') as tablefile:
        # Start w/ fiducial fit
        if sample.lower() == 'lowlow':
            printline= 'low [Fe/H]'
        elif sample.lower() == 'solar':
            printline= 'solar'
        elif sample.lower() == 'highfeh':
            printline= 'high [Fe/H]'
        elif sample.lower() == 'highalpha':
            printline= 'high [$\\alpha$/Fe]'
        printline+= delimiter
        # Fiducial
        printline+= densmodels[0]+delimiter
        printline+= extmaps[0]+delimiter
        printline+= _format_results(types[0],extmaps[0])
        tablefile.write(printline+'\\\\\n')
        # Green
        printline= delimiter+delimiter+extmaps[1]+delimiter
        printline+= _format_results_noerr(types[0],extmaps[1])       
        tablefile.write(printline+'\\\\\n')
        # Sale
        printline= delimiter+delimiter+extmaps[2]+delimiter
        printline+= _format_results_noerr(types[0],extmaps[2])       
        tablefile.write(printline+'\\\\\n')
        # Drimmel
        printline= delimiter+delimiter+extmaps[3]+delimiter
        printline+= _format_results_noerr(types[0],extmaps[3])       
        tablefile.write(printline+'\\\\\n')
        # Zero
        printline= delimiter+delimiter+extmaps[4]+delimiter
        printline+= _format_results_noerr(types[0],extmaps[4])       
        tablefile.write(printline+'\\\\\n')
        # Exp.
        printline= delimiter+densmodels[1]+delimiter+extmaps[0]+delimiter
        printline+= _format_results(types[1],extmaps[0])       
        tablefile.write(printline+'\\\\\n')
        # twoexp
        printline= delimiter+densmodels[2]+delimiter+extmaps[0]+delimiter
        printline+= _format_results(types[2],extmaps[0])       
        tablefile.write(printline+'\\\\\n')
    return None

def _format_results(type,extmap):
    if type.lower() == 'tribrokenexpflare':
        tsamples= samples_brexp
        out= {'hr1':numpy.median(1./tsamples[0]),
              'hr1err':numpy.std(1./tsamples[0]),
              'hr2':numpy.median(1./tsamples[2]),
              'hr2err':numpy.std(1./tsamples[2]),
              'rmax':numpy.median(numpy.exp(tsamples[3])),
              'rmaxerr':numpy.std(numpy.exp(tsamples[3])),
              'hz':numpy.median(1./tsamples[1]),
              'hzerr':numpy.std(1./tsamples[1]),
              'rf':numpy.median(tsamples[4]),
              'rferr':numpy.std(tsamples[4]),
              'ml':0.}
        if out['hr1'] > 10.:
            out['hr1']= sorted(1./tsamples[0])[int(round(0.05*tsamples.shape[1]))]
            return "$>{hr1:.1f}$&${hr2:.1f}\pm{hr2err:.1f}$&${rmax:.1f}\pm{rmaxerr:.1f}$&${hz:.2f}\pm{hzerr:.2f}$&${rf:.2f}\pm{rferr:.2f}$&\ldots&{ml:.0f}".format(**out)
        else:
            return "${hr1:.1f}\pm{hr1err:.1f}$&${hr2:.1f}\pm{hr2err:.1f}$&${rmax:.1f}\pm{rmaxerr:.1f}$&${hz:.2f}\pm{hzerr:.2f}$&${rf:.2f}\pm{rferr:.2f}$&\ldots&{ml:.0f}".format(**out)
    elif type.lower() == 'expplusconst':
        tsamples= samples_exp
        out= {'hr':numpy.median(1./tsamples[0]),
              'hrerr':numpy.std(1./tsamples[0]),
              'hz':numpy.median(1./tsamples[1]),
              'hzerr':numpy.std(1./tsamples[1]),
              'ml':-2.*(ml_exp-ml_brexp)}
        if out['hr'] > 10.:
            out['hr']= sorted(1./tsamples[0])[int(round(0.05*tsamples.shape[1]))]
            return "$>{hr:.1f}$&\ldots&\ldots&${hz:.2f}\pm{hzerr:.2f}$&\ldots&\ldots&{ml:.0f}".format(**out)
        else:
            return "${hr:.1f}\pm{hrerr:.1f}$&\ldots&\ldots&${hz:.2f}\pm{hzerr:.2f}$&\ldots&&{ml:.0f}".format(**out)
    if type.lower() == 'tribrokentwoexp':
        tsamples= samples_twoexp
        out= {'hr1':numpy.median(1./tsamples[0]),
              'hr1err':numpy.std(1./tsamples[0]),
              'hr2':numpy.median(1./tsamples[2]),
              'hr2err':numpy.std(1./tsamples[2]),
              'rmax':numpy.median(numpy.exp(tsamples[3])),
              'rmaxerr':numpy.std(numpy.exp(tsamples[3])),
              'hz1':numpy.median(1./tsamples[1]),
              'hz1err':numpy.std(1./tsamples[1]),
              'amp':numpy.median(densprofiles.ilogit(tsamples[4])),
              'amperr':numpy.std(densprofiles.ilogit(tsamples[4])),
              'hz2':numpy.median(1./tsamples[5]),
              'hz2err':numpy.std(1./tsamples[5]),
              'ml':-2.*(ml_twoexp-ml_brexp)}
        return "${hr1:.1f}\pm{hr1err:.1f}$&${hr2:.1f}\pm{hr2err:.1f}$&${rmax:.1f}\pm{rmaxerr:.1f}$&${hz1:.2f}\pm{hz1err:.2f}$&${amp:.2f}\pm{amperr:.2f}$&${hz2:.2f}\pm{hz2err:.2f}$&{ml:.0f}".format(**out)
                
def _format_results_noerr(type,extmap):
    if type.lower() == 'tribrokenexpflare':
        if 'Mar' in extmap:
            tbf= bf_brexp
            tml= ml_brexp
        elif 'Gre' in extmap:
            tbf= bf_brexp_g15
            tml= ml_brexp_g15
        elif 'Dri' in extmap:
            tbf= bf_brexp_drim
            tml= ml_brexp_drim
        elif 'Sal' in extmap:
            tbf= bf_brexp_sale
            tml= ml_brexp_sale
        elif 'zero' in extmap:
            tbf= bf_brexp_zero
            tml= ml_brexp_zero
        out= {'hr1':numpy.median(1./tbf[0]),
              'hr2':numpy.median(1./tbf[2]),
              'rmax':numpy.median(numpy.exp(tbf[3])),
              'hz':numpy.median(1./tbf[1]),
              'rf':numpy.median(tbf[4]),
              'ml':-2.*(tml-ml_brexp)}
        return "{hr1:.1f}&{hr2:.1f}&{rmax:.1f}&{hz:.2f}&{rf:.2f}&\ldots&{ml:.0f}".format(**out)
                
def plotCompareData(sample,savename,plotname):
    locs= ['highb','outdisk','meddisk','indisk']
    indices= [highbIndx,outDiskIndx,betwDiskIndx,inDiskIndx]
    data_indices= [data_highbIndx,data_outDiskIndx,
                   data_betwDiskIndx,data_inDiskIndx]
    # Full prediction for numbers
    Xs,pd= compareDataModel.predict_spacedist(bf_brexp,
                                              numpy.array(locations),
                                              copy.deepcopy(effsel_mar),
                                              distmods,
                                              type='tribrokenexpflare',coord='dm')
    for loc, index, data_index in zip(locs,indices,data_indices):
        bovy_plot.bovy_print()
        # High |b|
        # Marshall is fiducial
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_mar)[index],
                                                   distmods,type='tribrokenexpflare',
                                                   coord='dm')
        bovy_plot.bovy_hist(ldata['RC_DM_H'][data_index],
                            histtype='step',normed=True,
                            lw=_LW,
                            range=[7.,15.5],
                            bins=round(numpy.sqrt(numpy.sum(data_index))*2.),
                            yrange=[0.,
                                    1.2*numpy.amax(pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]))],
                            color='k',
                            xlabel=r'$\mu$')
        line_mar= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                      color='r',
                                      lw=2.*_LW,overplot=True,zorder=12)
        bovy_plot.bovy_text(r'$%i / %i\ \mathrm{stars}$' \
                                % (numpy.sum(data_index),
                                   len(ldata))
                            +'\n'+r'$%i / %i\ \mathrm{predicted}$' \
                                % (numpy.sum(pdt)/numpy.sum(pd)*len(ldata),len(ldata)),
                            top_left=True,size=14.)
        # Green
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_g15,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel)[index],
                                                   distmods,type='tribrokenexpflare',
                                                   coord='dm')
        line_g15= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                      color='b',
                                      lw=_LW,overplot=True,zorder=13)
        # Drimmel
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_drim,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_drim)[index],
                                                   distmods,type='tribrokenexpflare',
                                                   coord='dm')
        line_drim= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                       color='gold',
                                       lw=_LW,overplot=True,zorder=12)
        # Sale
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_sale,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_sale)[index],
                                                   distmods,type='tribrokenexpflare',
                                                   coord='dm')
        line_sale= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                       color='c',
                                       lw=_LW,overplot=True,zorder=12)
        # Zero
        Xs,pdt= compareDataModel.predict_spacedist(bf_brexp_zero,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_zero)[index],
                                                   distmods,type='tribrokenexpflare',
                                                   coord='dm')
        line_zero= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                       color='k',
                                       ls='-',lw=_LW,overplot=True,zorder=10)
        # Marshall + exp
        Xs,pdt= compareDataModel.predict_spacedist(bf_exp,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_mar)[index],
                                                   distmods,type='expplusconst',
                                                   coord='dm')
        line_exp= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                      color='r',
                                      lw=2*_LW,overplot=True,zorder=10,ls=':')
        # Marshall + twoexp
        Xs,pdt= compareDataModel.predict_spacedist(bf_twoexp,
                                                   numpy.array(locations)[index],
                                                   copy.deepcopy(effsel_mar)[index],
                                                   distmods,type='tribrokentwoexp',
                                                   coord='dm')
        line_twoexp= bovy_plot.bovy_plot(Xs,pdt/numpy.sum(pdt)/(Xs[1]-Xs[0]),
                                         color='r',
                                         lw=2*_LW,overplot=True,zorder=11,
                                         ls='--')
        if sample.lower() == 'lowlow' or sample.lower() == 'highalpha':
            if loc.lower() == 'highb':
                bovy_plot.bovy_text(r'$|b| > 10^\circ$',title=True,size=14.)
            elif loc.lower() == 'indisk':
                bovy_plot.bovy_text(r'$l < 70^\circ, |b| \leq 10^\circ$',
                                    title=True,size=14.)
            elif loc.lower() == 'meddisk':
                bovy_plot.bovy_text(r'$70^\circ \leq l \leq 140^\circ, |b| \leq 10^\circ$',
                                    title=True,size=14.)
            elif loc.lower() == 'outdisk':
                bovy_plot.bovy_text(r'$140^\circ < l < 250^\circ, |b| \leq 10^\circ$',
                                    title=True,size=14.)
            # Legend
            if loc.lower() == 'meddisk':
                pyplot.legend((line_mar[0],line_exp[0],line_twoexp[0]),
                              (r'$\mathrm{Marshall\ et\ al.\ (2006)}$'
                               +'\n'+r'$\mathrm{broken\ exp.\ w/\ flare}$',
                               r'$\mathrm{single\ exp.}$',
                               r'$\mathrm{broken\ exp.\ w/\ 2}\ h_Z$'),
                              loc='lower right',bbox_to_anchor=(.66,.42),
                              numpoints=8,
                              prop={'size':14},
                              frameon=False)
            elif loc.lower() == 'outdisk':
                pyplot.legend((line_g15[0],line_sale[0],line_drim[0],
                               line_zero[0]),
                              (r'$\mathrm{Green\ et\ al.\ (2015)}$',
                               r'$\mathrm{Sale\ et\ al.\ (2014)}$',
                               r'$\mathrm{Drimmel\ et\ al.\ (2003)}$',
                               r'$\mathrm{zero\ extinction}$'),
                              loc='lower right',bbox_to_anchor=(.66,.42),
                              numpoints=8,
                              prop={'size':14},
                              frameon=False)
        if loc.lower() == 'highb':
            if sample.lower() == 'lowlow':
                bovy_plot.bovy_text(r'$\mathrm{low\ [Fe/H]}$',
                                    top_right=True,size=18.)
            elif sample.lower() == 'solar':
                bovy_plot.bovy_text(r'$\mathrm{solar}$',
                                    top_right=True,size=18.)
            elif sample.lower() == 'highfeh':
                bovy_plot.bovy_text(r'$\mathrm{high\ [Fe/H]}$',
                                    top_right=True,size=18.)
            elif sample.lower() == 'highalpha':
                bovy_plot.bovy_text(r'$\mathrm{high}\ [\alpha/\mathrm{Fe]}$',
                                    top_right=True,size=18.)
        bovy_plot.bovy_end_print(plotname.replace('LOC',loc))
    return None

if __name__ == '__main__':
    # Input:
    #  - sample: 'lowlow', 'solar', 'highfeh', 'highalpha'
    #  - savename: name of the file for the pickle
    #  - tablename: name of the file for the table
    #  - plotname: name of the file for the plot
    sample= sys.argv[1]
    savename= sys.argv[2]
    tablename= sys.argv[3]
    plotname= sys.argv[4]
    # First, do the fits
    fitBroadSubsamples(sample,savename)
    # Then write the table
    writeTable(sample,savename,tablename)
    # And make the plot comparing data and model
    plotCompareData(sample,savename,plotname)
