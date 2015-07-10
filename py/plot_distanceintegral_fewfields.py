###############################################################################
# plot_distanceintegral_fewfields.py: make a useful plot of the distance 
#                                     integration in the effective volume
#                                     for a few APOGEE-sized fields
###############################################################################
import sys
import os, os.path
import pickle
import numpy
from scipy import signal, integrate, interpolate
import healpy
import matplotlib
matplotlib.use('Agg')
from galpy.util import save_pickles, bovy_plot
import densprofiles
import dust
import gaia_rc
# magnitude limits for the survey
_GMIN= 3.
_GMAX= 20.
_NGSAMPLES= 10000
_NSIDE= 2048
_DEGTORAD= numpy.pi/180.
def plot_distanceintegral_fewfields(savename,plotname,apogee=False):
    #ls and bs of the patches
    ls= [280.,30.,12.,37.5,90.,127.5,30.,60.,20.,-20.,-40.]
    bs= [45.,4.,2.,0.,0.,0.,0.,0.,0.,0.,0.]
    radius= 1.49
    locs= [4214,4244,4378,4312,4242,4318,4240,4241] # for apogee
    # If we are doing APOGEE selection, load the APOGEE selection function
    if apogee:
        selectFile= '../savs/selfunc-nospdata.sav'
        if os.path.exists(selectFile):
            with open(selectFile,'rb') as savefile:
                apo= pickle.load(savefile)
    if os.path.exists(savename):
        with open(savename,'rb') as savefile:
            area= pickle.load(savefile)
            if apogee:
                area_cdens= pickle.load(savefile)
    else: 
        if not apogee:
            # For samping over the absolute magnitude distribution
            iso= gaia_rc.load_iso()
            Gsamples= gaia_rc.sample_Gdist(iso,n=_NGSAMPLES)
        # Loop through each distance slice
        area= numpy.zeros((len(ls),len(dust._GREEN15DISTS)))
        if apogee:
            area_cdens= numpy.zeros((len(ls),len(dust._GREEN15DISTS)))
        for ii,dist in enumerate(dust._GREEN15DISTS):
            print "Working on distance %i: %f" % (ii,dist)
            G0= 0.68+dust.dist2distmod(dist)
            H0= -1.49+dust.dist2distmod(dist)
            # Load the dust map
            combinedmap= dust.load_combined(dist,nest=False,nside_out=_NSIDE)
            # Loop through a few patches
            theta, phi= healpy.pixelfunc.pix2ang(_NSIDE,
                                                 numpy.arange(healpy.pixelfunc.nside2npix(_NSIDE)),
                                                 nest=False)
            cosb= numpy.sin(theta)
            for ff,(ll,bb) in enumerate(zip(ls,bs)):
                vec= healpy.pixelfunc.ang2vec((90.-bb)*_DEGTORAD,ll*_DEGTORAD)
                ipixs= healpy.query_disc(_NSIDE,vec,radius*_DEGTORAD,
                                         inclusive=False,nest=False)
                ipixIndx= numpy.in1d(numpy.arange(healpy.pixelfunc.nside2npix(_NSIDE)),
                                     ipixs,assume_unique=True)
                # Calculate the density
                densmap= densprofiles.expdisk(phi[ipixIndx],
                                              numpy.pi/2.-theta[ipixIndx],
                                              dist*numpy.ones(len(ipixs)),
                                              glon=True,
                                              params=[1./3.,1./0.3],nest=False)
                combinedmask= numpy.zeros(len(ipixs))
                tcombinedmap= combinedmap[ipixIndx]
                if apogee:
                    # Loop through the three cohorts
                    for c in ['short','medium','long']:
                        hmin= apo.Hmin(locs[ff],cohort=c)
                        hmax= apo.Hmax(locs[ff],cohort=c)
                        combinedmask+= ((tcombinedmap > (hmin-H0))\
                                            *(tcombinedmap < (hmax-H0))).astype('float')*apo(locs[ff],hmin+0.01)
                else:
                    for jj in range(_NGSAMPLES):
                        combinedmask+= ((tcombinedmap > (_GMIN-G0-Gsamples[jj]+0.68))\
                                            *(tcombinedmap < (_GMAX-G0-Gsamples[jj]+0.68))).astype('float')
                combinedmask/= _NGSAMPLES
                # Compute cross correlation
                area[ff,ii]= float(numpy.sum(cosb[ipixIndx]*densmap\
                                                 *combinedmask))
                if apogee:
                    # Also store the approximation that the density is constant
                    cdens= densprofiles.expdisk(ll*_DEGTORAD,bb*_DEGTORAD,dist,
                                                glon=True,
                                                params=[1./3.,1./0.3],
                                                nest=False)
                    area_cdens[ff,ii]= float(numpy.sum(cosb[ipixIndx]*cdens\
                                                           *combinedmask))
        if apogee:
            save_pickles(savename,area,area_cdens)
        else:
            save_pickles(savename,area)
    # Plot the power spectrum
    if True:
        bovy_plot.bovy_print(fig_height=3.)
        matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
        colors=['b','g','r','c','gold','m','k','orange','y','b','g']
        lss=['-','-','-','-','-','-','-','-','-','--','--']
        for ff in range(len(ls)):
            psdx, psd= signal.periodogram(area[ff]*dust._GREEN15DISTS**3./numpy.sum(area[ff]*dust._GREEN15DISTS**3.),
                                          fs=1./(dust._GREEN15DISTMODS[1]-dust._GREEN15DISTMODS[0]),
                                          detrend=lambda x: x,scaling='spectrum')
            bovy_plot.bovy_plot(psdx[1:],psd[1:],
                                '-',loglog=True,
                                xlabel=r'$2\pi\,k_\mu\,(\mathrm{mag}^{-1})$',
                                ylabel=r'$P_k$',
                                xrange=[0.04,4.],
                                color=colors[ff],ls=lss[ff],
                                overplot=ff > 0)
            if ff == 0:
                bovy_plot.bovy_text(r'$\mathrm{normalized}\ D^3\,\nu_*(\mu|\theta)\,\textswab{S}(\mu)$',
                                    bottom_left=True,size=16.)
        bovy_plot.bovy_end_print(plotname)               
    else:
        bovy_plot.bovy_print(fig_height=3.)
        matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{yfonts}"]
        for ff in range(len(ls)):
            bovy_plot.bovy_plot(dust._GREEN15DISTMODS,
                                area[ff]*dust._GREEN15DISTS**3.,
                                'k-',
                                xlabel=r'$\mu\,(\mathrm{mag}^{-1})$',
                                ylabel=r'$D^3\,\nu_*(\mu|\theta)\,\textswab{S}(\mu)$')
        bovy_plot.bovy_end_print(plotname)
    for ff in range(len(ls)):
        spl= interpolate.InterpolatedUnivariateSpline(dust._GREEN15DISTMODS,
                                                      area[ff]*dust._GREEN15DISTS**3.,
                                                  k=5)
        fthder= [spl.derivatives(dm)[4] for dm in dust._GREEN15DISTMODS]
        print "Simpson error (%f,%f)= %g, volume = %g" % (ls[ff],bs[ff],0.5**4./180.*numpy.mean(numpy.fabs(fthder))/integrate.simps(area[ff]*dust._GREEN15DISTS**3.,dx=0.5),numpy.sum(area[ff]*dust._GREEN15DISTS**3.))
    return None

if __name__ == '__main__':
    plot_distanceintegral_fewfields(sys.argv[1], # savefilename
                                    sys.argv[2], # plotfilename
                                    len(sys.argv) > 3) # APOGEE selection?
