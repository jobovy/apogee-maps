import math
import numpy
import statsmodels.api as sm
lowess= sm.nonparametric.lowess
import esutil
from galpy.util import bovy_coords, bovy_plot
from scipy.interpolate import interp1d,UnivariateSpline
import apogee.tools.read as apread
import isodist
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle 
import apogee.tools.read as apread
from apogee.select import apogeeSelect
from astropy.io import fits
from astropy.table import Table, join
_R0= 8. # kpc
_Z0= 0.025 # kpc
_FEHTAG= 'FE_H'
_AFETAG= 'AVG_ALPHAFE'
_AFELABEL= r'$[\left([\mathrm{O+Mg+Si+S+Ca}]/5\right)/\mathrm{Fe}]$'

catpath = '../catalogues/'

selectFile= '../savs/selfunc-nospdata.sav'

if os.path.exists(selectFile):
    with open(selectFile,'rb') as savefile:
        apo= pickle.load(savefile)


def get_rgbsample(loggcut = [1.8, 3.0],
				  teffcut = [0, 10000], 
				  add_ages = False, 
				  agetype='Martig', 
				  apply_corrections=False, 
				  distance_correction=False,
				  verbose = False):
	"""
	Get a clean sample of dr12 APOGEE data with Michael Haydens distances
	---
	INPUT:
		None
	OUTPUT:
		Clean rgb sample with added distances
	HISTORY:
		Started - Mackereth 02/06/16 
	"""
	#get the allStar catalogue using apogee python (exlude all bad flags etc)
	allStar = apread.allStar(rmcommissioning=True, 
							 exclude_star_bad=True, 
							 exclude_star_warn=True, 
							 main=True,
							 ak=True, 
							 adddist=False)
	#cut to a 'sensible' logg range (giants which are not too high on the RGB)
	allStar = allStar[(allStar['LOGG'] > loggcut[0])&(allStar['LOGG'] < loggcut[1])&
					  (allStar['TEFF'] > teffcut[0])&(allStar['TEFF'] < teffcut[1])]
	if verbose == True:
		print str(len(allStar))+' Stars before Distance catalogue join (after Log(g) cut)'
	#load the distance VAC
	dists = fits.open(catpath+'DR12_DIST_R-GC.fits')[1].data
	#convert to astropy Table
	allStar_tab = Table(data=allStar)
	dists_tab = Table(data=dists)
	#join table
	tab = join(allStar_tab, dists_tab, keys='APOGEE_ID', uniq_col_name='{col_name}{table_name}', table_names=['','2'])
	data = tab.as_array()
	data= esutil.numpy_util.add_fields(data,[('M_J', float),
                                             ('M_H', float),
                                             ('M_K', float),
                                             ('MH50_DIST', float),
                                             ('MH50_GALR', float),
                                             ('MH50_GALZ', float),
                                             ('MH50_GALPHI', float),
                                             ('AVG_ALPHAFE', float)])
	data['MH50_DIST'] = (10**((data['HAYDEN_DISTMOD_50']+5)/5))/1e3
	
	if distance_correction == True:	
		data['MH50_DIST'] *= 1.05
	XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['MH50_DIST'],
                                degree=True)
	R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.,Zsun=0.025)
	data['MH50_GALR']= R
	data['MH50_GALPHI']= phi
	data['MH50_GALZ']= Z
	data['M_J'] = data['J0']-data['HAYDEN_DISTMOD_50']
	data['M_H'] = data['H0']-data['HAYDEN_DISTMOD_50']
	data['M_K'] = data['K0']-data['HAYDEN_DISTMOD_50']
	data['AVG_ALPHAFE'] = avg_alphafe_dr12(data)
	data[_FEHTAG] += -0.1
	#remove locations not in the apogee selection function (FIND OUT WHATS UP HERE)
	data = data[np.in1d(data['LOCATION_ID'], apo.list_fields())]
	# Remove locations outside of the Pan-STARRS dust map
	# In the Southern hemisphere
	data= data[data['LOCATION_ID'] != 4266] #240,-18
	data= data[data['LOCATION_ID'] != 4331] #5.5,-14.2
	data= data[data['LOCATION_ID'] != 4381] #5.2,-12.2
	data= data[data['LOCATION_ID'] != 4332] #1,-4
	data= data[data['LOCATION_ID'] != 4329] #0,-5
	data= data[data['LOCATION_ID'] != 4351] #0,-2
	data= data[data['LOCATION_ID'] != 4353] #358,0
	data= data[data['LOCATION_ID'] != 4385] #358.6,1.4
	# Close to the ecliptic pole where there's no data (is it the ecliptic pole?
	data= data[data['LOCATION_ID'] != 4528] #120,30
	data= data[data['LOCATION_ID'] != 4217] #123,22.4
	#remove any non-finite magnitudes
	data = data[np.isfinite(data['M_H'])]
	if verbose == True:
		print str(len(data))+' Stars with distance measures (and in good fields...)'
	if add_ages == True:
		if agetype == 'Martig':
			ages = fits.open(catpath+'DR12_martigages.fits')[1].data
			idtag = '2MASS_ID'
		if agetype == 'Cannon':
			ages = fits.open(catpath+'RGB_Cannon_Ages.fits')[1].data
			ages = esutil.numpy_util.add_fields(ages,[('Age', float)])
			ages['Age'] = np.exp(ages['ln_age'])
			idtag = 'ID'
		ages_tab = Table(data=ages)
		ages_tab.rename_column(idtag, 'APOGEE_ID')
		tab = join( ages_tab,data, keys='APOGEE_ID', uniq_col_name='{col_name}{table_name}', table_names=['','2'])
		allStar_full = tab.as_array()
		data = allStar_full
		if verbose == True:
			print str(len(data))+' Stars with ages'
	if apply_corrections == True:
		martig1 = np.genfromtxt(catpath+'martig2016_table1.txt', dtype=None, names=True, skip_header=2)
		fit = lowess(np.log10(martig1['Age_out']),np.log10(martig1['Age_in']))
		xs = np.linspace(-0.3,1.2,100)
		xsinterpolate = interp1d(xs,xs)
		fys = fit[:,0]-xsinterpolate(fit[:,1])
		interp = UnivariateSpline(fit[:,1], fys)
		corr_age = np.log10(data['Age'])+(interp(np.log10(data['Age'])))
		corr_age = 10**corr_age
		data['Age'] = corr_age
	return data
	
	
def avg_alphafe_dr12(data):    
    weight_o= np.ones(len(data))
    weight_s= np.ones(len(data))
    weight_si= np.ones(len(data))
    weight_ca= np.ones(len(data))
    weight_mg= np.ones(len(data))
    weight_o[data['O_H'] == -9999.0]= 0.
    weight_s[data['S_H'] == -9999.0]= 0.
    weight_si[data['SI_H'] == -9999.0]= 0.
    weight_ca[data['CA_H'] == -9999.0]= 0.
    weight_mg[data['MG_H'] == -9999.0]= 0.
    return (weight_o*data['O_H']+weight_s*data['S_H']
            +weight_si*data['SI_H']+weight_ca*data['CA_H']
            +weight_mg*data['MG_H'])/(weight_o+weight_s
                                      +weight_si+weight_ca
                                      +weight_mg)\
                                      -data['FE_H']-0.05

                                      
    
    
# Define the low-alpha, low-iron sample
def _lowlow_lowfeh(afe):
    # The low metallicity edge
    return -0.6
def _lowlow_highfeh(afe):
    # The high metallicity edge
    return -0.25
def _lowlow_lowafe(feh):
    # The low alpha edge (-0.15,-0.075) to (-0.5,0)
    return (0--0.075)/(-0.5--0.15)*(feh+0.1--0.15)-0.075
def _lowlow_highafe(feh): 
    # The high alpha edge (-0.15,0.075) to (-0.5,0.15)
    return (0.15-0.075)/(-0.5--0.15)*(feh+0.1--0.15)+0.075

def get_lowlowsample():
    """
    NAME:
       get_lowlowsample
    PURPOSE:
       get the RGB sample at low alpha, low iron
    INPUT:
       None so far
    OUTPUT:
       sample
    HISTORY:
       2015-03-18 - Started - Bovy (IAS)
       2016-07-02 - modification - Mackereth (LJMU)
    """
    # Get the full sample first
    data= get_rgbsample()
    # Now cut it
    lowfeh= _lowlow_lowfeh(0.)
    highfeh= _lowlow_highfeh(0.)
    indx= (data[_FEHTAG] > lowfeh)*(data[_FEHTAG] <= highfeh)\
        *(data[_AFETAG] > _lowlow_lowafe(data[_FEHTAG]))\
        *(data[_AFETAG] <= _lowlow_highafe(data[_FEHTAG]))
    return data[indx]

# Define the high-alpha sample
def _highalpha_lowfeh(afe):
    # The low metallicity edge
    return -0.8
def _highalpha_highfeh(afe):
    # The high metallicity edge
    return -0.2
def _highalpha_lowafe(feh):
    # The low alpha edge (-0.125,0.115) to (-0.6,0.215)
    return (0.2-0.1)/(-0.6--0.125)*(feh+0.1--0.125)+0.115
def _highalpha_highafe(feh):
    # The high alpha edge (-0.125,0.19) to (-0.6,0.29)
    return (0.275-0.175)/(-0.6--0.125)*(feh+0.1--0.125)+0.19

def get_highalphasample():
    """
    NAME:
       get_highalphasample
    PURPOSE:
       get the RC sample at high alpha
    INPUT:
       None so far
    OUTPUT:
       sample
    HISTORY:
       2015-03-24 - Started - Bovy (IAS)
    """
    # Get the full sample first
    data= get_rcsample()
    # Now cut it
    lowfeh= _highalpha_lowfeh(0.)
    highfeh= _highalpha_highfeh(0.)
    indx= (data[_FEHTAG] > lowfeh)*(data[_FEHTAG] <= highfeh)\
        *(data[_AFETAG] > _highalpha_lowafe(data[_FEHTAG]))\
        *(data[_AFETAG] <= _highalpha_highafe(data[_FEHTAG]))
    return data[indx]

# Define the solar sample
def _solar_lowfeh(afe):
    # The low metallicity edge
    return -0.2
def _solar_highfeh(afe):
    # The high metallicity edge
    return 0.
def _solar_lowafe(feh):
    # The low alpha edge (0.1,-0.075) to (-0.1,-0.075)
    return -0.075
def _solar_highafe(feh):
    # The high alpha edge (-0.15,0.1) to (0.1,0.05)
    return (0.1-0.05)/(-0.15-0.1)*(feh+0.1-0.1)+0.05

def get_solarsample():
    """
    NAME:
       get_solarsample
    PURPOSE:
       get the RC sample at solar abundances
    INPUT:
       None so far
    OUTPUT:
       sample
    HISTORY:
       2015-03-18 - Started - Bovy (IAS)
       2016-07-02 - modification - Mackereth (LJMU)
    """
    # Get the full sample first
    data= get_rgbsample()
    # Now cut it
    lowfeh= _solar_lowfeh(0.)
    highfeh= _solar_highfeh(0.)
    indx= (data[_FEHTAG] > lowfeh)*(data[_FEHTAG] <= highfeh)\
        *(data[_AFETAG] > _solar_lowafe(data[_FEHTAG]))\
        *(data[_AFETAG] <= _solar_highafe(data[_FEHTAG]))
    return data[indx]

# Define the high metallicity sample
def _highfeh_lowfeh(afe):
    # The low metallicity edge
    return 0.05
def _highfeh_highfeh(afe):
    # The high metallicity edge
    return 0.3
def _highfeh_lowafe(feh):
    # The low alpha edge (0.1,-0.075) to (-0.1,-0.075)
    return -0.075
def _highfeh_highafe(feh):
    # The high alpha edge (-0.15,0.1) to (0.1,0.05)
    return 0.05

def get_highfehsample():
    """
    NAME:
       get_highfehsample
    PURPOSE:
       get the RC sample at high [Fe/H]
    INPUT:
       None so far
    OUTPUT:
       sample
    HISTORY:
       2015-03-18 - Started - Bovy (IAS)
       2016-07-02 - modification - Mackereth (LJMU)
    """
    # Get the full sample first
    data= get_rgbsample()
    # Now cut it
    lowfeh= _highfeh_lowfeh(0.)
    highfeh= _highfeh_highfeh(0.)
    indx= (data[_FEHTAG] > lowfeh)*(data[_FEHTAG] <= highfeh)\
        *(data[_AFETAG] > _highfeh_lowafe(data[_FEHTAG]))\
        *(data[_AFETAG] <= _highfeh_highafe(data[_FEHTAG]))
    return data[indx]

   
def alphaedge(fehs):
    edge = np.zeros(len(fehs))
    edge[fehs < 0] = (0.12/-0.6)*fehs[fehs < 0]+0.03
    edge[fehs >= 0] = 0.03
    return edge
    
def highalphaedge(fehs):
    edge = np.zeros(len(fehs))
    edge[fehs < 0] = (-0.13/0.6)*fehs[fehs < 0]+0.04
    edge[fehs >= 0] = 0.04
    return edge

def lowalphaedge(fehs):
    edge = np.zeros(len(fehs))
    edge[fehs < 0] = (-0.10/0.6)*fehs[fehs < 0]+0.01
    edge[fehs >= 0] = 0.01
    return edge

def get_fehage(agebin = [0.,1.], fehbin = [0.,0.2], afebin = 'low', dr=None, agetype='Martig', apply_corrections=False, distance_correction=False):
    data = get_rgbsample(add_ages=True, agetype=agetype, apply_corrections=apply_corrections, distance_correction=distance_correction)
    if afebin == 'low':
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])\
            *(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])*(data[_AFETAG] < alphaedge(data[_FEHTAG]))
    if afebin == 'high':
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])\
            *(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])*(data[_AFETAG] >= alphaedge(data[_FEHTAG]))
    if afebin == 'highclean':
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])\
            *(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])*(data[_AFETAG] >= highalphaedge(data[_FEHTAG]))
    if afebin == 'lowclean':
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])\
            *(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])*(data[_AFETAG] <= lowalphaedge(data[_FEHTAG]))
    if afebin == 'lownew':
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])\
            *(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])*(data[_AFETAG] <= alphaedge(data[_FEHTAG])-0.025)
    if afebin == 'highnew':
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])\
            *(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])*(data[_AFETAG] >= alphaedge(data[_FEHTAG])+0.025)		
    if afebin == None:
        indx = (data['Age'] >= agebin[0])*(data['Age'] < agebin[1])*(data[_FEHTAG] >= fehbin[0])*(data[_FEHTAG] < fehbin[1])
    return data[indx]

def highalphalocus():
    data= get_rgbsample()
    indx= (data[_AFETAG] > (0.2-0.1)/(-0.6--0.125)*(data[_FEHTAG]+0.1--0.125)+0.11)\
        *(data[_FEHTAG] < -0.225)\
        +(data[_AFETAG] > 0.05/(-0.6--0.125)*(data[_FEHTAG]+0.1--0.125)+0.11)\
        *(data[_FEHTAG] >= -0.225)*(data[_FEHTAG] < 0.125)\
        +(data[_FEHTAG] >= 0.125)
    return lowess(data[_AFETAG][indx],data[_FEHTAG][indx],frac=0.6)
def lowalphalocus():
    data= get_rgbsample()
    indx= (data[_AFETAG] > (0.2-0.1)/(-0.6--0.125)*(data[_FEHTAG]+0.1--0.125)+0.11)\
        *(data[_FEHTAG] < -0.025)\
        +(data[_AFETAG] > 0.05/(-0.6--0.125)*(data[_FEHTAG]+0.1--0.125)+0.11)\
        *(data[_FEHTAG] >= -0.225)*(data[_FEHTAG] < 0.125)
    return lowess(data[_AFETAG][True-indx],data[_FEHTAG][True-indx],frac=0.6)



class MAPs:
    """Class that pixelizes the data sample in [Fe/H] and [a/Fe]"""
    def __init__(self,data=None,dfeh=0.1,dafe=0.05,fehmin=-0.75,fehmax=0.35,
                afemin=-0.075,afemax=0.275):
        """
        NAME:
           __init__
        PURPOSE:
           initialize the MAPs
        INPUT:
           data= (None) the data sample; if None, whole stat. RC sample
           dfeh, dafe= pixel size
           fehmin, fehmax, afemin, afemax= minimum and maximum FeH and AFe
        OUTPUT:
           object with pixelized data
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        if data is None: data= get_rcsample()
        self.data= data
        self.dx= dfeh
        self.dy= dafe
        self.xmin= fehmin
        self.xmax= fehmax
        self.ymin= afemin
        self.ymax= afemax
        # edges in X and Y
        self.xedges= numpy.arange(self.xmin,self.xmax+0.01,self.dx)
        self.yedges= numpy.arange(self.ymin,self.ymax+0.01,self.dy)
        # X and Y
        self.x= data[_FEHTAG]
        self.y= data[_AFETAG]
        return None

    def __call__(self,*args,**kwargs):
        """
        NAME:
           __call__
        PURPOSE:
           return the part of the sample in a (feh,afe) pixel
        INPUT:
           [Fe/H]
           [a/Fe]
        OUTPUT:
           returns data recarray in the bin that feh and afe are in
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        #Find bin
        xbin= int(math.floor((args[0]-self.xmin)/self.dx))
        ybin= int(math.floor((args[1]-self.ymin)/self.dy))
        #Return data
        return self.data[(self.x > self.xedges[xbin])\
                             *(self.x <= self.xedges[xbin+1])\
                             *(self.y > self.yedges[ybin])\
                             *(self.y <= self.yedges[ybin+1])]

    def map(self):
        """
        NAME:
           map
        PURPOSE:
           yield a map
        INPUT:
           (none)
        OUTPUT:
           iterates over the MAPs
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        nx= int((self.xmax-self.xmin)/self.dx)
        ny= int((self.ymax-self.ymin)/self.dy)
        gx= numpy.linspace(self.xmin+self.dx/2.,self.xmax-self.dx/2.,nx)
        gy= numpy.linspace(self.ymin+self.dy/2.,self.ymax-self.dy/2.,ny)
        for ii in range(nx):
            for jj in range(ny):
                yield self(gx[ii],gy[jj])

    def callIndx(self,*args,**kwargs):
        """
        NAME:
           callIndx
        PURPOSE:
           return index of the part of the sample in an [Fe/H] and [a/Fe] pixel
        INPUT:
           [Fe/H]
           [a/Fe]
        OUTPUT:
           returns index into data recarray in the bin that [Fe/H] and [a/Fe] are in
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        #Find bin
        xbin= int(math.floor((args[0]-self.xmin)/self.dx))
        ybin= int(math.floor((args[1]-self.ymin)/self.dy))
        #Return data
        return (self.x > self.xedges[xbin])\
            *(self.x <= self.xedges[xbin+1])\
            *(self.y > self.yedges[ybin])\
            *(self.y <= self.yedges[ybin+1])

    def xindx(self,x):
        """
        NAME:
           xindx
        PURPOSE:
           return the index corresponding to a [Fe/H] value
        INPUT:
           [Fe/H]
        OUTPUT:
           index
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        return int(math.floor((x-self.xmin)/self.dx))

    def yindx(self,y):
        """
        NAME:
           yindx
        PURPOSE:
           return the index corresponding to a [a/Fe] value
        INPUT:
           [a/Fe]
        OUTPUT:
           index
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        return int(math.floor((y-self.ymin)/self.dy))

    def plot(self,quant,func=numpy.median,minnstar=20.,submediany=False,
             returnz=False,justcalc=False,
             **kwargs):
        """
        NAME:
           plot
        PURPOSE:
           make a plot of a quantity as a function of X and Y
        INPUT:
           quant - the quantity (string that returns the quantity, like 
           'METALS') or a function of the data
           func - function of quantity to plot
           minnstar= minimum number of stars (20)
           submeany= subtract the median y
           justcalc= (False) if True, do not plot
           bovy_plot.bovy_dens2d kwargs
        OUTPUT:
           plot to output device
        HISTORY:
           2015-04-06 - Written - Bovy (IAS)
        """
        #First create 2D
        nx= int((self.xmax-self.xmin)/self.dx)
        ny= int((self.ymax-self.ymin)/self.dy)
        gx= numpy.linspace(self.xmin+self.dx/2.,self.xmax-self.dx/2.,nx)
        gy= numpy.linspace(self.ymin+self.dy/2.,self.ymax-self.dy/2.,ny)
        z2d= numpy.empty((nx,ny))
        if isinstance(quant,numpy.ndarray):
            z2d= numpy.reshape(quant,(nx,ny))
            for ii in range(z2d.shape[0]):
                for jj in range(z2d.shape[1]):
                    tdata= self(gx[ii],gy[jj])
                    if len(tdata) < minnstar:
                        z2d[ii,jj]= numpy.nan
        else:
            nbins= 0
            for ii in range(z2d.shape[0]):
                for jj in range(z2d.shape[1]):
                    tdata= self(gx[ii],gy[jj])
                    if len(tdata) < minnstar:
                        z2d[ii,jj]= numpy.nan
                    else:
                        nbins+= 1
                        if hasattr(quant, '__call__'):
                            z2d[ii,jj]= func(quant(tdata))
                        else:
                            z2d[ii,jj]= func(tdata[quant])
                if submediany:
                    z2d[ii,:]-= \
                        numpy.median(z2d[ii,True-numpy.isnan(z2d[ii,:])])
        if justcalc:
            if returnz:
                return z2d
            else:
                return None
        #Now plot
        xrange= kwargs.pop('xrange',[self.xmin,self.xmax])
        yrange= kwargs.pop('yrange',[self.ymin,self.ymax])
        if not kwargs.has_key('colorbar'):
            kwargs['colorbar']= True
        if not kwargs.has_key('shrink'):
            kwargs['shrink']= 0.78
        if not kwargs.has_key('vmin'):
            kwargs['vmin']= numpy.nanmin(z2d)
        if not kwargs.has_key('vmax'):
            kwargs['vmax']= numpy.nanmax(z2d)
        xlabel= r'$[\mathrm{Fe/H}]$'
        ylabel= _AFELABEL
        cmap= kwargs.pop('cmap','coolwarm')
        out= bovy_plot.bovy_dens2d(z2d.T,origin='lower',cmap=cmap,
                                   interpolation='nearest',
                                   xlabel=xlabel,ylabel=ylabel,
                                   xrange=xrange,yrange=yrange,
                                   **kwargs)
        if returnz:
            return z2d
        else:
            return out
        
    
