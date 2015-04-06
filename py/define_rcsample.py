###############################################################################
# define_rcsample: definitions of the sample used
###############################################################################
import math
import numpy
import esutil
from galpy.util import bovy_coords, bovy_plot
import apogee.tools.read as apread
from apogee.samples.rc import rcdist
import isodist
_R0= 8. # kpc
_Z0= 0.025 # kpc
_FEHTAG= 'FE_H'
_AFETAG= 'AVG_ALPHAFE'
_AFELABEL= r'$[\left([\mathrm{O+Mg+Si+S+Ca}]/5\right)/\mathrm{Fe}]$'
def get_rcsample():
    """
    NAME:
       get_rcsample
    PURPOSE:
       get the RC sample
    INPUT:
       None so far
    OUTPUT:
       sample
    HISTORY:
       2015-02-10 - Started - Bovy (IAS@KITP)
    """
    data= apread.rcsample()
    # Cut to statistical sample
    data= data[data['STAT'] == 1]
    # Add the M_H-based distances
    data= esutil.numpy_util.add_fields(data,[('RC_DIST_H', float),
                                             ('RC_DM_H', float),
                                             ('RC_GALR_H', float),
                                             ('RC_GALPHI_H', float),
                                             ('RC_GALZ_H', float)])
    rcd= rcdist()
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    data['RC_DIST_H']= rcd(jk,z,appmag=data['H0'],mh=True)
    data['RC_DM_H']= 5.*numpy.log10(data['RC_DIST_H'])+10.
    XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['RC_DIST_H'],
                                degree=True)
    R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.,Zsun=0.025)
    data['RC_GALR_H']= R
    data['RC_GALPHI_H']= phi
    data['RC_GALZ_H']= Z
    # Add the average alpha/Fe
    data= esutil.numpy_util.add_fields(data,[('AVG_ALPHAFE', float)])
    weight_o= numpy.ones(len(data))
    weight_s= numpy.ones(len(data))
    weight_si= numpy.ones(len(data))
    weight_ca= numpy.ones(len(data))
    weight_mg= numpy.ones(len(data))
    weight_o[data['O_H'] == -9999.0]= 0.
    weight_s[data['S_H'] == -9999.0]= 0.
    weight_si[data['SI_H'] == -9999.0]= 0.
    weight_ca[data['CA_H'] == -9999.0]= 0.
    weight_mg[data['MG_H'] == -9999.0]= 0.
    data['AVG_ALPHAFE']= (weight_o*data['O_H']+weight_s*data['S_H']
                          +weight_si*data['SI_H']+weight_ca*data['CA_H']
                          +weight_mg*data['MG_H'])/(weight_o+weight_s
                                                    +weight_si+weight_ca
                                                    +weight_mg)\
                                                    -data['FE_H']-0.05
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
    # Remove stars w/ DM < 8.49, because for standard candle RC, these cannot be in the sample
    data= data[data['RC_DM_H'] > 8.49]
    return data
    
# Define the low-alpha, low-iron sample
def _lowlow_lowfeh(afe):
    # The low metallicity edge
    return -0.5
def _lowlow_highfeh(afe):
    # The high metallicity edge
    return -0.15
def _lowlow_lowafe(feh):
    # The low alpha edge (-0.15,-0.075) to (-0.5,0)
    return (0--0.075)/(-0.5--0.15)*(feh--0.15)-0.075
def _lowlow_highafe(feh):
    # The high alpha edge (-0.15,0.075) to (-0.5,0.15)
    return (0.15-0.075)/(-0.5--0.15)*(feh--0.15)+0.075

def get_lowlowsample():
    """
    NAME:
       get_lowlowsample
    PURPOSE:
       get the RC sample at low alpha, low iron
    INPUT:
       None so far
    OUTPUT:
       sample
    HISTORY:
       2015-03-18 - Started - Bovy (IAS)
    """
    # Get the full sample first
    data= get_rcsample()
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
    return -0.7
def _highalpha_highfeh(afe):
    # The high metallicity edge
    return -0.1
def _highalpha_lowafe(feh):
    # The low alpha edge (-0.125,0.115) to (-0.6,0.215)
    return (0.2-0.1)/(-0.6--0.125)*(feh--0.125)+0.115
def _highalpha_highafe(feh):
    # The high alpha edge (-0.125,0.19) to (-0.6,0.29)
    return (0.275-0.175)/(-0.6--0.125)*(feh--0.125)+0.19

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
    return -0.1
def _solar_highfeh(afe):
    # The high metallicity edge
    return 0.1
def _solar_lowafe(feh):
    # The low alpha edge (0.1,-0.075) to (-0.1,-0.075)
    return -0.075
def _solar_highafe(feh):
    # The high alpha edge (-0.15,0.1) to (0.1,0.05)
    return (0.1-0.05)/(-0.15-0.1)*(feh-0.1)+0.05

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
    """
    # Get the full sample first
    data= get_rcsample()
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
    return 0.15
def _highfeh_highfeh(afe):
    # The high metallicity edge
    return 0.4
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
    """
    # Get the full sample first
    data= get_rcsample()
    # Now cut it
    lowfeh= _highfeh_lowfeh(0.)
    highfeh= _highfeh_highfeh(0.)
    indx= (data[_FEHTAG] > lowfeh)*(data[_FEHTAG] <= highfeh)\
        *(data[_AFETAG] > _highfeh_lowafe(data[_FEHTAG]))\
        *(data[_AFETAG] <= _highfeh_highafe(data[_FEHTAG]))
    return data[indx]

###############################################################################
# pixelization in (feh,afe)
###############################################################################
class MAPs:
    """Class that pixelizes the data sample in [Fe/H] and [a/Fe]"""
    def __init__(self,data=None,dfeh=0.1,dafe=0.05,fehmin=-0.65,fehmax=0.45,
                afemin=-0.075,afemax=0.325):
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
        print len(data)
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
        
    
