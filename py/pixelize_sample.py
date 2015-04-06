###############################################################################
# pixelize_sample: 
###############################################################################
import math
import numpy
from galpy.util import bovy_plot
class pixelXY:
    """Class that pixelizes the data sample in X and Y"""
    def __init__(self,data,dx=1.,dy=1.,xmin=2.5,xmax=15.5,
                 ymin=-10.5,ymax=10.5,rphi=False):
        """
        NAME:
           __init__
        PURPOSE:
           initialize the pixelizeXY
        INPUT:
           data - the data sample
           dx, dy= pixel size in kpc
           xmin, xmax, ymin, ymax= minimum and maximum x and y
        OUTPUT:
           object with pixelized data
        HISTORY:
           2013-10-06 - Written - Bovy (IAS)
        """
        self.data= data
        self.dx= dx
        self.dy= dy        
        self.xmin= xmin
        self.xmax= xmax
        self.ymin= ymin
        self.ymax= ymax
        #edges in X and Y
        self.xedges= numpy.arange(self.xmin,self.xmax+0.01,self.dx)
        self.yedges= numpy.arange(self.ymin,self.ymax+0.01,self.dy)
        #Compute X and Y
        self.rphi= rphi
        if rphi:
            self.x= data['RC_GALR']
            self.y= data['RC_GALPHI']/numpy.pi*180.
            self.y[self.y > 180.]-= 360.
        else:
            self.x= data['RC_GALR']*numpy.cos(data['RC_GALPHI'])
            self.y= data['RC_GALR']*numpy.sin(data['RC_GALPHI'])        
        return None

    def __call__(self,*args,**kwargs):
        """
        NAME:
           __call__
        PURPOSE:
           return the part of the sample in a X,Y pixel
        INPUT:
           X,Y
        OUTPUT:
           returns data recarray in the bin that feh and afe are in
        HISTORY:
           2013-10-06 - Written - Bovy (IAS)
        """
        #Find bin
        xbin= int(math.floor((args[0]-self.xmin)/self.dx))
        ybin= int(math.floor((args[1]-self.ymin)/self.dy))
        #Return data
        return self.data[(self.x > self.xedges[xbin])\
                             *(self.x <= self.xedges[xbin+1])\
                             *(self.y > self.yedges[ybin])\
                             *(self.y <= self.yedges[ybin+1])]

    def callIndx(self,*args,**kwargs):
        """
        NAME:
           callIndx
        PURPOSE:
           return index of the part of the sample in an X and Y pixel
        INPUT:
           X, Y
        OUTPUT:
           returns index into data recarray in the bin that feh and afe are in
        HISTORY:
           2013-10-060 - Written - Bovy (IAS)
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
           return the index corresponding to a X value
        INPUT:
           x
        OUTPUT:
           index
        HISTORY:
           2012-10-06 - Written - Bovy (IAS)
        """
        return int(math.floor((x-self.xmin)/self.dx))

    def yindx(self,y):
        """
        NAME:
           yindx
        PURPOSE:
           return the index corresponding to a y value
        INPUT:
           Y
        OUTPUT:
           index
        HISTORY:
           2013-10-06 - Written - Bovy (IAS)
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
           2013-10-06 - Written - Bovy (IAS)
        """
        #First create 2D
        gx= numpy.linspace(self.xmin+self.dx/2.,self.xmax-self.dx/2.,
                           int((self.xmax-self.xmin)/self.dx))
        gy= numpy.linspace(self.ymin+self.dy/2.,self.ymax-self.dy/2.,
                           int((self.ymax-self.ymin)/self.dy))
        z2d= numpy.empty((int((self.xmax-self.xmin)/self.dx),
                          int((self.ymax-self.ymin)/self.dy)))
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
                z2d[ii,:]-= numpy.median(z2d[ii,True-numpy.isnan(z2d[ii,:])])
        if justcalc:
            if returnz:
                return z2d
            else:
                return None
        #Now plot
        xrange=[self.xmin,self.xmax]
        yrange=[self.ymin,self.ymax]
        if not kwargs.has_key('colorbar'):
            kwargs['colorbar']= True
        if not kwargs.has_key('shrink'):
            kwargs['shrink']= 0.78
        if not kwargs.has_key('vmin'):
            kwargs['vmin']= numpy.nanmin(z2d)
        if not kwargs.has_key('vmax'):
            kwargs['vmax']= numpy.nanmax(z2d)
        if self.rphi:
            xlabel=r'$R\,(\mathrm{kpc})$'
            ylabel=r'$\mathrm{Galactocentric\ azimuth}\,(\mathrm{deg})$'
        else:
            xlabel=r'$X_{\mathrm{GC}}\,(\mathrm{kpc})$'
            ylabel=r'$Y_{\mathrm{GC}}\,(\mathrm{kpc})$'
        out= bovy_plot.bovy_dens2d(z2d.T,origin='lower',cmap='jet',
                                   interpolation='nearest',
                                   xlabel=xlabel,ylabel=ylabel,
                                   xrange=xrange,yrange=yrange,
                                   **kwargs)
        if returnz:
            return z2d
        else:
            return out
        
    
