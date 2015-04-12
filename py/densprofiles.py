###############################################################################
# densprofiles: a collection of parameterized density profiles
###############################################################################
from functools import wraps
import numpy
import healpy
from galpy.util import bovy_coords
_R0= 8.
_Zsun= 0.025

# Input decorators
def glonDecorator(func):
    """Decorator to convert input in (l/rad,b/rad,D/kpc) to (R,z,phi)"""
    @wraps(func)
    def glon_wrapper(*args,**kwargs):
        if kwargs.pop('glon',False):
            XYZ= bovy_coords.lbd_to_XYZ(args[0],args[1],args[2],degree=False)
            R,phi,z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],XYZ[:,1],XYZ[:,2],
                                                  Xsun=_R0,Zsun=_Zsun)
        else:
            R,phi,z= args[0],args[1],args[2]
        return func(R,phi,z,*args[3:],**kwargs)
    return glon_wrapper

def scalarDecorator(func):
    """Decorator to deal with scalar input"""
    @wraps(func)
    def scalar_wrapper(*args,**kwargs):
        if numpy.array(args[0]).shape == ():
            scalarOut= True
            newargs= ()
            for ii in range(len(args)):
                newargs= newargs+(numpy.array([args[ii]]),)
            args= newargs
        else:
            scalarOut= False
        result= func(*args,**kwargs)
        if scalarOut:
            return result[0]
        else:
            return result
    return scalar_wrapper

############################### LOGIT FOR AMPLITUDES ##########################
def logit(p):
    """The logit functions"""
    return numpy.log(p/(1.-p))
def ilogit(x):
    """The reverse logit"""
    return 1./(1.+numpy.exp(-x))

################################# HEALPIX MAPS ################################
def healpixelate(dist,densprofile,params=None,nside=512,nest=True):
    """
    NAME:
       healpixelate
    PURPOSE:
       Pixelate a density profile at a given distance from the Sun on a HEALPIX grid
    INPUT:
       dist - distance in kpc
       densprofile - density profile function from this module
       params= parameter array of the density profile
       nside= (512) HEALPIX nside
       nest= (True) True: nested pixelation; False: ring pixelation
    OUTPUT:
       density on grid (1D array)
    HISTORY:
       2015-03-04 - Written - Bovy (IAS)
    """
    npix= healpy.pixelfunc.nside2npix(nside)
    theta,phi= healpy.pixelfunc.pix2ang(nside,numpy.arange(npix),nest=nest)
    return densprofile(phi,numpy.pi/2.-theta,dist*numpy.ones(npix),glon=True,
                       params=params)   

def powspec(dist,densprofile,params=None,nside=512):
    """
    NAME:
       powspec
    PURPOSE:
       calculate the angular power spectrum of a density profile at a given distance
    INPUT:
       dist - distance in kpc
       densprofile - density profile function from this module
       params= parameter array of the density profile
       nside= (512) HEALPIX nside
    OUTPUT:
       (l,C_l)
    HISTORY:
       2015-03-04 - Written - Bovy (IAS)
    """
    dmap= healpixelate(dist,densprofile,params=params,nside=nside,nest=False)
    cl= healpy.sphtfunc.anafast(dmap-numpy.mean(dmap),pol=False)
    return (numpy.arange(len(cl)),cl)

############################### DENSITY PROFILES ##############################
@scalarDecorator
@glonDecorator
def expdisk(R,phi,z,glon=False,
            params=[1./3.,1./0.3],log=False):
    """
    NAME:
       expdisk
    PURPOSE:
       density of an exponential disk
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz]
       log= (False) if True, return the log of the density
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-04 - Written - Bovy (IAS)
    """
    if log:
        return -params[0]*(R-_R0)-params[1]*numpy.fabs(z)
    else:
        return numpy.exp(-params[0]*(R-_R0)-params[1]*numpy.fabs(z))

@scalarDecorator
@glonDecorator
def expdiskplusconst(R,phi,z,glon=False,
                     params=[1./3.,1./0.3,0.1]):
    """
    NAME:
       expdiskplusconst
    PURPOSE:
       density of an exponential disk plus a constant
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,log(amp)]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    return numpy.exp(-params[0]*(R-_R0)-params[1]*numpy.fabs(z))/2.\
        *numpy.fabs(params[1])\
        +numpy.exp(params[2])/24.

@scalarDecorator
@glonDecorator
def twoexpdisk(R,phi,z,glon=False,
               params=[1./3.,1./0.3,1./4.,1./0.5,logit(0.1)]):
    """
    NAME:
       twoexpdisk
    PURPOSE:
       density of a sum of two exponential disks
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,1/hz2,logit(amp2)]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    amp= ilogit(params[4])
    return (1.-amp)/2.*numpy.fabs(params[1])\
        *numpy.exp(-params[0]*(R-_R0)-params[1]*numpy.fabs(z))\
        +amp/2.*params[3]*numpy.exp(-params[2]*(R-_R0)-params[3]*numpy.fabs(z))

@scalarDecorator
@glonDecorator
def brokenexpdisk(R,phi,z,glon=False,
                  params=[1./3.,1./0.3,1./4.,numpy.log(10.)]):
    """
    NAME:
       brokenexpdisk
    PURPOSE:
       density of a broken exponential disk (two scale lengths)
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    sz= z[R <= Rb]
    bz= z[R > Rb]
    out[R <= Rb]= \
        numpy.fabs(params[1])/2.*numpy.exp(-params[0]*(sR-_R0)-params[1]*numpy.fabs(sz))
    out[R > Rb]=\
        numpy.exp(-params[2]*(bR-_R0)-params[1]*numpy.fabs(bz))\
        *numpy.fabs(params[1])/2.*numpy.exp(params[2]*(Rb-_R0)-params[0]*(Rb-_R0))
    return out

@scalarDecorator
@glonDecorator
def tribrokenexpdisk(R,phi,z,glon=False,
                     params=[1./3.,1./0.3,1./4.,numpy.log(10.)]):
    """
    NAME:
       tribrokenexpdisk
    PURPOSE:
       density of a broken exponential disk (two scale lengths, one negative, one positive)
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    sz= z[R <= Rb]
    bz= z[R > Rb]
    out[R <= Rb]= \
        numpy.fabs(params[1])/2.*numpy.exp(numpy.fabs(params[0])*(sR-_R0)-params[1]*numpy.fabs(sz))
    out[R > Rb]=\
        numpy.exp(-numpy.fabs(params[2])*(bR-_R0)-params[1]*numpy.fabs(bz))\
        *numpy.fabs(params[1])/2.*numpy.exp(numpy.fabs(params[2])*(Rb-_R0)+numpy.fabs(params[0])*(Rb-_R0))
    return out

@scalarDecorator
@glonDecorator
def symbrokenexpdisk(R,phi,z,glon=False,
                     params=[1./3.,1./0.3,numpy.log(10.)]):
    """
    NAME:
       symbrokenexpdisk
    PURPOSE:
       density of a broken exponential disk, symmetric around a break
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-04-06 - Written - Bovy (IAS)
    """
    return brokenexpdisk(R,phi,z,glon=False,
                         params=[-params[0],params[1],2.*params[0],params[2]])

@scalarDecorator
@glonDecorator
def brokenexpflaredisk(R,phi,z,glon=False,
                       params=[1./3.,1./0.3,1./4.,numpy.log(10.),1./5.]):
    """
    NAME:
       brokenexpflaredisk
    PURPOSE:
       density of a broken exponential disk (two scale lengths), with a flaring vertical profile
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    tinvhz= params[1]*numpy.exp((R-_R0)*params[4])
    out[R <= Rb]= numpy.exp(-params[0]*(sR-_R0))
    out[R > Rb]= numpy.exp(-params[2]*(bR-_R0))\
        *numpy.exp(params[2]*(Rb-_R0)-params[0]*(Rb-_R0))
    return numpy.fabs(tinvhz)/2.*out*numpy.exp(-tinvhz*numpy.fabs(z))

@scalarDecorator
@glonDecorator
def tribrokenexpflaredisk(R,phi,z,glon=False,
                          params=[1./3.,1./0.3,1./4.,numpy.log(10.),-1./5.]):
    """
    NAME:
       tribrokenexpflaredisk
    PURPOSE:
       density of a broken exponential disk (two scale lengths, one negative, one positive)
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    tinvhz= params[1]*numpy.exp((R-_R0)*params[4])
    out[R <= Rb]= numpy.exp(numpy.fabs(params[0])*(sR-_R0))
    out[R > Rb]= numpy.exp(-numpy.fabs(params[2])*(bR-_R0))\
        *numpy.exp(numpy.fabs(params[2])*(Rb-_R0)\
                       +numpy.fabs(params[0])*(Rb-_R0))
    return numpy.fabs(tinvhz)/2.*out*numpy.exp(-tinvhz*numpy.fabs(z))

@scalarDecorator
@glonDecorator
def brokentwoexpdisk(R,phi,z,glon=False,
                     params=[1./3.,1./0.3,1./4.,numpy.log(10.),
                             1./0.8,logit(0.1)]):
    """
    NAME:
       brokentwoexpdisk
    PURPOSE:
       density of a broken exponential disk (two scale lengths), with a 
       vertical density profile consisting of two scale heights
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak],1/hz2,logit(amp2)]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    out[R <= Rb]= numpy.exp(-params[0]*(sR-_R0))
    out[R > Rb]=\
        numpy.exp(-params[2]*(bR-_R0))\
        *numpy.exp(params[2]*(Rb-_R0)-params[0]*(Rb-_R0))
    amp= ilogit(params[4])
    return out*((1.-amp)/2.*numpy.fabs(params[1])\
                    *numpy.exp(-params[1]*numpy.fabs(z))\
                    +amp/2.*numpy.fabs(params[5])\
                    *numpy.exp(-params[5]*numpy.fabs(z)))

@scalarDecorator
@glonDecorator
def tribrokentwoexpdisk(R,phi,z,glon=False,
                        params=[1./3.,1./0.3,1./4.,numpy.log(10.),logit(0.1),
                             1./0.8]):
    """
    NAME:
       tribrokentwoexpdisk
    PURPOSE:
       density of a broken exponential disk (two scale lengths, one negative, one positive)
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-24 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    out[R <= Rb]= \
        numpy.exp(numpy.fabs(params[0])*(sR-_R0))
    out[R > Rb]= numpy.exp(-numpy.fabs(params[2])*(bR-_R0))\
        *numpy.exp(numpy.fabs(params[2])*(Rb-_R0)+numpy.fabs(params[0])*(Rb-_R0))
    amp= ilogit(params[4])
    return out*((1.-amp)/2.*numpy.fabs(params[1])\
                    *numpy.exp(-params[1]*numpy.fabs(z))\
                    +amp/2.*numpy.fabs(params[5])\
                    *numpy.exp(-params[5]*numpy.fabs(z)))

@scalarDecorator
@glonDecorator
def gaussexpdisk(R,phi,z,glon=False,
                 params=[1./3.,1./0.3,numpy.log(10.)]):
    """
    NAME:
       gaussexpdisk
    PURPOSE:
       density as a Gaussian in radius and an exponential vertically
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,log[Rmax]]
    OUTPUT:
       density
    HISTORY:
       2015-03-28 - Written - Bovy (IAS)
    """
    Rm= numpy.exp(params[2])
    return numpy.fabs(params[1])/2.*numpy.exp(-params[1]*numpy.fabs(z))\
        *numpy.exp(-params[0]**2./2.*((R-Rm)**2.-(_R0-Rm)**2.))

@scalarDecorator
@glonDecorator
def brokenquadexpdisk(R,phi,z,glon=False,
                      params=[1./3.,1./0.3,1./4.,numpy.log(10.)]):
    """
    NAME:
       brokenquadexpdisk
    PURPOSE:
       density of a broken exponential squared disk (two scale lengths)
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-04-06 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    sz= z[R <= Rb]
    bz= z[R > Rb]
    out[R <= Rb]= \
        numpy.fabs(params[1])/2.*numpy.exp(-params[0]**2./2.*(sR-_R0)**2.-params[1]*numpy.fabs(sz))
    out[R > Rb]=\
        numpy.exp(-params[2]**2./2.*(bR-_R0)**2.-params[1]*numpy.fabs(bz))\
        *numpy.fabs(params[1])/2.*numpy.exp(params[2]**2./2.*(Rb-_R0)**2.-params[0]**2./2.*(Rb-_R0)**2.)
    return out

@scalarDecorator
@glonDecorator
def symbrokenquadexpdisk(R,phi,z,glon=False,
                         params=[1./3.,1./0.3,numpy.log(10.)]):
    """
    NAME:
       symbrokenquadexpdisk
    PURPOSE:
       density of a broken exponential squared disk, symmetric around a break
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,log[Rbreak]]
    OUTPUT:
       density or log density
    HISTORY:
       2015-04-06 - Written - Bovy (IAS)
    """
    return brokenquadexpdisk(R,phi,z,glon=False,
                             params=[params[0],params[1],
                                     -2.*params[0],params[2]])

def logspiral(R,phi,tanp=numpy.tan(9.4/180.*numpy.pi),
              Rref=9.9,phiref=14.2/180.*numpy.pi,
              width=0.38,phimin=-21./180.*numpy.pi,phimax=88./180.*numpy.pi):
    """
    NAME:
       logspiral
    PURPOSE:
       return the overdensity due to a logarithmic spiral (default parameters: Perseus arm)
    INPUT:
       R - Galactocentric radius (/kpc)
       phi - Galactocentric azimuth (/rad)
       tanp= tan of the pitch angle
       Rref= reference radius
       phiref= reference phi
       width= width in radius (/kpc)
       phimin, phimax= extent of the arm in azimuth
    OUTPUT:
       overdensity (unit amplitude)
    HISTORY:
       2015-04-06 - Written - Bovy (IAS)
    """
    phi[phi > 180.]-= 360.
    Rs= Rref*numpy.exp(-(phi-phiref)*tanp)
    w= width*numpy.sqrt(1.+tanp**2.)
    out= numpy.zeros_like(R)
    gindx= (phi > phimin)*(phi < phimax)
    out[gindx]= 1./numpy.sqrt(2.*numpy.pi)/w\
        *numpy.exp(-0.5*(R-Rs)[gindx]**2./w**2.)
    return out

@scalarDecorator
@glonDecorator
def brokenexpdiskfixedspiral(R,phi,z,glon=False,
                             params=[1./3.,1./0.3,1./4.,numpy.log(10.),
                                     numpy.log(0.1)]):
    """
    NAME:
       brokenexpdiskfixedspiral
    PURPOSE:
       density of a broken exponential disk (two scale lengths) + a fixed spiral
    INPUT:
       R,phi,z - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz,1/hR2,log[Rbreak],log[spamp]]
    OUTPUT:
       density
    HISTORY:
       2015-04-06 - Written - Bovy (IAS)
    """
    Rb= numpy.exp(params[3])
    spamp= numpy.exp(params[4])
    out= numpy.empty_like(R)
    sR= R[R <= Rb]
    bR= R[R > Rb]
    out[R <= Rb]= \
        numpy.exp(-params[0]*(sR-_R0))
    out[R > Rb]=\
        numpy.exp(-params[2]*(bR-_R0))\
        *numpy.exp(params[2]*(Rb-_R0)-params[0]*(Rb-_R0))
    return numpy.fabs(params[1])/2.*(out+spamp*logspiral(R,phi))\
        *numpy.exp(-params[1]*numpy.fabs(z))

