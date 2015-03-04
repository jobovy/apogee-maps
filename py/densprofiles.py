###############################################################################
# densprofiles: a collection of parameterized density profiles
###############################################################################
from functools import wraps
import numpy
from galpy.util import bovy_coords
_R0= 8.
_Zsun= 0.025

def glonDecorator(func):
    """Decorator to convert input in (l/rad,b/rad,D/kpc) to (R,z,phi)"""
    @wraps(func)
    def glon_wrapper(*args,**kwargs):
        if kwargs.get('glon',False):
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

