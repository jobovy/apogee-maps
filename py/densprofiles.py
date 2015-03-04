###############################################################################
# densprofiles: a collection of parameterized density profiles
###############################################################################
import numpy
from galpy.util import bovy_coords



@glonDecorator
def expdisk(R,z,phi,glon=False,
            params=[1./3.,1./0.3],log=False):
    """
    NAME:
       expdisk
    PURPOSE:
       density of an exponential disk
    INPUT:
       R,z,phi - Galactocentric cylindrical coordinates or (l/rad,b/rad,D/kpc)
       glon= (False) if True, input coordinates above are (l,b,D)
       params= parameters [1/hR,1/hz]
       log= (False) if True, return the log of the density
    OUTPUT:
       density or log density
    HISTORY:
       2015-03-04 - Written - Bovy (IAS)
    """
    if log:
        return -params[0]*R-params[1]*numpy.fabs(z)
    else:
        return numpy.exp(-params[0]*R-params[1]*numpy.fabs(z))

