###############################################################################
# define_rcsample: definitions of the sample used
###############################################################################
import numpy
import esutil
from galpy.util import bovy_coords
import apogee.tools.read as apread
from apogee.samples.rc import rcdist
import isodist
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
    return data
    
