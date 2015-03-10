###############################################################################
# define_rcsample: definitions of the sample used
###############################################################################
import numpy
import esutil
from galpy.util import bovy_coords
import apogee.tools.read as apread
from apogee.samples.rc import rcdist
import isodist
_FEHTAG= 'FE_H'
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
                                                    -data['FE_H']
    return data
    
