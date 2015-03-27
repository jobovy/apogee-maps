###############################################################################
# define_rcsample: definitions of the sample used
###############################################################################
import numpy
import esutil
from galpy.util import bovy_coords
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

