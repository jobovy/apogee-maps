#############################################################################
# calc_mass.py: calculate the total stellar mass in an observed sample
#############################################################################
import numpy
import fitDens
import densprofiles
def calcDiskMass(data,samples,
                 locations,effsel,distmods,
                 type='tribrokenexpflare'):
    """
    NAME:
       calcDiskMass
    PURPOSE:
       calculate the local surface density for a set of density profiles
    INPUT:
       data - the data array
       samples - an array [nparam,ndens] of density-profile parameters
       locations - locations of the APOGEE effective selection function
       effsel - array (nloc,nD) of the effective selection function, includes area of the field
       distmods - grid of distance moduli on which the effective selection function is pre-computed
       type= ('exp') type of density profile to use
    OUTPUT:
       local surface density in Msol/pc^2
    HISTORY:
       2015-04-29 - Written - Bovy (IAS)
    """
    # Setup the density function and its initial parameters
    densfunc= fitDens._setup_densfunc(type)
    # Setup the integration of the effective volume
    effsel, Rgrid, phigrid, zgrid= \
        fitDens._setup_effvol(locations,effsel,distmods)
    out= []
    for sample in samples.T:
        # Setup the density function, fix the normalization for Rb < R0
        if 'tribroken' in type and numpy.exp(sample[3]) < densprofiles._R0:
            norm= numpy.exp(-(sample[0]+sample[2])\
                                 *(numpy.exp(sample[3])-densprofiles._R0))
        else:
            norm= 1.
        tdensfunc= lambda x,y,z: densfunc(x,y,z,params=sample)*norm
        out.append(calcDiskMass_single(data,tdensfunc,
                                       effsel,Rgrid,phigrid,zgrid))
    return numpy.array(out)*12500.

def calcDiskMass_single(data,
                        densfunc,effsel,Rgrid,phigrid,zgrid):
    """
    NAME:
       calcDiskMass_single
    PURPOSE:
       calculate the local surface density for a single density profile
    INPUT:
       data - the data array
       densfunc - function that returns the density when called with R,phi,z
       effsel - array (nloc,nD), includes D^2*Delta D factor
       Rgrid, phigrid, zgrid - array (nloc,nD) of the cylindrical Galactocentric coordinates corresponding to the (l_loc,b_loc,D) of the effective selection function      
    OUTPUT:
       local surface density in Msol/pc^2
    HISTORY:
       2015-04-29 - Written - Bovy (IAS)
    """
    # Calculate the effective volume
    pred= fitDens.effvol(densfunc,effsel,Rgrid,phigrid,zgrid)
    # Area included in effvol is in deg^2
    return len(data)/pred/10.**6.*(180./numpy.pi)**2.
