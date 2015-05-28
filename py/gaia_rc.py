###############################################################################
# gaia_rc.py: functions for the RC in Gaia
###############################################################################
import numpy
from scipy import optimize
import isodist
from apogee.samples import rc
from apogee.util import localfehdist
def G_jordi(g,gz):
    """G in terms of g and g-z From Jordi et al. (2010)"""
    return g-0.1154-0.4175*gz-0.0497*gz**2.+0.0016*gz**3. 
def G(age,Z,iso):
    """G for a given age and metallicity"""
    Zs= iso.Zs()
    tZ= Zs[numpy.argmin(numpy.fabs(Z-Zs))]
    p= iso(age,tZ)
    jk= p['J']-p['Ks']
    indx= (jk < 0.8)*(jk > 0.5)\
        *(tZ <= 0.06)\
        *(tZ <= rc.jkzcut(jk,upper=True))\
        *(tZ >= rc.jkzcut(jk))\
        *(p['logg'] >= rc.loggteffcut(10.**p['logTe'],tZ,upper=False))\
        *(p['logg'] <= rc.loggteffcut(10.**p['logTe'],tZ,upper=True))
    outG= G_jordi(p['g'],p['g']-p['z'])
    # Average over the IMF
    sindx= numpy.argsort(p['M_ini'])
    outG= outG[sindx]
    int_IMF= p['int_IMF'][sindx]
    w= (numpy.roll(int_IMF,-1)-int_IMF)/(int_IMF[-1]-int_IMF[0])
    w= w[indx]
    outG= outG[indx]
    return (numpy.nansum(w*outG)/numpy.nansum(w),
            numpy.sqrt(numpy.nansum(w*outG**2.)/numpy.nansum(w)
                       -(numpy.nansum(w*outG)/numpy.nansum(w))**2.))
def Gdist(tG,iso):
    """The distribution of G for the local metallicity distribution
    REWRITE TO TAKE A Z(G) FUNCTION OBTAINED THROUGH INTERPOLATION"""
    try:
        tZ= optimize.brenth(lambda x: G(9.7,x,iso)[0]-tG,0.0001,0.06)
    except ValueError:
        return 0.
    # Add Jacobian
    return localfehdist(isodist.Z2FEH(tZ,zsolar=0.017))*1.
def load_iso():
    Zs= [0.0035,0.017,0.035]
    return isodist.PadovaIsochrone(type='sdss-2mass',parsec=True,
                                   Z=Zs)
