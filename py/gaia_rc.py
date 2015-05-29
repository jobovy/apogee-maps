###############################################################################
# gaia_rc.py: functions for the RC in Gaia
###############################################################################
import numpy
from scipy import interpolate
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
def Gdist(tG,ZG):
    """The distribution of G for the local metallicity distribution
    REWRITE TO TAKE A Z(G) FUNCTION OBTAINED THROUGH INTERPOLATION"""
    try:
        tZ= ZG(tG)
        tjac= numpy.fabs(ZG.derivatives(tG)[0]/tZ/numpy.log(10.))
    except ValueError:
        return 0.
    # Add Jacobian
    return localfehdist(isodist.Z2FEH(tZ,zsolar=0.017))*tjac
def sample_Gdist(iso,n=1000):
    """Sample from the distribution of MG"""
    # First calculate the ditribution
    Gmin, Gmax= 0.1, 0.9
    Gs= numpy.linspace(Gmin,Gmax,201)
    ZG= load_ZG(iso)
    pG= numpy.array([Gdist(g,ZG) for g in Gs])
    pGmax= numpy.nanmax(pG)
    # Now rejection sample
    out= []
    while len(out) < n:
        # generate uniformly
        new= numpy.random.uniform(size=2*(n-len(out)))*(Gmax-Gmin)+Gmin
        pnew= numpy.array([Gdist(nn,ZG) for nn in new])
        out.extend(list(new[numpy.random.uniform(size=2*(n-len(out)))*pGmax \
                                < pnew]))
    return numpy.array(out[:n])
def load_iso():
    Zs= numpy.arange(0.0005,0.0605,0.0005)
    return isodist.PadovaIsochrone(type='sdss-2mass',parsec=True,
                                   Z=Zs)
def load_ZG(iso,s=0.00005):
    """Return a function that gives Z as a function of G"""
    Zs= iso.Zs()
    Zs= Zs[Zs < 0.05]
    tage= iso.logages()[62]
    Gs= numpy.array([G(tage,z,iso)[0] for z in Zs])
    sindx= numpy.argsort(Gs)
    Zs= Zs[sindx]
    Gs= Gs[sindx]
    goodIndx= True-numpy.isnan(Gs)
    return interpolate.UnivariateSpline(Gs[goodIndx],Zs[goodIndx],k=3,
                                        s=s)
