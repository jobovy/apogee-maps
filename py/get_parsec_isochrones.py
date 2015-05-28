from mechanize import Browser
import sys
import os
import time
import numpy
_ERASESTR= "                                                                                "
def get_parsec_isochrone_onez(z,
                              savefilename,
                              ages=[6.6,10.15,0.05],
                              photsys="tab_mag_odfnew/tab_mag_sloan_2mass.dat",
                              eta_reimers=0.2):
    """
    NAME:
       get_parsec_isochrone_onez
    PURPOSE:
       download PARSEC isochrone tables for a single metallicity and a bunch of ages
    INPUT:
       z - metallicity Z
       savefilename - filename to save the isochrone tables to
       ages= [6.6,1.15,0.05] log10 of minage, maxage, dage
       photsys= photometric system
       eta_reimers= (0.2) mass-loss efficiency parameter
    OUTPUT:
       saves the isochrone table as a gzipped file in savefilename
    HISTORY:
       2015-05-28 - Written based on old version - Bovy (IAS)
    """
    br= Browser()
    br.open('http://stev.oapd.inaf.it/cgi-bin/cmd_2.5')
    form= br.forms().next() #There is only one, hopefully!
    br.form= form
    br["photsys_file"]=[photsys]
    br["eta_reimers"]= str(eta_reimers)
    br["isoc_val"]= ["1"]
    br["isoc_zeta0"]= str(z)
    br["isoc_lage0"]=str(ages[0])
    br["isoc_lage1"]=str(ages[1])
    br["isoc_dlage"]=str(ages[2])
    br.form.find_control(name='output_gzip').items[0].selected = True
    br.submit()
    link=br.find_link()
    filename= link.text
    os.system('wget -q http://stev.oapd.inaf.it/~lgirardi/tmp/%s -O %s' % (filename,savefilename))
    return None

def get_parsec_isochrones():
    #Setup
    wait= 60. #s
    zs= numpy.arange(0.0005,0.06005,0.0005)
    ages= [6.6,10.1,0.05] #upper has to be <= 10.13
    photsys= "tab_mag_odfnew/tab_mag_sloan_2mass.dat"
    eta_reimers= 0.2
    if eta_reimers == 0.2:
        savedir= os.path.join(os.getenv('ISODIST_DATA'),
                              'parsec-sdss-2mass')
        basefilename= 'parsec-sdss-2mass-Z-'
    else:
        savedir= os.path.join(os.getenv('ISODIST_DATA'),
                              'parsec-%.1f-sdss-2mass' % eta_reimers)
        basefilename= 'parsec-%.1f-sdss-2mass-Z-' % eta_reimers
    for z in zs:
        sys.stdout.write('\r'+"Working on Z = %.4f ...\r" % z)
        sys.stdout.flush()
        savefilename= os.path.join(savedir,
                                   basefilename+'%.4f.dat.gz' % z)
        if os.path.exists(savefilename): continue
        get_parsec_isochrone_onez(z,
                                  savefilename,
                                  ages=ages,
                                  photsys=photsys,
                                  eta_reimers=eta_reimers)
        if wait != 0:
            time.sleep(wait)
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    return None

if __name__ == '__main__':
    get_parsec_isochrones()
