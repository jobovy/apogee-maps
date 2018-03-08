import os, os.path
import pexpect
import subprocess
from astropy.io import fits
import numpy as np
from optparse import OptionParser
import sys
import tempfile
from ftplib import FTP
import shutil

_ERASESTR= "                                                                                "

def get_agetables(args,options):
    cat = 'J/MNRAS/456/3655/'
    tab1name = 'table1.dat'
    tab2name = 'table2.dat.gz'
    table2out = 'DR12_martigages_vizier.fits'
    table1out = 'martig_table1.fits'
    out = args[0]
    _download_file_vizier(cat, out, catalogname=tab1name)
    _download_file_vizier(cat, out, catalogname=tab2name)
    subprocess.call(['gunzip', os.path.join(out,tab2name)])
    ztab2name = 'table2.dat'
    tab1_colnames = '2MASS_ID, Teff, Logg, [M/H], [C/M], [N/M], Massin, e_Massin, Massout, Age_in, e_Agein, E_Agein, Age_out'
    tab2_colnames = '2MASS_ID, Teff, Logg, [M/H], [C/M], [N/M], Massout, Age'
    tab1 = np.genfromtxt(os.path.join(out,tab1name), dtype=None, names=tab1_colnames)
    tab2 = np.genfromtxt(os.path.join(out,ztab2name), dtype=None, names=tab2_colnames)
    hdu1 = fits.BinTableHDU.from_columns(tab1)
    hdu2 = fits.BinTableHDU.from_columns(tab2)
    hdu1.writeto(os.path.join(out,table1out))
    hdu2.writeto(os.path.join(out,table2out))


def _download_file_vizier(cat,filePath,catalogname='catalog.dat'):
    '''
    Stolen from Jo Bovy's gaia_tools package!
    '''
    sys.stdout.write('\r'+"Downloading file %s ...\r" \
                         % (os.path.basename(filePath)))
    sys.stdout.flush()
    try:
        # make all intermediate directories
        os.makedirs(os.path.dirname(filePath)) 
    except OSError: pass
    # Safe way of downloading
    downloading= True
    interrupted= False
    file, tmp_savefilename= tempfile.mkstemp()
    os.close(file) #Easier this way
    ntries= 1
    while downloading:
        try:
            ftp= FTP('cdsarc.u-strasbg.fr')
            ftp.login('anonymous', 'test')
            ftp.cwd(os.path.join('pub','cats',cat))
            with open(tmp_savefilename,'wb') as savefile:
                ftp.retrbinary('RETR %s' % catalogname,savefile.write)
            shutil.move(tmp_savefilename,os.path.join(filePath,catalogname))
            downloading= False
            if interrupted:
                raise KeyboardInterrupt
        except:
            raise
            if not downloading: #Assume KeyboardInterrupt
                raise
            elif ntries > _MAX_NTRIES:
                raise IOError('File %s does not appear to exist on the server ...' % (os.path.basename(filePath)))
        finally:
            if os.path.exists(tmp_savefilename):
                os.remove(tmp_savefilename)
        ntries+= 1
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()        
    return None
    



def get_options():
	usage = "usage: %prog [options] <outpath>"
	parser = OptionParser(usage=usage)
	# Distances at which to calculate the effective selection function
	parser.add_option("--convert",dest='convert',default=True ,action='store_true', help="convert to fits?")
	return parser

if __name__ == '__main__':
    parser = get_options()
    options, args= parser.parse_args()
    get_agetables(args,options)
