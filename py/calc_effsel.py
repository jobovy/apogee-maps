###############################################################################
# calc_effsel: calculate the effective selection function for various setups, 
#              **includes the area of the plate**
###############################################################################
import os, os.path
import pickle
import multiprocessing
import numpy
from optparse import OptionParser
from galpy.util import save_pickles, multi
import apogee.select.apogeeSelect
import mwdust
from define_rcsample import get_rcsample
def calc_effsel(args,options,sample=None):
    # Work-horse function to compute the effective selection function, 
    # sample is a data sample of stars to consider for the (JK,Z) sampling
    # Setup selection function
    selectFile= '../savs/selfunc-nospdata.sav'
    if os.path.exists(selectFile):
        with open(selectFile,'rb') as savefile:
            apo= pickle.load(savefile)
    else:
        # Setup selection function
        apo= apogee.select.apogeeSelect()
        # Delete these because they're big and we don't need them
        del apo._specdata
        del apo._photdata
        save_pickles(selectFile,apo)
    # Get the full data sample for the locations (need all locations where 
    # stars could be observed, so the whole sample, not just the subsample
    # being analyzed)
    data= get_rcsample()
    locations= list(set(list(data['LOCATION_ID'])))
    # Load the dust map and setup the effective selection function
    if options.dmap.lower() == 'green15':
        dmap3d= mwdust.Green15(filter='2MASS H')
    elif options.dmap.lower() == 'marshall06':
        dmap3d= mwdust.Marshall06(filter='2MASS H')
    elif options.dmap.lower() == 'drimmel03':
        dmap3d= mwdust.Drimmel03(filter='2MASS H')
    elif options.dmap.lower() == 'sale14':
        dmap3d= mwdust.Sale14(filter='2MASS H')
    elif options.dmap.lower() == 'zero':
        dmap3d= mwdust.Zero(filter='2MASS H')
    # Sample the M_H distribution
    if options.samplemh:
        if sample is None: sample= data
        MH= sample['H0']-sample['RC_DM']
        MH= numpy.random.permutation(MH)[:1000] # do 1,000 max
    else:
        MH= -1.49
    apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=dmap3d,MH=MH)
    # Distances at which to calculate the effective selection function
    distmods= numpy.linspace(options.dm_min,options.dm_max,options.ndm)
    ds= 10.**(distmods/5-2.)
    # Now compute all selection functions
    out= multi.parallel_map((lambda x: _calc_effsel_onelocation(\
                locations[x],apof,apo,ds)),
                            range(len(locations)),
                            numcores=numpy.amin([len(locations),
                                                 multiprocessing.cpu_count(),options.multi]))
    # Save out
    out= numpy.array(out)
    save_pickles(args[0],locations,out,distmods,ds)
    return None

def _calc_effsel_onelocation(locid,apof,apo,ds):
    # Calculate the effective selection function for a given location
    try:
        esf= apof(locid,ds)*apo.area(locid)
    except (IndexError, TypeError,ValueError):
        esf= -numpy.ones_like(ds)
    return esf

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the effective selection function will be saved to"
    parser = OptionParser(usage=usage)
    # Distances at which to calculate the effective selection function
    parser.add_option("--dm_min",dest='dm_min',default=7.,type='float',
                      help="Minimum distance modulus")
    parser.add_option("--dm_max",dest='dm_max',default=15.5,type='float',
                      help="Maximum distance modulus")
    parser.add_option("--ndm",dest='ndm',default=301,type='int',
                      help="Number of distance moduli to calculate the function at")
    # Dust map to use
    parser.add_option("--dmap",dest='dmap',default='green15',
                      help="Dust map to use ('Green15', 'Marshall03', 'Drimmel03', 'Sale14', or 'zero'")
    # Sample over the M_H of the sample?
    parser.add_option("--samplemh",action="store_true", dest="samplemh",
                      default=False,
                      help="If set, sample the M_H distribution of the sub-sample (default= full sample)")
    # Multiprocessing?
    parser.add_option("-m","--multi",dest='multi',default=1,type='int',
                      help="number of cpus to use")
    return parser

if __name__ == '__main__':
    parser= get_options()
    options, args= parser.parse_args()
    calc_effsel(args,options)
