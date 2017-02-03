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
import isodist
from define_rcsample import get_rcsample
from define_rgbsample import get_rgbsample
from scipy.interpolate import interp1d

def calc_effsel(args,options,sample=None, fehbin=[0.,0.1], agebin=[0.,13.]):
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
    if options.samp == 'rc':
    	data= get_rcsample(dr = options.dr)
    if options.samp == 'rgb':
    	data= get_rgbsample(add_ages = True)
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
        if options.samp == 'rc':
        	MH= sample['H0']-sample['RC_DM']
        	MH= numpy.random.permutation(MH)[:1000] # do 1,000 max
        if options.samp == 'rgb':
        	if options.modelmh == True:
        		MH = sample_iso_MH(fehbin = fehbin, n=1000)
        	if options.modelmh != True:
        		MH = numpy.random.permutation(sample['M_H'])[:1000]
    else:
        MH= -1.49
    apof= apogee.select.apogeeEffectiveSelect(apo,dmap3d=dmap3d,MH=MH)
    # Distances at which to calculate the effective selection function
    if options.distgrid == False:
        distmods= numpy.linspace(options.dm_min,options.dm_max,options.ndm)
        ds= 10.**(distmods/5-2.)
    if options.distgrid == True:
        ds = numpy.linspace(0.01,20.,300)
        distmods= 5*numpy.log(ds)-2.
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
    
def sample_iso_MH(fehbin = [-0.6, 0.2], agebin=[0.,13.], n=1000, agebias=False, imftype='chabrier2003', isochrones='Padova'):
	#Load pre-computed parsec isochrone grid
	iso_file = open('../savs/'+isochrones+'_grid_'+imftype+'.sav')
	iso_grid = pickle.load(iso_file)
	#add rgb age bias if agebias = True
	if agebias == True:
		iso_grid[:,6] *= agebias(10**iso_grid[:,0])
	#Perform sample cuts (logg and J-K)
	gridcuts = (iso_grid[:,3] > 1.8)&(iso_grid[:,3] < 3.0)&(iso_grid[:,5] > 0.5)
	cutgrid = iso_grid[gridcuts]
	n_weights = cutgrid[:,6]*(10**cutgrid[:,0]/cutgrid[:,1])
	#make [Fe/H], age cut
	fehcut = (isodist.Z2FEH(cutgrid[:,1])>=fehbin[0])&(isodist.Z2FEH(cutgrid[:,1])<fehbin[1])\
			&(10**cutgrid[:,0] >= agebin[0])&(10**cutgrid[:,0] < agebin[1])
	# compute CDF of M_H
	sorter_H = numpy.argsort(cutgrid[:,4][fehcut])
	cdf_H = numpy.cumsum(n_weights[fehcut][sorter_H])/numpy.sum(n_weights[fehcut])
	#Interpolate CDF and take n samples
	intercdf_H = interp1d(cdf_H, cutgrid[:,4][fehcut][sorter_H])
	rand = numpy.random.uniform(size=n, low=0.0001, high=0.999999)
	model_MH = intercdf_H(rand)
	return model_MH
	

def agebias(age):
    return (age+1)**-0.7

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
    parser.add_option("--distancegrid",action="store_true", dest="distgrid",
                      default=False,
                      help="if set, use distance grid rather than distmod")
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
    # Alternate Data releases?
    parser.add_option("--dr",dest='dr',default=None,
    				  help="Data release to use")
    # RGB or RC?
    parser.add_option("--samp", dest='samp', default='rc', help = "rc or rgb")
    # model M_H distribution?
    parser.add_option("--modelmh", dest='modelmh', default=False, help="If True, sample a model M_H distribution from parsec isochrones.")
    return parser

if __name__ == '__main__':
    parser= get_options()
    options, args= parser.parse_args()
    calc_effsel(args,options)
