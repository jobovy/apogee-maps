import numpy as np
import define_rgbsample
import calc_effsel as ce
from optparse import OptionParser

class EffselOptions(object):
	def __init__(self, dm_min = 7.,
				 dm_max = 15.5,
				 ndm = 301,
				 dmap = 'marshall06',
				 samplemh = False,
				 m = 32,
				 dr = '12',
				 samp = 'rc',
				 modelmh = False,
				 distgrid = False):
		self.dm_min = dm_min
		self.dm_max = dm_max
		self.ndm = ndm
		self.dmap = dmap
		self.samplemh = samplemh
		self.multi = m
		self.dr = dr
		self.samp = samp
		self.modelmh = modelmh
		self.distgrid = distgrid

opt = EffselOptions(dmap = 'green15', samplemh = True, samp='rgb', modelmh=True, distgrid=False)

def effsel_bins(opt):
	agebins = np.arange(1.,14.,2.)
	fehbins = np.arange(-0.6,0.3,0.1)
	for i in range(0, len(fehbins)-1):
		for j in range(0,len(agebins)-1):
			print 'Calculating Effective Selection function for [Fe/H] = '+str(round(fehbins[i],1))+' and age = '+str(round(agebins[j],1))+''
			if opt.distgrid == False:
				filename = '../essf/maps/essf_'+opt.samp+'_'+opt.dmap+'_modelmh_feh'+str(round(fehbins[i],1))+'_age'+str(round(agebins[j],1))+'.sav'
			if opt.distgrid == True:
				filename = '../essf/maps/essf_'+opt.samp+'_'+opt.dmap+'_distgrid_modelmh_feh'+str(round(fehbins[i],1))+'.sav'   
			ce.calc_effsel([filename,], opt, fehbin=[fehbins[i],fehbins[i+1]], agebin=[agebins[j], agebins[j+1]])


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
                      default=True,
                      help="If set, sample the M_H distribution of the sub-sample (default= full sample)")
    # Multiprocessing?
    parser.add_option("-m","--multi",dest='multi',default=1,type='int',
                      help="number of cpus to use")
    # Alternate Data releases?
    parser.add_option("--dr",dest='dr',default=None,
    				  help="Data release to use")
    # RGB or RC?
    parser.add_option("--samp", dest='samp', default='rgb', help = "rc or rgb")
    # model M_H distribution?
    parser.add_option("--modelmh", dest='modelmh', default=True, help="If True, sample a model M_H distribution from parsec isochrones.")
    return parser

if __name__ == '__main__':
    parser= get_options()
    options, args = parser.parse_args()
    effsel_bins(options)