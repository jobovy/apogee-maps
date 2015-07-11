###############################################################################
# check-abundanceScatter-clusters.py: check the precision of the abundances
#                                     used to define the MAPs, by looking
#                                     at internal scatter in open clusters
###############################################################################
import numpy
import read_clusterdata
import define_rcsample
from galpy.util import bovy_plot
# Minimum and maximum Teff to consider
_TEFFMIN= 4500.
_TEFFMAX= 5200.
def checkAbundanceScatterClusters():
    # First read the cluster data
    cldata= read_clusterdata.read_caldata()
    # Read the allStar data to match
    # For each of the calibration open clusters, calculate the offset from the 
    # mean in our FEHTAG and AFETAG
    clusters= ['M71','N2158','N2420','N188','M67','N7789','N6819',
               'N6791']
    fehoffset= []
    afeoffset= []
    for cluster in clusters:
        tdata= cldata[cldata['CLUSTER'] == cluster.upper()]
        tdata= tdata[(tdata['TEFF'] < _TEFFMAX)\
                         *(tdata['TEFF'] > _TEFFMIN)\
                         *(tdata['LOGG'] < 3.5)]
        # Compute the average feh and afe and save the offsets
        medianfeh= numpy.median(tdata['FE_H'])
        medianafe= numpy.median(tdata[define_rcsample._AFETAG])
        fehoffset.extend(tdata['FE_H']-medianfeh)
        afeoffset.extend(tdata[define_rcsample._AFETAG]-medianafe)
        if cluster == 'M67': print medianfeh, medianafe, len(tdata)
    fehoffset= numpy.array(fehoffset)
    afeoffset= numpy.array(afeoffset)
    print 'FE_H scatter %g' % (numpy.nanstd(fehoffset[numpy.fabs(fehoffset) < 0.3]))
    print 'A_FE scatter %g' % (numpy.nanstd(afeoffset[numpy.fabs(afeoffset) < 0.3]))
    gindx= (numpy.fabs(fehoffset) < 0.3)*(numpy.fabs(afeoffset) < 0.3)
    print 'FE_H/A_FE correlation %g' % (numpy.mean(afeoffset[gindx]*fehoffset[gindx])/numpy.nanstd(fehoffset[numpy.fabs(fehoffset) < 0.3])/numpy.nanstd(afeoffset[numpy.fabs(afeoffset) < 0.3]))
    print 'FE_H robust scatter %g' % (1.4826*numpy.median(numpy.fabs(fehoffset)))
    print 'A_FE robust scatter %g' % (1.4826*numpy.median(numpy.fabs(afeoffset)))
    bovy_plot.bovy_print()
    bovy_plot.bovy_hist(fehoffset,range=[-0.3,0.3],bins=31,histtype='step')
    bovy_plot.bovy_hist(afeoffset,range=[-0.3,0.3],bins=31,histtype='step',
                        overplot=True)
    bovy_plot.bovy_end_print('test.png')
    return None

if __name__ == '__main__':
    checkAbundanceScatterClusters()
