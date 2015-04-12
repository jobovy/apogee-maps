import numpy
import apogee.tools.read as apread
import astropy.io.ascii
import define_rcsample
def read_caldata(filename='../cldata/aj485195t4_mrt.txt'):
    data= astropy.io.ascii.read(filename)
    data.rename_column('Cluster','CLUSTER')
    data.remove_column('Teff')
    data.rename_column('TeffC','TEFF')
    data.remove_column('logg')
    data.rename_column('loggC','LOGG')
    data.remove_column('[M/H]')
    data.rename_column('[M/H]C','FEH')
    data.rename_column('2MASS','ID')
    # Now match to allStar to get the location_ids, and abundances
    alldata= apread.allStar(raw=True)
    locids= numpy.empty(len(data),dtype='int')
    fehs= numpy.empty(len(data),dtype='float')
    afes= numpy.empty(len(data),dtype='float')
    for ii in range(len(data)):
        if 'Pleiades' in data['CLUSTER'][ii]: continue
        indx= alldata['APOGEE_ID'] == data['ID'][ii]
        if numpy.sum(indx) == 0:
            raise ValueError('allStar match for %s not found ...' % (data['ID'][ii]))
        if len(list(set(alldata['LOCATION_ID'][indx]))) > 1:
            raise ValueError('Multiple matches found for for %s ...' % (data['ID'][ii]))
        locids[ii]= alldata['LOCATION_ID'][indx][0]
        fehs[ii]= alldata['FE_H'][indx][0]
        afes[ii]= define_rcsample.avg_alphafe(alldata[indx])[0]
    data['LOCATION_ID']= locids
    data['FE_H']= fehs
    data[define_rcsample._AFETAG]= afes
    return data
