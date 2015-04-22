###############################################################################
# triangleMAPs.py: make triangle plots for all of the MAPs
###############################################################################
import sys
import pickle
import numpy
from galpy.util import bovy_plot
import triangle
import define_rcsample
def triangleMAPs(savefilename,basename):
    with open(savefilename,'rb') as savefile:
        bf= numpy.array(pickle.load(savefile))
        samples= numpy.array(pickle.load(savefile))
        bf_g15= numpy.array(pickle.load(savefile))
        samples_g15= numpy.array(pickle.load(savefile))
        bf_zero= numpy.array(pickle.load(savefile))
        samples_zero= numpy.array(pickle.load(savefile))
    labels= []
    for jj in range(samples.shape[2]):
        labels.append(r"$\mathrm{param}\ %i$" % jj)
    maps= define_rcsample.MAPs()
    for ii, map in enumerate(maps.map()):
        if ii >= len(bf): break
        tfeh= numpy.nanmedian(map['FE_H'])
        tafe= numpy.nanmedian(map[define_rcsample._AFETAG])
        for tbf,tsamples,ext in zip([bf,bf_g15,bf_zero],
                                    [samples,samples_g15,samples_zero],
                                    ['fid','g15','zero']):
            triangle.corner(tsamples[ii,].T,quantiles=[0.16, 0.5, 0.84],
                            labels=labels,
                            show_titles=True,title_args={"fontsize": 12},
                            bins=21)
            bovy_plot.bovy_text(r'$[\mathrm{{Fe/H}}] = {feh:.1f},$'\
                                    .format(feh=tfeh)+'\n'
                                +r'$[\alpha/\mathrm{{Fe}}] = {afe:.2f}$'\
                                    .format(afe=tafe),
                                top_left=True,size=16.)
            bovy_plot.bovy_end_print(basename+"_%i_%s.png" % (ii,ext))
    return None

if __name__ == '__main__':
    triangleMAPs(sys.argv[1],sys.argv[2])
