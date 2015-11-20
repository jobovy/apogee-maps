apogee-maps
============

This repository contains the code for the papers "The stellar
population structure of the Galactic disk" `Bovy, Rix, Schlafly et
al. (2015) <http://arxiv.org/abs/1509.05796>`__ and "On Galactic
density modeling in the presence of dust extinction" `Bovy, Rix, Green
et al. (2015) <http://arxiv.org/abs/1509.06751>`__. This code is
provided for illustration purposes only, but it should allow you to
perform the same analysis as in these papers if you want. Much of this
code runs for quite a long time and no guidance is given here as to
how long each command takes (none of it takes *ridiculously* long
though).

Prerequisites
--------------

This code heavily depends on the `apogee
<https://github.com/jobovy/apogee>`__ and `mwdust
<https://github.com/jobovy/mwdust>`__ general-use packages. You need to
install these (including downloading all of the necessary dust data)
first. This code also makes use of the numpy, scipy, matplotlib
stack. Parts of the code require healpy. The code will complain if you
do not have other required packages.

The stellar population structure of the Galactic disk
-------------------------------------------------------

This code is mainly concerned with fitting spatial-density profiles to
subsets of APOGEE data. Almost all of the commands below require the
effective-selection function (see the papers above for a description
of this), which you should first calculate as::

   python calc_effsel.py ../essf/essf_green15.sav --dmap=green15 -m 32
   python calc_effsel.py ../essf/essf_marshall06.sav --dmap=marshall06 -m 32
   python calc_effsel.py ../essf/essf_sale14.sav --dmap=sale14 -m 32
   python calc_effsel.py ../essf/essf_drimmel03.sav --dmap=drimmel03 -m 32
   python calc_effsel.py ../essf/essf_zero.sav --dmap=zero -m 32

The first two of these are the most important, as they are used in
most parts of the code. The option ``-m 32`` calculates this in
parallel (for different lines of sight) on 32 CPUs. There is some
overhead in first computing and saving the overall APOGEE selection
function (see ``apogee.select.apogeeSelect``), so you want to run
these in series (at least the first two).

The heart of the code are the ``densprofiles.py`` and ``fitDens.py``
files. These contain all of the different density profiles and the
code to fit and MCMC-explore a given density profile to a subset of
APOGEE data.

The notebook `SpatialFitExample
<https://github.com/jobovy/apogee-maps/blob/master/py/SpatialFitExample.ipynb>`__
has general explorations of the fit provided by different types of
density profiles. It formed the backbone of the exploratory phase of
this project.

The following sections list the commands that generate the figures in
the paper.

Data
++++

Figure 1::

    python plot_dustwapogee.py ./dustwapogee.ps

Figure 2::

       python plot_vs_jkz.py -o ./rcmodel_mode_jkz_h_parsec.ps rcmodel_mode_jkz_h_parsec_newlogg.sav -b H -t mode --parsec

(this code is in `jobovy/apogee-rc
<https://github.com/jobovy/apogee-rc>`__ as it was part of an earlier
paper).

Figure 3::

       python plot_rcdistancecomparison.py ./rcdistcomp.ps

Figure 4::

       python plot_afe_spectra.py afe_feh-0.3_residuals.sav ./spectra_ELEM.ps

Figure 5::

       python plot_afefeh.py ./overall_afefeh.ps 

Figure 6::

       python plot_spatial_broad.py ./spatial_broad_SAMPLE.ps
       
Checking the abundance scatter of open clusters::

	 python check-abundanceScatter-clusters.py


Fits to broad abundance-selected subsamples
++++++++++++++++++++++++++++++++++++++++++++

Figure 7::

       python fitBroadSubsamples.py lowlow ../broadfits/lowlow.sav ./lowlow.txt ./lowlow_LOC.ps
       python fitBroadSubsamples.py solar ../broadfits/solar.sav ./solar.txt ./solar_LOC.ps
       python fitBroadSubsamples.py highfeh ../broadfits/highfeh.sav ./highfeh.txt ./highfeh_LOC.ps       

Figure 8::

       python fitBroadSubsamples.py highalpha ../broadfits/highalpha.sav ./highalpha.txt ./highalpha_LOC.ps 

Figure 9::

       python plot_broadsurfdens.py ./broadfits-radial.ps

Fits to mono-abundance populations (MAPs)
++++++++++++++++++++++++++++++++++++++++++

Fits are done as follows::

     python fitMAPs.py tribrokenexp ../mapfits/tribrokenexp.sav
     python fitMAPs.py tribrokenexpfixedflare ../mapfits/tribrokenexpfixedflare.sav
     python fitMAPs.py tribrokenexpflare ../mapfits/tribrokenexpflare.sav
     python fitMAPs.py tribrokenexpinvlinflare ../mapfits/tribrokenexpinvlinflare.sav
     python fitMAPs.py tribrokenexplinflare ../mapfits/tribrokenexplinflare.sav
     python fitMAPs.py tribrokentwoexp ../mapfits/tribrokentwoexp.sav

that is, ``fitMAPs.py`` just fits a given density profile to all
MAPs. It also runs MCMC for each MAP and saves all of these results
(in a rather large file!).

Figure 10::

       python plot_maprmax.py ../mapfits/tribrokenexpflare.sav ./mapfits-rpeak.png

Figure 11::
       
       python plot_mapsurfdens.py ./mapfits-radial.ps
       python plot_mapsurfdens_highalpha.py ./mapfits-radial-highalpha.ps 

The resulting output profiles for the low- and high-alpha MAPs displayed in this figure can be found under `<csv/mapsurfdens.csv>`__ and `<csv/mapsurfdens_highalpha.csv>`__.

Figure 12::

       python plot_mapflarepdf.py flare_lowalpha.sav ./mapfits-flare-lowalpha.ps
       python plot_mapflarepdf.py flare_highalpha.sav ./mapfits-flare-highalpha.ps

Figure 13::

       python plot_mapflare.py ./mapfits-radialflare.ps
       python plot_mapflare_highalpha.py ./mapfits-radialflare-highalpha.ps 

The resulting output profiles for the low- and high-alpha MAPs displayed in this figure can be found under `<csv/mapflare.csv>`__ and `<csv/mapflare_highalpha.csv>`__.

Figure 14::

       python plot_maphz.py ./mapfits-hz.png

Figure 15::

       python plot_maptwohz.py ./mapfits-twohz.ps


On Galactic density modeling in the presence of dust extinction
-----------------------------------------------------------------

Many fewer figures in this paper (phew!), but here we go. Figure 1::

     python plot_dust_gaia.py 5.0 ./dust_5.0kpc.ps 

Figure 2::

       python plot_gaia_rcmag.py ./gaia_mg.ps

(see ``gaia-rc.py`` for some code to get the RC's properties in the
Gaia passband). Figure 3::

     python plot_powspec.py 5.0 ../savs/PowspecDensAndDustAndESSF ./powspec_dens_dust_essf_D5.0.ps
     python plot_powspec.py 6.3 ../savs/PowspecDensAndDustAndESSF ./powspec_dens_dust_essf_D6.3.ps

Figure 4::

       python plot_distanceintegral.py ../savs/distInt.sav /dev/null

(and similar for subsets of the sky, see options in
``plot_distanceintegral.py``).

Figure 5::

       python plot_ah_location.py 4240 ./ah_4240.png

and similar for other locations (like 4240). Figure 6::

    python plot_effsel_location.py 4240 ./effsel_4240.ps 

also similar for other locations.