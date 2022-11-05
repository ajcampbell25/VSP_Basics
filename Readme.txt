VSP processing in jupyter notebooks

So far I have:

1. Generate spectra of various types, look at wavelets
2. Process a VSP using median filters to separate wavefields
3. Deconvolve a VSP using traditional down-wave decon. This notebook uses the outputs from
   the median filter notebook.
4. A simple AVO modeling notebook. It covers the basics in a basic way!
5. A Q estimation notebook, using spectral ratios to estimate Q from a ZVSP.
6. An anisotropy notebook in which I look at the simpler methods of estimating anisotropy from 
   walkaway VSPs.
   
   
iovsp folder has segy reading and writing modules
plotvsp folder has trace plotting modules
procvsp folder has some basic seismic tool modules
nb_images folder has graphics files for the notebook markdown sections
data folder contains 4 VSP data sets, 3 of them are finite difference synthetics from a fictional model