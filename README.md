# Absolute near-infared radial velocities

This code measures absolute radial velocities for low-resolution NIR spectra. I use telluric features to provide absolute wavelength calibration, and then cross-correlate with a standard star. To use the code for your star, you will need observations of a standard star (included) and your science target (examples included). You will need both the telluric and non-telluric-corrected spectra and the FITS headers, but no other spectra is required. Orders do not need to be combined. It should work for all extant spectra taken with SpeX on IRTF or similar instruments. 

Please see nir_rv and check_standards for an example of usage. 

Run check_standards to test code.

## Altenrative sub-routines

### Cross-correlation
There are three options for cross-correlation routines: xcorl, c_correlate, and cross_correlate. These are used in ern_rv, and can be selected by keyword in the top-level routine nir_rv. c_correlate is the default since I believe it is the most commonly available, and ought to be included exist in any modern IDL distribution. xcorl is not my code, but may be available to you. I prefer the routine cross_correlate which may be available in your installation. xcorl and cross_correlate produce consistent results; I have improved the peak ginding in c_correlate but values may still differ by up to 1 km/s.

### Continuum fitting
There are two options for continuum fitting: a spline-based routine and the routine contf. contf is not my code, but may be available to you. The default spline-based routine is based on contf and is included in this code. The choice is not important in most cases.

# Reference
If you use this code in published research, please cite Newton et al. (2014): http://adslabs.org/adsabs/abs/2014AJ....147...20N/

[![DOI](https://zenodo.org/badge/4705/ernewton/tellrv.svg)](https://zenodo.org/badge/latestdoi/4705/ernewton/tellrv)

# License
Copyright 2015 Elisabeth R. Newton. Licensed under the terms of the MIT License (see LICENSE).
