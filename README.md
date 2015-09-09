# nir_rv
Radial velocities for low-resolution NIR spectra using telluric features to provide absolute wavelength calibration. Requires a standard star (M dwarf spectra included), and both the telluric and non-telluric-corrected spectra. Orders do not need to be combined. 

See check_standards for an example of usage. 

Run check_standards to test code.

Choose cross-correlation routine. xcorl is not my code, but may be available to you. c_correlate should exist in any modern IDL distribution. I prefer the routine cross_correlate which may be available in your installation. xcorl and cross_correlate produce consistent results; values from c_correlate may differ by up to 1 km/s.

If you use this code in published research, please cite Newton et al. (2014): http://adslabs.org/adsabs/abs/2014AJ....147...20N/
