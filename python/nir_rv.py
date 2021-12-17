import numpy as np
from astropy.io import fits
from get_helio import GET_HELIO
from tell_model import TELL_MODEL
from ern_rv import ERN_RV

def DO_RVSHIFTS(rv0, hdr, shdr, rv_std=None, quiet=0):
    """
    Get and apply heliocentric correction for both the object and standard.
    Only applies if the standard spectrum is not already in the rest frame.
    Args:
        rv0: Relative RV between the object and standard.
        hdr: FITS header of the object.
        shdr: FITS header of the standard.
        rv_std: Radial velocity of the standard.
        quiet: 0 for verbose, 1 for quiet output.
    Returns:
        rv: Absolute radial velocity of the object.
        hcorr: Heliocentric correction for the object.
    """

    # barycentric correction
    hcorr = GET_HELIO(hdr, quiet=quiet)

    # offset from standard star
    hcorr_std = GET_HELIO(shdr, quiet=quiet)
    if rv_std is None:
        rv_std = 0
    rvoff = rv_std - hcorr_std

    # heliocentric corrected absolute RV
    rv = rv0 + hcorr + rvoff

    return rv, hcorr

def NIR_RV(mydata_tc, hdr, mydata, mystd_tc, shdr, mystd=None,
           wlcal=None, atrest=False, stdrv=None,
           atrans=None,
           plorder=None, pixscale=None,
           s_plorder=None, spixscale=None,
           trange=0, wrange=None,
           maxshft=None,
           shift_algo=None, quiet=0):

    """
    Calculates the radial velocity of the object.
    Args:
        mydata_tc: Telluric-corrected object data.
        hdr: FITS header for object data.
        mydata: Non-telluric-corrected object data.
        mystd_tc: Tellucir-corrected standard data.
        mystd: Non-telluric-corrected standard data.
        wlcal: 1 if the standard spectrum has already been telluric-caliburated, 0 if it has not.
        atrest: 1 if the standard spectrum is already shifted to rest wavelengths, 0 if it is not.
        sdtrv: Known RV of the standard star.
        atrans: Path to atmopsheric transmission spectrum file.
        plorder: Order of Legendre polynomial for the object in TELL_MODEL.
        s_plorder: Order of Legendre polynomial for the standard in TELL_MODEL.
        trange: Max and min wavelength to interpolate over in TELL_MODEL.
        wrange: Max and min wavelength to interpolate over in ERN_RV.
        maxshft: Maximum allowable wavelength shift (in microns) between the atmospheric transmission spectrum
                 and the object/standard in TELL_MODEL.
        shift_algo: Algorithm to use when calculating the RV shift between the object and standard in ERN_RV.
                    "LM" for Levenberg-Marquardt fit.
                    "CC" for cross-correlation.
        quiet: 0 for verbose, 1 for quiet output.
    Returns:
        rv: Absolute radial velocity of the object.
        rv0: Radial velocity of the object relative to the standard.
        torest: Radial velocity correction necessary to shift the object to the rest frame.
        mshft: Wavelength shift required to match the object with the telluric spectrum.
    """

    if pixscale is None:
        pixscale = np.abs(np.median(mydata_tc[0:-2,0]-mydata_tc[1:-1,0]))
    if spixscale is None:
        spixscale = np.abs(np.median(mystd_tc[0:-2,0]-mystd_tc[1:-1,0]))
    if plorder is None:
        plorder = 5
    if s_plorder is None:
        s_plorder = 5
    if shift_algo == None:
        shift_algo = "LM"

    data = mydata
    data_tc = mydata_tc
    if mystd is not None:
        std = mystd
    std_tc = mystd_tc

    if atrans is None:
        if quiet == 0:
            print("NIR_RV: Reading atrans from spex directory.")
            with fits.open('../atrans.fits') as hdul:
                atrans = hdul[0].data
    elif type(atrans) is str:
        if quiet == 0:
            print("NIR_RV: Reading atrans from file name supplied.")
            with fits.open(atrans) as hdul:
                atrans = hdul[0].data
    else:
        if quiet == 0:
            print("NIR_RV: atrans array supplied.")

    # parameter for oversampling the spectrum; value set by minimum needed for spex echelle
    oversamp=6.0

    # maximum shift to consider; in spex echelle, no shifts beyond 0.0006 (and very few >0.0004)
    if maxshft is None:
        maxshft=0.0008
    if quiet == 0:
        print("NIR_RV: max shift in microns is ", maxshft)

    #
    # SCIENCE TARGET
    #

    # data: get wavelength calibration by modeling the telluric features
    data_new, atrans_new, mshft = TELL_MODEL(atrans, data, plorder, trange, maxshft, oversamp, pixscale)

    # shift telluric corrected data to absolute wavelength
    data_tc_new = data_tc
    data_tc_new[0] = data_tc[0]+mshft

    #
    # STANDARD STAR
    #

    # standard: get wavelength calibration by modeling the telluric features
    if wlcal is None:
        maxshft=maxshft/s_pixscale*oversamp
        std_new, atrans_new, s_shft = TELL_MODEL(atrans, std, plorder, trange, maxshft, oversamp, pixscale)
        # shift telluric corrected standard to absolute wavelength
        std_tc_new = std_tc
        std_tc_new[0] = std_tc[0]+s_mshft
        if quiet == 0:
            print("NIR_RV: Wavecal done for standard.")
    else: # user says this has already been done
        std_tc_new = std_tc
        if quiet == 0:
            print("NIR_RV: wlcal keyword set. No calibration done.")

    #
    # GET RADIAL VELOCITY
    #

    # now cross-correlate data and standard to get relative radial velocity
    rv0 = ERN_RV(data_tc_new, std_tc_new, wrange=wrange, pixscale=pixscale, shift_algo=shift_algo, quiet=quiet)

    if atrest is True:
        bc = GET_HELIO(hdr, quiet=quiet)
        rv = rv0 + bc.value # actual RV
        torest = rv0 # the required correction
        if quiet == 0:
            print("NIR_RV: atrest keyword set. No RV for standard star.")
    else:
        rv, bc = DO_RVSHIFTS(rv0, hdr, shdr, rv_std=stdrv, quiet=quiet)
        torest = rv - bc
        if quiet == 0:
            print("NIR_RV: Adjusting for RV of standard star using velocity provided and barycentric velocity.")

    shft = mshft - (mshft+data_tc[0])*torest/(3e5) # this is what you add to get to rest wavelengths! Fixed to treat linear shift correctly.

    return rv, rv0, torest, mshft