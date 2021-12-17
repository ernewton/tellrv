import numpy as np
from astropy.io import fits
from order_variables import order_variables
from nir_rv import NIR_RV

def check_standards(method = 3):

    quiet = 1

    # read in atmospheric transmission spectrum (Lord, 1992)
    # file is big, best to read it in only once
    # you can find atrans in the spextool libary at Spextool/data/atrans.fits
    with fits.open('../atrans.fits') as hdul:
        atrans = hdul[0].data

    # list of standard stars and their RVs
    standards = np.loadtxt("../standards.txt",
                           dtype={"names": ("stars", "sptypes", "rvs"),
                                  "formats": ("U10", int, float)})
    stars, sptypes, rvs = standards["stars"], standards["sptypes"], standards["rvs"]
    spec_dir = "../spec/"

    # RV standard: non-telluric corrected spectrum
    # absolute wavelength calibration is set by the telluric features

    # RV standard: telluric-corrected spectrum
    # cross-correlation with standard star is done on final data product

    # METHOD 1: use the standard spectrum as-is
    if method == 1:
        file = spec_dir+stars[0]+".fits"
        with fits.open(file) as hdul:
            std0 = hdul[0].data
            shdr = hdul[0].header
        file_tc = spec_dir+stars[0]+"_tc.fits"
        with fits.open(file_tc) as hdul:
            std_tc0 = hdul[0].data
        std_tc0[:,0,:]=std0[:,0,:]  # want to use original wavelength array
        stdrv = rvs[0]

    # METHOD 2: use the wavelength-calibrated standard star
    elif method == 2:
        wlcal = 1
        atrest = False
        stdrv = rvs[0]
        file_tc = spec_dir+'J0727+0513_wlcal.fits'
        with fits.open(file_tc) as hdul:
            std_tc0 = hdul[0].data
            shdr = hdul[0].header
        std0 = std_tc0

    # METHOD 3: use the at-rest standard star created by makemeastandwich
    elif method == 3:
        wlcal = 1
        atrest = True
        stdrv = rvs[0]
        file_tc = spec_dir+'J0727+0513_rest.fits'
        with fits.open(file_tc) as hdul:
            std_tc0 = hdul[0].data
            shdr = hdul[0].header
        std0 = std_tc0

    for i in range(len(stars)):
        std = std0
        std_tc = std_tc0

        # read in data for this star
        star = stars[i]
        file = spec_dir+star+".fits"
        with fits.open(file) as hdul:
            data = hdul[0].data
            hdr = hdul[0].header
        file_tc = spec_dir+star+"_tc.fits"
        with fits.open(file_tc) as hdul:
            data_tc = hdul[0].data

        print("===========")
        print("Star is: ", star)
        print("File is: ", file)

        data_tc[:,0,:]=data[:,0,:] # wavelength was modified during reduction

        # J, H, and K bands independently
        temp = np.zeros(3)
        for order in range(3):
            # get pertinant variables for SpeX
            wrange, trange, pixscale, plorder = order_variables(hdr, order, instrument="spex")
            # calculate the RV
            myrv, x1, x2, x3 = NIR_RV(data_tc[order], hdr, data[order],
                               std_tc[order], shdr, std[order],
                               wlcal=wlcal, atrest=atrest, stdrv=stdrv, # already wavelength calibrated?
                               atrans=atrans,
                               pixscale=pixscale, plorder=plorder,
                               spixscale=pixscale, s_plorder=plorder,  # standard is from same set-up
                               wrange=wrange, trange=trange,
                               quiet=quiet)
            if (order < 3):
                temp[order] = myrv

        # Newton et al. (2014) adopts the median of the first three orders as the RV
        print("RVs from K, H, J bands:", temp)
        print("RV from NIR spectrum is:", np.median(temp))
        print("Actual RV is:           ", rvs[i])
        print("===========")

    print("NOTE: RVs from Newton et al. (2014) include a -2.6 km/s systematic RV correction, not included here.")

check_standards()