import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

def LEGPOLY(x, p):
    """
    Legendre polynomial of arbitrary order
    Args:
        x: vector of pixel positions
        p: array of polynomial coefficients
    Returns:
        Legendre polynomial for given coefficients.
    """
    poly = np.polynomial.legendre.Legendre(p, [x[0],x[-1]])
    return poly(x)

def TELL_FUNC_fit(p, lambd, atrans, data, pixscale, oversamp, maxshft):
    """
    Function for fitting atmospheric transmission spectrum to either the object or standard spectrum.
    Input to Levenberg-Marquardt fit in TELL_FUNC.
    Args:
        p: Array containing the model parameters.
           p[0] is the shift in pixels.
           p[1] is the flux scaling factor.
           p[2:] is the Legendre polynomial coefficients.
        lambd: The wavelength vector for the spectra.
        atrans: The flux vector for the asmopsheric transmission spectrum.
        data: The flux vector for the object or standard spectrum.
        pixscale: Microns/pixel of the data.
        oversamp: Oversampling factor.
        maxshft: Maximum number of microns atrans can shift.
    Returns:
        Difference in flux between model (shifted, scaled, and curved atrans) and data at each wavelength.
        To be minimized.
    """

    # scale atrans by a constant to account for precipital water vapor and airmass differnces
    atrans_new=atrans**(p[1])

    # shift wavelengths of atrans spectrum
    shft = p[0]*pixscale/oversamp
    wl_shift = lambd + shft
    atrans_new=np.interp(x=wl_shift, xp=lambd, fp=atrans_new)

    # curve atrans spectrum with Legendre polynomial
    x = np.linspace(-1, 1, len(data))
    poly = LEGPOLY(x, p[2:])
    atrans_curved = atrans_new*poly

    # apply maximum shift condition
    if abs(shft) > maxshft:
        diff = np.ones(len(data))*np.inf
    else:
        diff=(data-atrans_curved)

    return diff

def TELL_FUNC(p, lambd, atrans, data, pixscale, oversamp):
    """
    Function for returning diagnostics relating to the fit above.
    Args:
        p: Array containing the model parameters.
           p[0] is the shift in pixels.
           p[1] is the flux scaling factor.
           p[2:] is the Legendre polynomial coefficients.
        lambd: The wavelength vector for the spectra.
        atrans: The flux vector for the asmopsheric transmission spectrum.
        data: The flux vector for the object or standard spectrum.
        pixscale: Microns/pixel of the data.
        oversamp: Oversampling factor.
    Returns:
        shft: Shift in microns between the two spectra.
        model: Altered atmospheric transmission spectrum.
        cont: Curve generated from Legendre polynomial.
        diff: Difference between model and data at each point.
    """

    # scale atrans by a constant to account for precipital water vapor and airmass differnces
    atrans_new=atrans**(p[1])

    # shift wavelengths of atrans spectrum
    shft = p[0]*pixscale/oversamp
    wl_shift = lambd + shft
    atrans_new=np.interp(x=wl_shift, xp=lambd, fp=atrans_new)

    # curve atrans spectrum with Legendre polynomial
    x = np.linspace(-1, 1, len(data), dtype=np.float32)
    poly = LEGPOLY(x, p[2:])
    atrans_curved = atrans_new*poly

    diff=(data-atrans_curved)
    model=atrans_curved
    cont=poly

    return shft, model, cont, diff

def TELLSPEC_INTERP(data, atrans, pixscale, oversamp, trange):
    """
    Interpolate data and atrans onto supersampled, uniformly spaced grids.
    Args:
        data: The flux vector for the object or standard spectrum.
        atrans: The flux vector for the asmopsheric transmission spectrum.
        pixscale: Microns/pixel of the data.
        oversamp: Oversampling factor.
        trange: Max and min wavelength to interpolate over.
    Returns:
        wl_vector: Oversampled wavelength vector.
        atrans_interp: Interpolated flux vector of atmospheric transmission spectrum.
        data_interp: Interpolated flux vector of object or standard spectrum.
    """

    # wavelength range for new data
    # NOTE! This differs slightly from what was used in Newton et al.
    # (2014,2015), for which the region of interest was selected in
    # 'tell_func'
    if len(trange) != 2:
        start_wl = np.min(data[0])
        end_wl = np.max(data[0])
    else:
        start_wl = np.min(trange)
        end_wl = np.max(trange)

    # new oversampled wavelength vector on which to interpolate all data
    N_wl = (end_wl-start_wl)*oversamp/pixscale
    wl_vector = np.linspace(start_wl, end_wl, int(N_wl))

    # interpolate atrans and object flux onto wl_vector
    roi = np.argwhere((atrans[0] >= start_wl-0.1) & (atrans[0] < end_wl+0.1) & ~np.isnan(atrans[1]))[:,0]
    f = interp1d(atrans[0,roi], atrans[1,roi], kind='cubic')
    atrans_interp = f(wl_vector)

    roi = np.argwhere((data[0] >= start_wl-0.1) & (data[0] <= end_wl+0.1) & ~np.isnan(data[1]) & (data[1] >= 0))[:,0]
    f = interp1d(data[0,roi], data[1,roi], kind='cubic')
    data_interp = f(wl_vector)

    return wl_vector, atrans_interp, data_interp

def TELL_MODEL(atrans, data, plorder=5, trange=0, maxshft=0, oversamp=1, pixscale=0):
    """
    Modify the atmospheric transmission spectrum until it matches the observation to find the necessary wavelength shift.
    Args:
        data: The flux vector for the object or standard spectrum.
        atrans: The flux vector for the asmopsheric transmission spectrum.
        plorder: Order of Legendre polynomial.
        trange: Max and min wavelength to interpolate over.
        maxshft: Maximum number of microns atrans can shift.
        oversamp: Oversampling factor.
        pixscale: Microns/pixel of the data.
    Returns:
        data_new: Flattened object/standard spectrum.
        atrans_new: Shifted atmospheric transmission spectrum.
        shft: Shift in microns between the atmospheric transmission sectrum and the object/standard spectrum.
    """

    if trange == 0:
        trange = [data[0,0], data[0,-1]]
    if pixscale == 0:
        pixscale = np.mean(data[0,1:-1]-data[0,0:-2])
    if maxshft == 0:
        maxshft = 0.0008

    # Interpolate spectra onto oversampled grid.
    lambda_interp, atrans_interp, data_interp = TELLSPEC_INTERP(data, atrans, pixscale, oversamp, trange)

    # fit to find shift in wavelength
    p0 = np.concatenate([np.array([0,2]), np.ones(plorder)])
    res = least_squares(TELL_FUNC_fit,
                        x0=p0,
                        args=(lambda_interp,atrans_interp,data_interp,pixscale,oversamp,maxshft),
                        method="lm")

    # res.x[0] is shift in pixels
    # res.x[1] is atrans flux scaling
    # res.x[2:] are Legendre polynomial coefficients
    shft, model, cont, diff = TELL_FUNC(res.x, lambd=lambda_interp, atrans=atrans_interp, data=data_interp, pixscale=pixscale, oversamp=oversamp)

    data_new = [[lambda_interp+shft],[data_interp/cont]]
    atrans_new = [[lambda_interp],[atrans_interp**res.x[1]]]

    return data_new, atrans_new, shft