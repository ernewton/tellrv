import numpy as np
from scipy.interpolate import CubicSpline, interp1d
from scipy.signal import correlate
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

def find_rv0(shft, lambd, flat_obj, flat_std, pixscale, oversamp):
    """
    Function for finding the wavelength shift between the object and standard.
    Input to Levenberg-Marquardt fit in ERN_RV.
    Args:
        shft: Pixel shift between object and standard spectra.
        lambd: The wavelength vector for the spectra.
        flat_obj: Flattened object flux.
        flat_std: Flattened standard flux.
        pixscale: Microns/pixel of the data.
        oversamp: Oversampling factor.
    Returns:
        Difference in flux between shifted object and standard at each wavelength.
        To be minimized.
    """

    delta_wl = shft*pixscale/oversamp
    wl_shift = lambd + delta_wl
    obj_new=np.interp(x=wl_shift, xp=lambd, fp=flat_obj)
    diff=(obj_new-flat_std)
    return diff

def ERN_RV(data, std, pixscale=None, wrange=None, oversamp=None, shift_algo="LM", quiet=0):
    """
    Calculates the radial velocity of the object relative to the standard.
    Args:
        data: Object data.
        std: Standard data.
        pixscale: Microns/pixel of the data.
        wrange: Wavelength range to calculate shift over.
        oversamp: Oversampling factor.
        shift_algo: Algorithm to use when calculating the RV shift.
                    "LM" for Levenberg-Marquardt fit.
                    "CC" for cross-correlation.
        quiet: 0 for verbose, 1 for quiet output.
    Returns:
        Relative radial velocity shift, between the object and standard.
    """

    # remove nans from data and std
    data = data[:,~np.isnan(data[1])]
    std = std[:,~np.isnan(std[1])]

    if pixscale is None:
        pixscale = np.abs(np.median(data[0,0:-2]-data[0,1:-1]))
    if wrange is None:
        wrange = [data[0,0], data[0,-1]]
    if oversamp is None:
        oversamp = 6

    if shift_algo == "CC":
        # select wavelength range (log lambda)
        start_wl = np.log(wrange[0])
        end_wl = np.log(wrange[1])
        roi = np.argwhere((data[0] >= wrange[0]) & (data[0] <= wrange[1]))[:,0]

        # create oversampled grid uniformly spaced in log lambda
        N_wl = (end_wl-start_wl)*oversamp*np.mean(data[0,roi])/pixscale
        # wl_vector = np.interp(np.arange(N_wl), (0, N_wl-1), (start_wl, end_wl))
        wl_vector = np.linspace(start_wl, end_wl, int(N_wl))

        # interpolate object and standard onto new grid
        f_obj = interp1d(np.log(data[0]), data[1], kind='cubic')
        int_obj = f_obj(wl_vector)
        sroi = np.argwhere((std[0] >= wrange[0]) & (std[0] <= wrange[1]))[:,0]
        f_std = interp1d(np.log(std[0]), std[1], kind='cubic')
        int_std = f_std(wl_vector)

        # custom flattening routine
        idxs = np.arange(len(int_obj))
        bin_edges = np.linspace(0, len(int_obj), 10)
        binned_idxs = np.array([np.mean([bin_edges[i], bin_edges[i+1]]) for i in range(len(bin_edges)-1)])
        bin_idxs = np.array([np.argwhere((idxs > bin_edges[i]) & (idxs < bin_edges[i+1]))[:,0] for i in range(len(bin_edges)-1)])
        binned_idxs = np.concatenate([[0], binned_idxs, [len(int_obj)]])
        binned_obj = np.array([np.mean(int_obj[idx]) for idx in bin_idxs])
        binned_obj = np.concatenate([[int_obj[0]], binned_obj, [int_obj[-1]]])
        binned_std = np.array([np.mean(int_std[idx]) for idx in bin_idxs])
        binned_std = np.concatenate([[int_std[0]], binned_std, [int_std[-1]]])

        cubic_obj = interp1d(binned_idxs, binned_obj, kind='cubic')
        interpolated_obj = cubic_obj(idxs)
        cubic_std = interp1d(binned_idxs, binned_std, kind='cubic')
        interpolated_std = cubic_std(idxs)

        flat_obj = int_obj/interpolated_obj - 1
        flat_std = int_std/interpolated_std - 1

        # cross-correlate to find shift
        padded_obj = np.concatenate([np.zeros(len(flat_obj)//2),flat_obj,np.zeros(len(flat_obj)//2)])
        corr = np.correlate(padded_obj, flat_std, "same")
        if len(flat_obj)%2 != 0:
            lag = np.linspace(-len(flat_obj)+1, len(flat_obj)-1, len(corr))
        else:
            lag = np.linspace(-len(flat_obj), len(flat_obj)-1, len(corr))
        p = np.argmax(corr)
        shft = -lag[p]
        # quadratic offset
        aa = corr[p]
        bb = 0.5*(corr[p+1] - corr[p-1])
        cc = 0.5*(corr[p+1] + corr[p-1] - 2.0*aa)
        quad_offset = -0.5*bb/cc
        shft = -lag[p] - quad_offset

    elif shift_algo == "LM":
        # select wavelength range
        start_wl = wrange[0]
        end_wl = wrange[1]
        roi = np.argwhere((data[0] >= wrange[0]) & (data[0] <= wrange[1]))[:,0]

        # create oversampled grid uniformly spaced in lambda
        N_wl = (end_wl-start_wl)*oversamp*np.mean(data[0,roi])/pixscale
        wl_vector = np.linspace(start_wl, end_wl, int(N_wl))

        # interpolate object and standard onto new grid
        f_obj = interp1d(data[0], data[1], kind='cubic')
        int_obj = f_obj(wl_vector)
        sroi = np.argwhere((std[0] >= wrange[0]) & (std[0] <= wrange[1]))[:,0]
        f_std = interp1d(std[0], std[1], kind='cubic')
        int_std = f_std(wl_vector)

        # custom flattening routine
        idxs = np.arange(len(int_obj))
        bin_edges = np.linspace(0, len(int_obj), 10)
        binned_idxs = np.array([np.mean([bin_edges[i], bin_edges[i+1]]) for i in range(len(bin_edges)-1)])
        bin_idxs = np.array([np.argwhere((idxs > bin_edges[i]) & (idxs < bin_edges[i+1]))[:,0] for i in range(len(bin_edges)-1)])
        binned_idxs = np.concatenate([[0], binned_idxs, [len(int_obj)]])
        binned_obj = np.array([np.mean(int_obj[idx]) for idx in bin_idxs])
        binned_obj = np.concatenate([[int_obj[0]], binned_obj, [int_obj[-1]]])
        binned_std = np.array([np.mean(int_std[idx]) for idx in bin_idxs])
        binned_std = np.concatenate([[int_std[0]], binned_std, [int_std[-1]]])

        cubic_obj = interp1d(binned_idxs, binned_obj, kind='cubic')
        interpolated_obj = cubic_obj(idxs)
        cubic_std = interp1d(binned_idxs, binned_std, kind='cubic')
        interpolated_std = cubic_std(idxs)

        flat_obj = int_obj/interpolated_obj - 1
        flat_std = int_std/interpolated_std - 1

        # use Levenberg-Marquartd fit to find shift
        res = least_squares(find_rv0,
                            x0=[0],
                            args=(wl_vector,flat_obj,flat_std,pixscale,oversamp),
                            method="lm")
        shft = -res.x[0]

    # offset and relative RV
    # if shift is positive, shifting object to redder wavelengths = it was blue shifted = negative RV
    offset = -shft*pixscale/(oversamp*np.mean(data[0,roi])) # delta lambda/lambda
    rv0 = 3.0e5*offset                                      # relative obj-to-std radial velocity

    return rv0