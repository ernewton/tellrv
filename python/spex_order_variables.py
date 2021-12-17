def spex_order_variables(hdr, order):
    """
    Finds variables for the Legendre polynomial in TELL_FUNC for spex observations, depending on spectral order.
    Args:
        hdr: FITS header for the data.
        order: Spectral order (usually 0 is K band, 1 is H band, 2 is J band).
    Returns:
        wrange: Wavelength range for ERN_RV.
        trange: Wavelength range for TELL_MODEL.
        pixscale: Microns/pixel for the data.
        polydegree: Order of the Legendre polynomial.
    """

    if order == 0:
        wrange = [2.18, 2.41]
        trange = [1.995, 2.4]
        pixscale = hdr['DISPO03']
        polydegree = 5
    elif order == 1:
        wrange = [1.49, 1.73]
        trange = [1.43, 1.8]
        pixscale = hdr['DISPO04']
        polydegree = 4
    elif order == 2:
        wrange = [1.15,1.32]
        trange = [1.142,1.35]
        pixscale = hdr['DISPO05']
        polydegree = 5
    elif order == 3:
        wrange = [1.0, 1.1]
        trange = [0.94, 1.17]
        pixscale = hdr['DISPO06']
        polydegree = 5
    elif order == 4:
        wrange = [0.82, 0.92]
        trange = [0.89, 0.99]
        pixscale = hdr['DISPO07']
        polydegree = 4
    elif order == 5:
        wrange = [0.81, 0.9]
        trange = [0.81, 0.9]
        pixscale = hdr[hdr, 'DISPO08']
        polydegree = 5
    else: print("Order is", order, "This is an invalid order.")

    return wrange, trange, pixscale, polydegree