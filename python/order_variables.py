from spex_order_variables import spex_order_variables

def order_variables(hdr, order, instrument=None):
    """
    Finds the variables for the Legendre polynomial in TELL_FUNC, which depends on the iunstrument and spectral order.
    Args:
        hdr: FITS header for the data.
        order: Spectral order (usually 0 is K band, 1 is H band, 2 is J band).
        instrument: What instrument was used to collect the observations (only spex compatible at the moment.)
    Returns:
        Variables for the Legendre polynomial fit.
    """

    if instrument is not None:
        if instrument == 'spex':
            return spex_order_variables(hdr, order)
        # elif instrument == 'fire':
        #     return fire_order_variables(hdr, order, wrange, trange, pixscale, polydegree)
        else:
            print('Invalid telescope')
    else:
        return spex_order_variables(hdr, order)

