import numpy as np
from astropy.io import fits
import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
import datetime

def GET_HELIO(hdr, degrees=0, barycentric=0, quiet=0):
    """
    Calculated heliocentric/barycentric correction based on the time and location at whch the observation was taken.
    Args:
        hdr: FITS header of the data.
        degrees: Whether or not RA and Dec are in degrees.
        barycentric: Whether or not to calculate the barycentric correction (opposed to heloicentric).
        quiet: 0 for verbose, 1 for quiet output.
    Returns:
        Heliocentric/Barycentric correction.
    """

    # step 1: define telescope longitude and latitude
    tel = hdr['TELESCOP']
    if tel == 'NASA IRTF':
        lat = 19. + 49./60. + 34.38594/3600. # north
        lon = 360. - (155. + 28./60. + 19.19564/3600.) # west
        alt = 13674.76*0.3048 # feet to meters
        date = hdr['DATE_OBS']
        time = hdr['TIME_OBS']
        degrees = 0

    # step 2: get time of observation
    date_ex = date.split('-')
    year, month, date = int(date_ex[0]), int(date_ex[1]), int(date_ex[2])
    time_ex = time.split(':')
    hour, minute, second = int(time_ex[0]), int(time_ex[1]), round(float(time_ex[2]))
    if second == 60:
        second = 0
        minute += 1
    if minute == 60:
        minute = 0
        hour += 1
    dt = datetime.datetime(year, month, date, hour, minute, second)
    obstime = astropy.time.Time(dt)

    # step 3: get ra and dec in degrees if they aren't already
    if 'TCS_RA' in list(hdr):
        ra = hdr['TCS_RA']
        dec = hdr['TCS_DEC']
    elif 'RA' in list(hdr):
        ra = hdr['RA']
        dec = hdr['DEC']
    else:
        print("Unable to find target coordinates in FITS header.")
    if degrees == 0:
        ra_ex = ra.split(':')
        ra_deg = 15 * (int(ra_ex[0]) + int(ra_ex[1])/60 + float(ra_ex[2])/3600)
        dec_ex = dec.split(':')
        dec_deg = abs( int(dec_ex[0]) + int(dec_ex[1])/60 + float(dec_ex[2])/3600)
        #if the first character of the dec is negative, then make negative
        if dec[0] == '-':
            dec_deg = -dec_deg
    elif degrees != 0:
        ra_deg = ra
        dec_deg = dec

    # step 4: caculate barycentric or heliocentric correction, depending on keyword
    # see here: https://docs.astropy.org/en/stable/coordinates/velocities.html
    loc = EarthLocation.from_geodetic(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)
    sc = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
    if barycentric == 0:
        heliocorr = sc.radial_velocity_correction('heliocentric', obstime=obstime, location=loc)  
        corr = heliocorr.to(u.km/u.s)  
        if quiet == 0:
            print("GET_HELIO: returning heliocentric.")
    else:
        barycorr = sc.radial_velocity_correction('barycentric', obstime=obstime, location=loc)  
        corr = barycorr.to(u.km/u.s)  
        if quiet == 0:
            print("GET_HELIO: returning barycentric.")

    return corr

