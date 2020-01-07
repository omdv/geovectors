import numpy as np
from geographiclib.geodesic import Geodesic

from geovectorslib.utils import wrap90deg, wrap180deg, wrap360deg


def inverse(lats1: 'list', lons1: 'list', lats2: 'list', lons2: 'list') -> 'dict':
    """
    Inverse geodesic using Vincenty approach.

    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.

    Ref:
    https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """

    eps = np.finfo(float).eps

    lon1, lon2 = map(wrap180deg, [np.array(lons1), np.array(lons2)])
    lat1, lat2 = map(wrap90deg, [np.array(lats1), np.array(lats2)])
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # WGS84
    a = 6378137.0
    b = 6356752.314245
    f = 1.0 / 298.257223563

    # L = difference in longitude
    L = lon2 - lon1

    # U = reduced latitude, defined by tan U = (1-f)·tanφ.
    tan_U1 = (1 - f) * np.tan(lat1)
    cos_U1 = 1 / np.sqrt((1 + tan_U1 * tan_U1))
    sin_U1 = tan_U1 * cos_U1
    tan_U2 = (1 - f) * np.tan(lat2)
    cos_U2 = 1 / np.sqrt((1 + tan_U2 * tan_U2))
    sin_U2 = tan_U2 * cos_U2

    # checks for antipodal points
    antipodal = (np.abs(L) > np.pi / 2) | (np.abs(lat2 - lat1) > np.pi / 2)

    # delta = difference in longitude on an auxiliary sphere
    delta = L
    sin_delta = np.sin(delta)
    cos_delta = np.cos(delta)

    # sigma = angular distance P₁ P₂ on the sphere
    sigma = np.zeros(lat1.shape)
    sigma[antipodal] = np.pi

    cos_sigma = np.ones(lat1.shape)
    cos_sigma[antipodal] = -1

    # sigmam = angular distance on the sphere from the equator
    # to the midpoint of the line
    # azi = azimuth of the geodesic at the equator
    delta_prime = np.zeros(lat1.shape)
    iterations = 0

    # init before loop to allow mask indexing
    sin_Sq_sigma = np.zeros(lat1.shape)
    sin_sigma = np.zeros(lat1.shape)
    cos_sigma = np.zeros(lat1.shape)
    sin_azi = np.zeros(lat1.shape)
    cos_Sq_azi = np.ones(lat1.shape)
    cos_2_sigma_m = np.ones(lat1.shape)
    C = np.ones(lat1.shape)

    # init mask
    m = np.ones(lat1.shape, dtype=bool)
    conv_mask = np.zeros(lat1.shape, dtype=bool)

    while (np.abs(delta[m] - delta_prime[m]) > 1e-12).any():
        sin_delta = np.sin(delta)
        cos_delta = np.cos(delta)
        sin_Sq_sigma = (cos_U2 * sin_delta) * (cos_U2 * sin_delta) + (
            cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_delta
        ) * (cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_delta)

        # co-incident/antipodal points mask - exclude from the rest of the loop
        # the value has to be about 1.e-4 experimentally
        # otherwise there are issues with convergence for near antipodal points
        m = (np.abs(sin_Sq_sigma) > eps) & (sin_Sq_sigma != np.nan)

        sin_sigma[m] = np.sqrt(sin_Sq_sigma[m])
        cos_sigma[m] = sin_U1[m] * sin_U2[m] + cos_U1[m] * cos_U2[m] * cos_delta[m]
        sigma[m] = np.arctan2(sin_sigma[m], cos_sigma[m])
        sin_azi[m] = cos_U1[m] * cos_U2[m] * sin_delta[m] / sin_sigma[m]
        cos_Sq_azi[m] = 1 - sin_azi[m] * sin_azi[m]

        # on equatorial line cos²azi = 0
        cos_2_sigma_m[m] = cos_sigma[m] - 2 * sin_U1[m] * sin_U2[m] / cos_Sq_azi[m]
        cos_2_sigma_m[cos_Sq_azi == 0] = 0

        C[m] = f / 16 * cos_Sq_azi[m] * (4 + f * (4 - 3 * cos_Sq_azi[m]))

        delta_prime = delta
        delta = L + (1 - C) * f * sin_azi * (
            sigma
            + (C * sin_sigma)
            * (cos_2_sigma_m + C * cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m))
        )

        # Exceptions
        iteration_check = np.abs(delta)
        iteration_check[antipodal] = iteration_check[antipodal] - np.pi
        if (iteration_check > np.pi).any():
            raise Exception('delta > np.pi')

        iterations += 1
        if iterations >= 50:
            conv_mask = np.abs(delta - delta_prime) > 1e-12
            break

    uSq = cos_Sq_azi * (a * a - b * b) / (b * b)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    d_sigma = (B * sin_sigma) * (
        cos_2_sigma_m
        + (B / 4)
        * (
            cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m)
            - (B / 6 * cos_2_sigma_m)
            * (-3 + 4 * sin_sigma * sin_sigma)
            * (-3 + 4 * cos_2_sigma_m * cos_2_sigma_m)
        )
    )
    # length of the geodesic
    s12 = b * A * (sigma - d_sigma)

    # note special handling of exactly antipodal points where sin²sigma = 0
    # (due to discontinuity atan2(0, 0) = 0 but atan2(eps, 0) = np.pi/2 / 90°)
    # in which case bearing is always meridional, due north (or due south!)

    azi1 = np.arctan2(cos_U2 * sin_delta, cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_delta)
    azi1[np.abs(sin_Sq_sigma) < eps] = 0

    azi2 = np.arctan2(
        cos_U1 * sin_delta, -sin_U1 * cos_U2 + cos_U1 * sin_U2 * cos_delta
    )
    azi2[np.abs(sin_Sq_sigma) < eps] = np.pi

    azi1 = wrap360deg(azi1 * 180 / np.pi)
    azi2 = wrap360deg(azi2 * 180 / np.pi)

    # if distance is too small
    mask = s12 < eps
    azi1[mask] = None
    azi2[mask] = None

    # use geographiclib for points which didn't converge
    if conv_mask[conv_mask].shape[0] > 0:
        vInverse = np.vectorize(Geodesic.WGS84.Inverse)
        _ps = vInverse(
            lats1[conv_mask], lons1[conv_mask], lats2[conv_mask], lons2[conv_mask]
        )
        s12[conv_mask] = [_p['s12'] for _p in _ps]
        azi1[conv_mask] = [_p['azi1'] for _p in _ps]
        azi2[conv_mask] = [_p['azi2'] for _p in _ps]

    return {'s12': s12, 'azi1': azi1, 'azi2': azi2, 'iterations': iterations}


def direct(lats1: 'list', lons1: 'list', brgs: 'list', dists: 'list') -> 'dict':
    """
    Direct geodesic using Vincenty approach.

    Given a start point, initial bearing (0° is North, clockwise)
    and distance in meters, this will calculate the destination point
    and final bearing travelling along a (shortest distance) great circle arc.

    Bearing in 0-360deg starting from North clockwise.

    all variables must be of same length

    Ref:
    https://www.movable-type.co.uk/scripts/latlong.html
    https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """

    lon1, lat1, brg = map(np.radians, [lons1, lats1, brgs])

    # WGS84
    a = 6378137.0
    b = 6356752.314245
    f = 1.0 / 298.257223563

    azi1 = brg
    s = dists

    sin_azi1 = np.sin(azi1)
    cos_azi1 = np.cos(azi1)

    tan_U1 = (1 - f) * np.tan(lat1)
    cos_U1 = 1 / np.sqrt((1 + tan_U1 * tan_U1))
    sin_U1 = tan_U1 * cos_U1

    # sigma1 = angular distance on the sphere from the equator to P1
    sigma1 = np.arctan2(tan_U1, cos_azi1)

    sin_azi = cos_U1 * sin_azi1  # azi = azimuth of the geodesic at the equator
    cos_Sq_azi = 1 - sin_azi * sin_azi
    uSq = cos_Sq_azi * (a * a - b * b) / (b * b)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))

    sigma = s / (b * A)
    # sigma = angular distance P₁ P₂ on the sphere
    # sigmam = angular distance on the sphere from the equator to
    # the midpoint of the line
    sin_sigma = np.sin(sigma)
    cos_sigma = np.cos(sigma)
    cos_2_sigma_m = np.cos(2 * sigma1 + sigma)

    iterations = 0
    sigma_prime = 0

    while (np.abs(sigma - sigma_prime) > 1e-12).any():
        cos_2_sigma_m = np.cos(2 * sigma1 + sigma)
        sin_sigma = np.sin(sigma)
        cos_sigma = np.cos(sigma)
        delta_sigma = (B * sin_sigma) * (
            cos_2_sigma_m
            + B
            / 4
            * (
                cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m)
                - (B / 6 * cos_2_sigma_m)
                * (-3 + 4 * sin_sigma * sin_sigma)
                * (-3 + 4 * cos_2_sigma_m * cos_2_sigma_m)
            )
        )
        sigma_prime = sigma
        sigma = s / (b * A) + delta_sigma
        iterations += 1
        if iterations >= 50:
            raise Exception('Vincenty formula failed to converge')

    x = sin_U1 * sin_sigma - cos_U1 * cos_sigma * cos_azi1
    lat2 = np.arctan2(
        sin_U1 * cos_sigma + cos_U1 * sin_sigma * cos_azi1,
        (1 - f) * np.sqrt(sin_azi * sin_azi + x * x),
    )
    lambda_ = np.arctan2(
        sin_sigma * sin_azi1, cos_U1 * cos_sigma - sin_U1 * sin_sigma * cos_azi1
    )
    C = f / 16 * cos_Sq_azi * (4 + f * (4 - 3 * cos_Sq_azi))
    L = lambda_ - (1 - C) * f * sin_azi * (
        sigma
        + C
        * sin_sigma
        * (cos_2_sigma_m + C * cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m))
    )
    lon2 = lon1 + L
    azi2 = np.arctan2(sin_azi, -x)

    lat2 = wrap90deg(lat2 * 180.0 / np.pi)
    lon2 = wrap180deg(lon2 * 180.0 / np.pi)
    azi2 = wrap360deg(azi2 * 180.0 / np.pi)

    return {'lat2': lat2, 'lon2': lon2, 'azi2': azi2, 'iterations': iterations}
