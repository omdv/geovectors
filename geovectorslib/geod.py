import numpy as np

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
    tanU1 = (1 - f) * np.tan(lat1)
    cosU1 = 1 / np.sqrt((1 + tanU1 * tanU1))
    sinU1 = tanU1 * cosU1
    tanU2 = (1 - f) * np.tan(lat2)
    cosU2 = 1 / np.sqrt((1 + tanU2 * tanU2))
    sinU2 = tanU2 * cosU2

    # checks for antipodal points
    antipodal = (np.abs(L) > np.pi / 2) | (np.abs(lat2 - lat1) > np.pi / 2)

    # delta = difference in longitude on an auxiliary sphere
    delta = L

    # sigma = angular distance P₁ P₂ on the sphere
    sigma = np.zeros(lat1.shape)
    sigma[antipodal] = np.pi

    cossigma = np.ones(lat1.shape)
    cossigma[antipodal] = -1

    # sigmam = angular distance on the sphere from the equator
    # to the midpoint of the line
    # azi = azimuth of the geodesic at the equator
    deltaprime = np.zeros(lat1.shape)
    iterations = 0

    # init before loop to allow mask indexing
    sinSqsigma = np.zeros(lat1.shape)
    sinsigma = np.zeros(lat1.shape)
    cossigma = np.zeros(lat1.shape)
    sinazi = np.zeros(lat1.shape)
    cosSqazi = np.ones(lat1.shape)
    cos2sigmam = np.ones(lat1.shape)
    C = np.ones(lat1.shape)

    # init mask
    m = np.ones(lat1.shape, dtype=bool)

    while (np.abs(delta[m] - deltaprime[m]) > 1e-12).any():
        sindelta = np.sin(delta)
        cosdelta = np.cos(delta)
        sinSqsigma = (cosU2 * sindelta) * (cosU2 * sindelta) + (
            cosU1 * sinU2 - sinU1 * cosU2 * cosdelta
        ) * (cosU1 * sinU2 - sinU1 * cosU2 * cosdelta)

        # co-incident/antipodal points mask - exclude from the rest of the loop
        m = (np.abs(sinSqsigma) > eps) & (sinSqsigma != np.nan)

        sinsigma[m] = np.sqrt(sinSqsigma[m])
        cossigma[m] = sinU1[m] * sinU2[m] + cosU1[m] * cosU2[m] * cosdelta[m]
        sigma[m] = np.arctan2(sinsigma[m], cossigma[m])
        sinazi[m] = cosU1[m] * cosU2[m] * sindelta[m] / sinsigma[m]
        cosSqazi[m] = 1 - sinazi[m] * sinazi[m]

        # on equatorial line cos²azi = 0
        cos2sigmam[m] = cossigma[m] - 2 * sinU1[m] * sinU2[m] / cosSqazi[m]
        cos2sigmam[cosSqazi == 0] = 0

        C[m] = f / 16 * cosSqazi[m] * (4 + f * (4 - 3 * cosSqazi[m]))

        deltaprime = delta
        delta = L + (1 - C) * f * sinazi * (
            sigma
            + C
            * sinsigma
            * (cos2sigmam + C * cossigma * (-1 + 2 * cos2sigmam * cos2sigmam))
        )

        # Exceptions
        iterationCheck = np.abs(delta)
        iterationCheck[antipodal] = iterationCheck[antipodal] - np.pi
        if (iterationCheck > np.pi).any():
            raise Exception('delta > np.pi')

        iterations += 1
        if iterations >= 20:
            raise Exception('Vincenty formula failed to converge')

    uSq = cosSqazi * (a * a - b * b) / (b * b)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    dsigma = (
        B
        * sinsigma
        * (
            cos2sigmam
            + B
            / 4
            * (
                cossigma * (-1 + 2 * cos2sigmam * cos2sigmam)
                - B
                / 6
                * cos2sigmam
                * (-3 + 4 * sinsigma * sinsigma)
                * (-3 + 4 * cos2sigmam * cos2sigmam)
            )
        )
    )
    # length of the geodesic
    s12 = b * A * (sigma - dsigma)

    # note special handling of exactly antipodal points where sin²sigma = 0
    # (due to discontinuity atan2(0, 0) = 0 but atan2(eps, 0) = np.pi/2 / 90°)
    # in which case bearing is always meridional, due north (or due south!)

    azi1 = np.arctan2(cosU2 * sindelta, cosU1 * sinU2 - sinU1 * cosU2 * cosdelta)
    azi1[np.abs(sinSqsigma) < eps] = 0

    azi2 = np.arctan2(cosU1 * sindelta, -sinU1 * cosU2 + cosU1 * sinU2 * cosdelta)
    azi2[np.abs(sinSqsigma) < eps] = np.pi

    azi1 = wrap360deg(azi1 * 180 / np.pi)
    azi2 = wrap360deg(azi2 * 180 / np.pi)

    # if distance is too small
    mask = s12 < eps
    azi1[mask] = None
    azi2[mask] = None

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

    sinazi1 = np.sin(azi1)
    cosazi1 = np.cos(azi1)

    tanU1 = (1 - f) * np.tan(lat1)
    cosU1 = 1 / np.sqrt((1 + tanU1 * tanU1))
    sinU1 = tanU1 * cosU1

    # sigma1 = angular distance on the sphere from the equator to P1
    sigma1 = np.arctan2(tanU1, cosazi1)

    sinazi = cosU1 * sinazi1  # azi = azimuth of the geodesic at the equator
    cosSqazi = 1 - sinazi * sinazi
    uSq = cosSqazi * (a * a - b * b) / (b * b)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))

    sigma = s / (b * A)
    # sigma = angular distance P₁ P₂ on the sphere
    # sigmam = angular distance on the sphere from the equator to
    # the midpoint of the line

    iterations = 0
    sigmaprime = 0

    while (np.abs(sigma - sigmaprime) > 1e-12).any():
        cos2sigmam = np.cos(2 * sigma1 + sigma)
        sinsigma = np.sin(sigma)
        cossigma = np.cos(sigma)
        deltasigma = (
            B
            * sinsigma
            * (
                cos2sigmam
                + B
                / 4
                * (
                    cossigma * (-1 + 2 * cos2sigmam * cos2sigmam)
                    - B
                    / 6
                    * cos2sigmam
                    * (-3 + 4 * sinsigma * sinsigma)
                    * (-3 + 4 * cos2sigmam * cos2sigmam)
                )
            )
        )
        sigmaprime = sigma
        sigma = s / (b * A) + deltasigma
        iterations += 1
        if iterations >= 100:
            raise Exception('Vincenty formula failed to converge')

    x = sinU1 * sinsigma - cosU1 * cossigma * cosazi1
    lat2 = np.arctan2(
        sinU1 * cossigma + cosU1 * sinsigma * cosazi1,
        (1 - f) * np.sqrt(sinazi * sinazi + x * x),
    )
    lambd = np.arctan2(
        sinsigma * sinazi1, cosU1 * cossigma - sinU1 * sinsigma * cosazi1
    )
    C = f / 16 * cosSqazi * (4 + f * (4 - 3 * cosSqazi))
    L = lambd - (1 - C) * f * sinazi * (
        sigma
        + C
        * sinsigma
        * (cos2sigmam + C * cossigma * (-1 + 2 * cos2sigmam * cos2sigmam))
    )
    lon2 = lon1 + L
    azi2 = np.arctan2(sinazi, -x)

    lat2 = wrap90deg(lat2 * 180.0 / np.pi)
    lon2 = wrap180deg(lon2 * 180.0 / np.pi)
    azi2 = wrap360deg(azi2 * 180.0 / np.pi)

    return {'lat2': lat2, 'lon2': lon2, 'azi2': azi2, 'iterations': iterations}
