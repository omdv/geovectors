from collections import namedtuple

def ellipsoid_constants():
    """
    Define available ellipsoid geometries.
    Ref: geopy/internet
    """
    e = namedtuple('Ellipsoid', ['major', 'minor', 'flat'])
    ellipsoids = {
        'ETRS89': e(6378137.0, 6356752.314140, 1. / 298.257222101),
        'WGS84': e(6378137.0, 6356752.3142, 1. / 298.257223563),
        'GRS80': e(6378137.0, 6356752.3141, 1. / 298.257222101),
        'GRS67': e(6378160.0, 6356774.719, 1. / 298.25),
        'Intl-1924': e(6378388.0, 6356911.946, 1. / 297.0), 
        'Clarke-1880': e(6378249.145, 6356514.86955, 1. / 293.465),
        'Airy-1830': e(6377563.396, 6356256.909, 1. / 299.3249646),
    }
    return ellipsoids
