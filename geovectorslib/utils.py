import numpy as np


def wrap90deg(deg: 'np.ndarray') -> 'np.ndarray':
    """
    Make sure that degrees (e.g. lats) is within -90;90 range.

    Ref: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """
    mask = (deg < -90) | (deg > 90)
    deg[mask] = np.abs((deg[mask] % 360 + 270) % 360 - 180) - 90
    return deg


def wrap180deg(deg: 'np.ndarray') -> 'np.ndarray':
    """
    Make sure that degrees (e.g. lon) is within +/-(0-180) range.

    Ref: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """
    mask = (deg <= -180) | (deg >= 180)
    deg[mask] = (deg[mask] + 540) % 360 - 180
    return deg


def wrap360deg(deg: 'np.ndarray') -> 'np.ndarray':
    """
    Make sure that degrees (e.g. heading) is within 360 range.

    Ref: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """
    mask = (deg < 0) | (deg > 360)
    deg[mask] = (deg[mask] % 360 + 360) % 360
    return deg
