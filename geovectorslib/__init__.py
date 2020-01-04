from pkg_resources import DistributionNotFound, get_distribution

from .geod import direct, inverse


try:
    __version__ = get_distribution('GeoVectors').version
except DistributionNotFound:
    __version__ = '(local)'

__all__ = ['direct', 'inverse']
