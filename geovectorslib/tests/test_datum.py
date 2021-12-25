# pylint: disable=redefined-outer-name,unused-variable,expression-not-assigned,singleton-comparison
from geovectorslib import datum

def test_ellipsoid_constants(expect):
    g = datum.ellipsoid_constants()
    expect(g['ETRS89'].minor) == 6356752.314140
    expect(g['WGS84'].major) == 6378137.
    expect(g['GRS80'].minor) == 6356752.3141
    expect(g['Airy-1830'].flat) == 1. / 299.3249646
    expect(g['Intl-1924'].major) == 6378388.
    expect(g['Clarke-1880'].major) == 6378249.145
    expect(g['GRS67'].minor) == 6356774.719
