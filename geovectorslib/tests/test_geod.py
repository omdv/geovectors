# pylint: disable=redefined-outer-name,unused-variable,expression-not-assigned,singleton-comparison
import numpy as np
from geographiclib.geodesic import Geodesic

from geovectorslib import geod


def describe_inverse():
    def fixed(expect):
        r = geod.inverse(
            [-90, -89, 0, 89, 90, 0, 0],
            [-180, -170, 0, 150, 180, 0, -180],
            [90, 89, 0, -89, -90, 0, 0],
            [180, 170, 0, -150, -180, -180, 180],
        )
        t = {
            's12': [
                20003931.4586233,
                19783938.13763718,
                0.0,
                19810477.03907933,
                20003931.4586233,
                20003931.4586233,
                0.0,
            ],
            'azi1': [0.0, 349.99759033, np.nan, 149.99358041, 0.0, 0.0, np.nan],
            'azi2': [180.0, 349.99759033, np.nan, 149.99358041, 180.0, 180.0, np.nan],
            'iterations': 3,
        }
        np.testing.assert_allclose(r['s12'], t['s12'], rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(r['azi1'], t['azi1'], rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(r['azi2'], t['azi2'], rtol=1e-9, atol=1e-9)
        expect(r['iterations']) == t['iterations']

    def vs_geographiclib():
        glib = Geodesic.WGS84

        lats = np.arange(-90, 90, 45)
        lons = np.arange(-180, 180, 30)
        vp1 = [(p1, p2) for p1 in lats for p2 in lons]
        vp2 = [(1, 180) for _ in range(len(vp1))]

        g1s = [glib.Inverse(*vp1[i], *vp2[i]) for i in range(len(vp1))]
        g2s = geod.inverse(*list(zip(*vp1)), *list(zip(*vp2)))

        s12 = [g1s[i]['s12'] for i in range(len(g1s))]

        # not testing azimuth, as it varies from 0-360deg.
        np.testing.assert_allclose(g2s['s12'], s12, rtol=1e-5, atol=1e-5)

    def large_size():
        np.random.seed(42)
        lats1 = np.random.uniform(-90, 90, 100000)
        lons1 = np.random.uniform(-180, 180, 100000)
        lats2 = np.random.uniform(-90, 90, 100000)
        lons2 = np.random.uniform(-180, 180, 100000)
        _ = geod.inverse(lats1, lons1, lats2, lons2)

    def near_antipodal():
        glib = Geodesic.WGS84

        lats1 = np.array([0, 0, 16.24568372, 0])
        lons1 = np.array([0, -180, 124.84613035, 0])
        lats2 = np.array([0, 0, -16.70728358, 0.5])
        lons2 = np.array([180, 180, -55.2234313, 179.7])

        vInverse = np.vectorize(glib.Inverse)

        g1s = vInverse(lats1, lons1, lats2, lons2)
        g2s = geod.inverse(lats1, lons1, lats2, lons2)

        s12 = [g1s[i]['s12'] for i in range(len(g1s))]
        np.testing.assert_allclose(g2s['s12'], s12, rtol=1e-10, atol=1e-10)


def describe_direct():
    def calculate(expect):
        r = geod.direct(
            [-90, -89, 0, 89, 90],
            [-180, -170, 0, 150, 180],
            [40, 50, 300, 360, 30],
            [5000.0, 3.0e6, 6.0e6, 10.0e6, 20.0e6],
        )
        t = {
            'lat2': [-89.95523483, -62.46780606, 23.95378158, 1.02790157, -89.96480152],
            'lon2': [-140.0, -121.47529512, -49.91980895, -30.0, -30.0],
            'azi2': [
                2.88637944e-12,
                1.65855758e00,
                2.88716163e02,
                1.80000000e02,
                1.80000000e02,
            ],
            'iterations': 4,
        }
        np.testing.assert_allclose(r['lat2'], t['lat2'], rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(r['lon2'], t['lon2'], rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(r['azi2'], t['azi2'], rtol=1e-8, atol=1e-8)

        expect(r['iterations']) == t['iterations']

    def vs_geographiclib():
        glib = Geodesic.WGS84

        npoints = 20
        lats = np.arange(-90, 90, 180 / npoints)
        lons = np.arange(-180, 180, 360 / npoints)
        dist = np.arange(0, 6000.0e3, 6000.0e3 / npoints)
        brgs = np.arange(-20, 360, 380 / npoints)

        g1s = [glib.Direct(lats[i], lons[i], brgs[i], dist[i]) for i in range(npoints)]
        g2s = geod.direct(lats, lons, brgs, dist)

        lats2 = np.array([g1s[i]['lat2'] for i in range(len(g1s))])

        # not testing lons as the -180;180 range is not consistent
        np.testing.assert_allclose(g2s['lat2'], lats2)

    def large_size():
        np.random.seed(42)
        lats1 = np.random.uniform(-90, 90, 100000)
        lons1 = np.random.uniform(-180, 180, 100000)
        brgs = np.random.uniform(0, 360, 100000)
        dist = np.random.uniform(0, 20e6, 100000)
        _ = geod.inverse(lats1, lons1, brgs, dist)
