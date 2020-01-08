# Overview

[![Unix Build Status](https://img.shields.io/travis/omdv/geovectors/master.svg?label=unix)](https://travis-ci.org/omdv/geovectors)
[![Coverage Status](https://img.shields.io/coveralls/omdv/geovectors/master.svg)](https://coveralls.io/r/omdv/geovectors)
[![Scrutinizer Code Quality](https://img.shields.io/scrutinizer/g/omdv/geovectors.svg)](https://scrutinizer-ci.com/g/omdv/geovectors/?branch=master)
[![PyPI Version](https://img.shields.io/pypi/v/GeoVectors.svg)](https://pypi.org/project/GeoVectors)
[![PyPI License](https://img.shields.io/pypi/l/GeoVectors.svg)](https://pypi.org/project/GeoVectors)

This library provides vectorized direct and inverse geodesic methods.

The motivation was to have the accurate and fast vectorized geodesic routines for sailboat routing project. There are few python libraries, with [geographiclib](https://geographiclib.sourceforge.io/html/python/index.html) being the most accurate and reliable. Haversine method, which is widely used as an example of fast inverse method can be vectorized rather easily, however the errors are expected to be at [least 0.5%](https://en.wikipedia.org/wiki/Haversine_formula#Formulation). There are no vectorized AND accurate options.

This library is based on `numpy` and uses [Vincenty's formulae](https://en.wikipedia.org/wiki/Vincenty's_formulae). It is heavily based on the [Movable Type Scripts blog](https://www.movable-type.co.uk/scripts/latlong-vincenty.html) and Javascript [Geodesy](https://www.npmjs.com/package/geodesy) code.

Vincenty's inverse algorithm is accurate, but sensitive to [nearly antipodal points](https://en.wikipedia.org/wiki/Vincenty%27s_formulae#Nearly_antipodal_points). One approach would be to return `NaN` for such points, with the assumption that they are not frequently observed in practical applications, however as [this discussion](https://gis.stackexchange.com/questions/84885/difference-between-vincenty-and-great-circle-distance-calculations) nicely pointed out the package cannot be complete if it cannot handle these situations. I found that the issue can be solved by relaxing one of convergence criteria, but it results in errors up to 0.25% vs geographiclib for these points.

So, instead, this library uses the vectorized Vincenty's formulae with **geographiclib as a fallback** for unconverged points.

See [notebook](https://github.com/omdv/geovectors/blob/master/notebooks/demo.ipynb) for execution time comparisons vs geographiclib.

```
Direct method for 100,000 points

94.9 ms ± 25 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
vs
9.79 s ± 1.4 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

```
Inverse method for 100,000 points

1.5 s ± 504 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
vs
24.2 s ± 3.91 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

# Setup

## Requirements

* Python 3.7+
* Numpy
* Geographiclib

## Installation

Install it directly into an activated virtual environment:

```text
$ pip install GeoVectors
```

# Usage

After installation, the package can imported:

```text
$ python
>>> from geovectorslib import direct, inverse
>>> direct(lats1, lon1, bearings, distances)
>>> inverse(lats1, lons1, lats2, lons2)
```

```text
Latitudes in decimal degrees [-90; +90].
Longitudes in decimal degrees [-180; +180].
Bearings in decimal degrees [0; 360].
Distances in meters.
```

# References

[Movable Type Scripts](https://www.movable-type.co.uk/scripts/latlong-vincenty.html)

[Geodesy](https://www.npmjs.com/package/geodesy)

[Geopy](https://pypi.org/project/geopy/)

[Geographiclib](https://geographiclib.sourceforge.io/html/python/index.html)

[Stackoverflow discussion](https://gis.stackexchange.com/questions/84885/difference-between-vincenty-and-great-circle-distance-calculations)