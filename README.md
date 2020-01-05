# Overview

[![Unix Build Status](https://img.shields.io/travis/omdv/geovectors/master.svg?label=unix)](https://travis-ci.org/omdv/geovectors)
[![Windows Build Status](https://img.shields.io/appveyor/ci/omdv/geovectors/master.svg?label=window)](https://ci.appveyor.com/project/omdv/geovectors)
[![Coverage Status](https://img.shields.io/coveralls/omdv/geovectors/master.svg)](https://coveralls.io/r/omdv/geovectors)
[![Scrutinizer Code Quality](https://img.shields.io/scrutinizer/g/omdv/geovectors.svg)](https://scrutinizer-ci.com/g/omdv/geovectors/?branch=master)
[![PyPI Version](https://img.shields.io/pypi/v/GeoVectors.svg)](https://pypi.org/project/GeoVectors)
[![PyPI License](https://img.shields.io/pypi/l/GeoVectors.svg)](https://pypi.org/project/GeoVectors)

Python library for vectorized geodesic calculations with both direct and inverse methods. While working on a project for sailboat weather routing I needed a fast and accurate geodesic methods for direct and inverse great circle route calculations.

I found few python libraries with such functionality, however none of them had a vectorized form.

The two most commonly used methods are Haversine and Vincenty. Haversine is approximate, so errors may accumulate over long distances, but very computationally effective. Vincenty is a standard method, delivering a millimeter-level accuracies.

This library is based on `numpy` and uses Vincenty implementation and in particular is based on this excellent [blog](https://www.movable-type.co.uk/scripts/latlong-vincenty.html) and javascript [geodesy](https://www.npmjs.com/package/geodesy) package.

### Speed improvement
O(n) for non-vectorized vs O(log n) for this one.
Details to be added.

### Known issues
The implementation of Inverse problem is sensitive to near-antipodal points. It is a common issue, shared by geodesy package as well, try for instance the inverse problem between [16.24568372,124.84613035] and [-16.70728358, -55.2234313]. I expect my package to converge, but I found the error for such near antipodal points to reach 0.25% vs Geographicslib.

# Setup

## Requirements

* Python 3.7+

## Installation

Install it directly into an activated virtual environment:

```text
$ pip install GeoVectors
```

or add it to your [Poetry](https://poetry.eustace.io/) project:

```text
$ poetry add GeoVectors
```

# Usage

After installation, the package can imported:

```text
$ python
>>> import geovectorslib
>>> geovectorslib.__version__
```

# References
[Movable Type Scripts](https://www.movable-type.co.uk/scripts/latlong-vincenty.html)
[Geodesy](https://www.npmjs.com/package/geodesy)
[Geopy](https://pypi.org/project/geopy/)
[Geographiclib](https://geographiclib.sourceforge.io/html/python/index.html)