#!/usr/bin/env python

try:
    from setuptools import setup, Extension
except:
    from distutils.core import setup, Extension

setup(name='SkyLinesPolyEncode',
    version='0.1.1',
    description="SkyLines Polyline encoding (C++ extension)",
    long_description="Multidimensional encoding of line & polygon coordinates for use in SkyLines. Based on cGPolyEncode by Robert Coup <robert.coup@koordinates.com>",
    author='Tobias Lohner',
    author_email='tobias@lohner-net.de',
    provides=['skylinespolyencode'],
    keywords='gis,geospatial,google-maps,gmaps,mapping',
    url='http://github.com/TobiasLohner/SkyLinesPolyEncode/',
    ext_modules=[
        Extension("skylinespolyencode", ["skylinespolyencode_py.cpp", "SkyLinesPolyEncoder.cpp"]),
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Web Environment',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Utilities'],
)

