#!/usr/bin/env python

import numpy

sourceFileList = ['PyWindProfile/PyWindProfile.cpp',\
                  '../Gaussian.cpp',\
                  '../Utilities.cpp',\
                  '../OpticalDepth.cpp',\
                  '../Porosity.cpp',\
                  '../AnalyticOpticalDepth.cpp',\
                  '../SmoothA1.cpp',\
                  '../IsotropicSeries.cpp',\
                  '../Series.cpp',\
                  '../NumericalOpticalDepth.cpp',\
                  '../NumericalOpticalDepthU.cpp',\
                  '../NumericalOpticalDepthZ.cpp',\
                  '../RAD_OpticalDepth.cpp',\
                  '../mal_Integration.cpp',\
                  '../Lx.cpp',\
                  '../ResonanceScattering.cpp',\
                  '../HeLikeRatio.cpp',\
                  '../mal_RootFinderNewton.cpp',\
                  '../UxRoot.cpp']

libraryDirList = ['/opt/local/lib/']
libraryNameList = ['gsl','gslcblas']
includeDirList = [numpy.get_include (), '/opt/local/include/', 'PyWindProfile']


def configuration (parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration ('PyWindProfile', parent_package, top_path)
    config.add_extension ('PyWindProfile', 
                          sources = sourceFileList, 
                          library_dirs = libraryDirList,
                          libraries = libraryNameList,
                          include_dirs = includeDirList)
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration (top_path='').todict()
    description_string\
        = 'Provides functions based on XSPEC local model package windprofile'
    setup (author='Maurice A Leutenegger',
           author_email='maurice.a.leutenegger@nasa.gov',
           description=description_string,
           version = '0.1',
           #name='PyWindProfile',
           #packages=[''],
           #package_dir={'': '.'},
           #package_data={'': ['PyWindProfile']},
           packages=['PyWindProfile'],
           **config)
