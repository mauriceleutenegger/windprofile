#!/usr/bin/env python

import numpy

def configuration (parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration ('PyWindProfile', parent_package, top_path)
    sourceFiles\
        = ['PyWindProfile.cpp',\
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
               '../mal_Integration.cpp',\
               '../Lx.cpp',\
               '../ResonanceScattering.cpp',\
               '../HeLikeRatio.cpp',\
               '../mal_RootFinderNewton.cpp',\
               '../UxRoot.cpp']
    libraryNames = ['gsl',]
    config.add_extension ('PyWindProfile', 
                          sources = sourceFiles, 
                          library_dirs = ['/opt/local/lib/'],
                          libraries = libraryNames,
#                          include_dirs = [numpy.get_include (),'/sw/include/','.'])
                          include_dirs = [numpy.get_include (), '/opt/local/include/', '.'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration (top_path='').todict()
    description_string\
        = 'Provides functions based on XSPEC local model package windprofile'
    setup (author='MAL',
           author_email='NA',
           description=description_string,
           version = '0.1',
           **config)
