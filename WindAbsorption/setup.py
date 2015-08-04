#!/usr/bin/env python
import numpy

def configuration (parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration ('windabsorption', parent_package, top_path)
    sourceFiles\
        = ['WindAbsorption.cpp',\
               '../Utilities.cpp',\
               '../Porosity.cpp',\
               '../NumericalOpticalDepth.cpp',\
               '../NumericalOpticalDepthZ.cpp',\
               '../NumericalOpticalDepthU.cpp',\
               '../AnalyticOpticalDepth.cpp',\
               '../Series.cpp',\
               '../IsotropicSeries.cpp',\
               '../SmoothA1.cpp',\
               '../OpticalDepth.cpp',\
               '../mal_integration.cpp',\
               '../AngleAveragedTransmission.cpp',\
               '../IntegratedLuminosity.cpp']
    libraryNames = ['gsl',]
    config.add_extension ('windabsorption', 
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
        = 'C/C++ function to calculate emission-measure weighted transmission of wind'
    setup (author='MAL',
           author_email='NA',
           description=description_string,
           version = '0.1',
           **config)
