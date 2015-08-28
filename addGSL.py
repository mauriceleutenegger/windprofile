#!/usr/bin/env python

import os
import sys
import platform
import re

def check_if_GSL_exists (gsl_dir) :
    if platform.system () == 'Darwin' :
        libfn = 'libgsl.dylib'
    else :
        libfn = 'libgsl.so' # check this
    if os.path.isfile (os.path.join (gsl_dir, 'lib', libfn)) :
        return True
    return False


def getGSLdir (user_gsl_dir="") :
    gsl_dirs = ['/usr', '/usr/local', '/opt/local', '/sw']
    if user_gsl_dir :
        if os.path.exists (user_gsl_dir) :
            gsl_dirs.insert (0, user_gsl_dir)
        else :
            print "User supplied GSL location %s does not exist." % user_gsl_dir
    for gsl_dir in gsl_dirs :
        if check_if_GSL_exists (gsl_dir) :
            return gsl_dir
    raise IOError ("No valid GSL installation found")

        
def getHEASOFTdir (user_HEASOFT_dir="") :
    HEASOFT_dirs = [os.getcwd ()]
    if user_HEASOFT_dir :
        if os.path.exists (user_HEASOFT_dir) :
            HEASOFT_dirs.insert (0, user_HEASOFT_dir)
        else :
            print "User supplied HEASOFT location %s does not exist."\
                % user_HEASOFT_dir
    initfile = 'headas-init.sh'
    HEADAS = ""
    try : 
        HEADAS = os.environ['HEADAS']
    except KeyError :
        pass
    if os.path.exists (HEADAS) :
        HEASOFT_dirs.append (HEADAS)
    for HEASOFT_dir in HEASOFT_dirs :
        if os.path.isfile (os.path.join (HEASOFT_dir, initfile)) :
            # this works in BUILD_DIR or $HEASOFT
            return HEASOFT_dir        
    raise IOError ("No valid HEASOFT installation found")

    
def addGSL (gsldir, heasoftdir, platform_system) :
    path_rel = '../Xspec/src/tools/initpackage/'
    fn = os.path.join (heasoftdir, path_rel, 'xspackage.tmpl')
    fn_old = os.path.join (heasoftdir, path_rel, 'xspackage.tmpl.old')
    tempfile = os.path.join (heasoftdir, path_rel, 'temp.txt')
    lib_add = " \\\n\t\t\t "
    if platform_system == "Darwin" :
        gsl_lib_dir = os.path.join (gsldir, 'lib')
        lib_add += "-L" + gsl_lib_dir + " "
    lib_add += "-lgsl -lgslcblas"
    f = open (fn, 'r')
    lines = f.readlines ()
    f.close ()
    # first test if it's been run yet
    for line in lines :
        if re.search ('gsl', line) :
            print "Detected previously modified Makefile template:"
            print fn
            print "Refusing to modify again."
            return
    g = open (tempfile, 'w')
    i = 0
    if platform_system == "Darwin" :
        # advance to HD_CXXFLAGS
        while not re.search ('HD_CXXFLAGS', lines[i]) :
            g.write (lines[i])
            i += 1
        #advance to whitespace
        j = i
        while lines[j].strip () :
            j += 1
        j -= 1
        gsl_include_dir = '-I' + os.path.join (gsldir, 'include')
        lines[j] = lines[j].rstrip (' \n\r\t\\')
        lines[j] += ' \\\n\t\t\t '
        lines[j] += gsl_include_dir + '\n'
    while not re.search ('HD_SHLIB_LIBS', lines[i]) :
        g.write (lines[i])
        i += 1
    j = i
    while lines[j].strip () :
        j += 1
    j -= 1
    if re.search ('###', lines[j]) :
        j -= 1
    lines[j] = lines[j].rstrip (' \n\r\t\\')
    lines[j] += lib_add + '\n'
    for k in range (i, len(lines)) :
        g.write (lines[k])
    g.close ()
    os.rename (fn, fn_old)
    os.rename (tempfile, fn)
    return
    
def main () :
    import argparse
    description = "\tModify Xspec initpackage Makefile to link local models "\
                  + "library to GSL,\n\tthe GNU scientific library.\n\t" \
                  + "This script will attempt to detect the GSL and HEASOFT "\
                  + "directories,\n\tor they can be supplied by the user.\n\t" \
                  + "The HEASOFT directory is checked in this order:\n\t"\
                  + "1. User supplied\n\t2. Current working directory (this"\
                  + "also works for BUILD_DIR)\n\t3. $HEADAS"
    
    parser = argparse.ArgumentParser\
             (description=description,\
              formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument ("-gsldir", help="path for gsl install", default="")
    parser.add_argument ("-heasoftdir", help="path for heasoft install",\
                         default="")
    args = parser.parse_args ()
    gsldir = getGSLdir (args.gsldir)
    heasoftdir = getHEASOFTdir (args.heasoftdir)
    platform_system = platform.system ()
    recognized_platforms = ['Darwin', 'Linux'] 
    if not platform_system in recognized_platforms :
        print "Platform/system %s not recognized" % platform_system
        print "This script is written for the following platforms:"
        for p in recognized_platforms :
            print p
        sys.exit (2)
    addGSL (gsldir, heasoftdir, platform_system)
    return

if __name__ == "__main__":
    main ()
