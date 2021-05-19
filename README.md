# windprofile
X-ray emission line profiles for OB stars

MAIN DOCUMENTATION IS HERE:
https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/windprof.html

Supplemental documentation for models using tabulated opacities and transmissions (e.g. windtabs):

These are the relevant keywords to set using xset:

keyword                default value

WINDTABSDIRECTORY      ./
set to give path where data files are found

KAPPAFILENAME          kappa.fits
file with opacities (fixed abundance)

KAPPAHEIIFILENAME      kappaHeII.fits
file with opacities with He II recombined

KAPPAZFILENAME         kappa.fits
file with opacities split by element (for variable abundance)

TRANSMISSIONFILENAME   tau_transmission.fits
file with tabulated transmission vs tau

TRANSMISSION2DFILENAME tau_transmission_HeII.fits
file with tabulated transmission vs both tau and kapparatio

HEII
if this is set to 1, will allow for HeII recombination using KAPPAHEIIFILENAME and TRANSMISSION2DFILENAME

SAVEKAPPAZ
if this is set to 1, will write out the kappa that is calculated from a variable abundance model to an ascii file

KAPPAZOUTFILE          kappaZ.txt
this is the file that it's written to