/***************************************************************************
    WindAbsorption.cpp   - Windabsorption: computes fraction of transmitted 
                           light over the whole wind, weighted by emission 
                           measure.
                           FractionalWindEmission: computes unrenormalized
                           emission from a fraction of the wind, from 
                           infinite radius up to R (allows to make a 
                           cumulative distribution of wind emission)
                           AngleAveragedTransmission: computes AAT as a
                           function of radius.
                           Implemented to be called from python.
                           There is a python script to run it and generate
                           a fits file for use in the XSPEC model windtabs.

                             -------------------
    begin				: Summer 2008
    copyright			: (C) 2009 by Maurice Leutenegger
    email				: Maurice.A.Leutenegger@nasa.gov
 ***************************************************************************/
 /* This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */


#include "xsTypes.h"
#include "../IntegratedLuminosity.h"
#include "../AngleAveragedTransmission.h"
#include "../OpticalDepth.h"
#include "../Utilities.h"
#include <gsl/gsl_const_cgsm.h>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Python.h"
#include "numpy/arrayobject.h"

static const Real CONST_HC_KEV_A = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 
    GSL_CONST_CGSM_SPEED_OF_LIGHT * 1.e5 / GSL_CONST_CGSM_ELECTRON_VOLT;

 // arguments are: 
 // q, u0, beta, h, isNumerical, isAnisotropic, isRosseland, taustar[]
 // return value is:
 // transmission[]

static PyObject* Py_WindAbsorption (PyObject* obj, PyObject* args)
{
  // pyobjects; taustar is input array, transmission is output array
  PyObject *oTAUSTAR;
  PyArrayObject *Taustar = NULL, *Transmission = NULL;
  // preset arguments
  Real q (0.);
  Real u0 (0.5);
  Real beta (1.0); 
  Real h (0.); 
  int isNumerical (0);
  int isAnisotropic (0);
  int isRosseland (0);
  // declare variables
  size_t Tsize = 0;
  // get arguments from python
  if (!PyArg_ParseTuple 
      (args, "Oddddiii", &oTAUSTAR, &q, &u0, &beta, &h, 
       &isNumerical, & isAnisotropic, &isRosseland)) {
    PyErr_SetString (PyExc_ValueError, 
		     "WindAbsorption: Invalid number of parameters.");
    goto _fail;
  }
 // make array object from python numpy array
  Taustar = (PyArrayObject*) PyArray_ContiguousFromAny
    (oTAUSTAR, NPY_FLOAT64, 1, 1);
  if (!Taustar) goto _fail; 
  Tsize = Taustar->dimensions[0];
  // and make an output array object with same size
  Transmission = (PyArrayObject*) PyArray_FromDims
    (1, (int*) &Tsize, NPY_FLOAT64);
  if (!Transmission) goto _fail; 
  // for your protection
  {
    Real* Taustar_Carray = (npy_float64*) Taustar->data;
    Real* Transmission_Carray = (npy_float64*) Transmission->data;

   // calculate answer; protected from goto by braces
    Velocity* V;
    OpticalDepth* Tau;
    AngleAveragedTransmission* T;
    IntegratedLuminosity* L;
    V = new Velocity (beta, 0.);
    bool isHeII = false; // fixed for now
    Tau = new OpticalDepth (0., h , beta, (bool) isNumerical, (bool) isAnisotropic, (bool) isRosseland, false, isHeII);
    T = new AngleAveragedTransmission (Tau);
    L = new IntegratedLuminosity (q, u0, 0., V, T);
    T->setEpsRel (1.e-4); // default accuracy
    L->setEpsRel (1.e-4);
    // do calculation:
    L->setTransparentCore (true); 
    Real IntrinsicLuminosity = L->getLuminosity ();
    L->setTransparentCore (false);
    for (size_t i = 0; i < Tsize; i++) {
      Tau->setParameters (Taustar_Carray[i], h);
      Real Luminosity = L->getLuminosity ();
      Transmission_Carray[i] = Luminosity / IntrinsicLuminosity;
    }
    // clean up
    delete L;
    delete T;
    delete Tau;
    delete V;
  } // end protect braces

  Py_DECREF (Taustar);
  return PyArray_Return (Transmission);
 _fail:
  Py_XDECREF (Taustar);
  Py_XDECREF (Transmission);
  return NULL;

}


static PyObject* Py_WindAbsorptionHeII (PyObject* obj, PyObject* args)
{
  // pyobjects; taustar is input array, transmission is output array
  PyObject *oTAUSTAR;
  PyArrayObject *Taustar = NULL, *Transmission = NULL;
  // pyobjects; taustar is input array, transmission is output array
  PyObject *oKAPPARATIO;
  PyArrayObject *kappaRatio = NULL;
  // preset arguments
  Real q (0.);
  Real u0 (0.5);
  Real beta (1.0); 
  Real h (0.); 
  int isNumerical (0);
  int isAnisotropic (0);
  int isRosseland (0);
  // declare variables
  size_t Tsize = 0;
  size_t Ksize = 0;
  size_t TrSize[2] = {0, 0};
  // get arguments from python
  if (!PyArg_ParseTuple 
      (args, "OOddddiii", &oTAUSTAR, &oKAPPARATIO, &q, &u0, &beta, &h,
       &isNumerical, & isAnisotropic, &isRosseland)) {
    PyErr_SetString (PyExc_ValueError, 
		     "WindAbsorption: Invalid number of parameters.");
    goto _fail;
  }
 // make array object from python numpy array
  Taustar = (PyArrayObject*) PyArray_ContiguousFromAny
    (oTAUSTAR, NPY_FLOAT64, 1, 1);
  if (!Taustar) goto _fail; 
  Tsize = Taustar->dimensions[0];
 // make array object from python numpy array
  kappaRatio = (PyArrayObject*) PyArray_ContiguousFromAny
    (oKAPPARATIO, NPY_FLOAT64, 1, 1);
  if (!kappaRatio) goto _fail; 
  Ksize = kappaRatio->dimensions[0];
  TrSize[0] = Tsize;
  TrSize[1] = Ksize;
  // and make a 2D output array object with correct size
  Transmission = (PyArrayObject*) PyArray_FromDims
    (2, (int*) &TrSize, NPY_FLOAT64);
  if (!Transmission) goto _fail; 
  // for your protection
  {

    /* The following code is supposed to give us a way to put elements of
     a Real[][] into a 2d array stored as a numpy object*/
    PyObject ** op = (PyObject**) &Transmission; // pointer to Transmission
    double **result; 
    int nrows, ncols;
    int datatype = NPY_FLOAT64;
    PyArray_As2D (op, (char ***) &result, &nrows, &ncols, datatype);

    Real* Taustar_Carray = (npy_float64*) Taustar->data;
    Real* kappaRatio_Carray = (npy_float64*) kappaRatio->data;
    //    Real* Transmission_Carray = (npy_float64*) Transmission->data;

   // calculate answer; protected from goto by braces
    Velocity* V;
    OpticalDepth* Tau;
    OpticalDepth* TauHeII;
    AngleAveragedTransmission* T;
    IntegratedLuminosity* L;
    V = new Velocity (beta, 0.);
    bool isHeII = false; // fixed for now
    Tau = new OpticalDepth (0., h , beta, (bool) isNumerical, (bool) isAnisotropic, (bool) isRosseland, false, false);
    TauHeII = new OpticalDepth (0., h , beta, true, (bool) isAnisotropic, (bool) isRosseland, false, true);
    T = new AngleAveragedTransmission (Tau, TauHeII, 0.); 
    // start with kappaRatio = 0., will change it later
    L = new IntegratedLuminosity (q, u0, 0., V, T);
    T->setEpsRel (1.e-4); // default accuracy
    L->setEpsRel (1.e-4);
    // do calculation:
    L->setTransparentCore (true);  
    Real IntrinsicLuminosity = L->getLuminosity ();
    L->setTransparentCore (false);
    for (size_t i = 0; i < Tsize; i++) {
      Tau->setParameters (Taustar_Carray[i], h); // set tau*
      TauHeII->setParameters (Taustar_Carray[i], h);
      for (size_t j = 0; j < Ksize; j++) {
	T->setKappaRatio (kappaRatio_Carray[j]); // set kappaRatio
	Real Luminosity = L->getLuminosity ();
	//	Transmission_Carray[i] = Luminosity / IntrinsicLuminosity;
	result[i][j] = Luminosity / IntrinsicLuminosity;
	std::cout << Taustar_Carray[i] << " " << kappaRatio_Carray[j] << "\n";
      }
    }
    // clean up
    delete L;
    delete T;
    delete Tau;
    delete TauHeII;
    delete V;
  } // end protect braces

  Py_DECREF (Taustar);
  return PyArray_Return (Transmission);
 _fail:
  Py_XDECREF (Taustar);
  Py_XDECREF (Transmission);
  return NULL;

}

static PyObject* Py_WindAbsorptionHeII_Fixed_Kappa 
(PyObject* obj, PyObject* args)
{
  // pyobjects; taustar is input array, transmission is output array
  PyObject *oTAUSTAR;
  PyArrayObject *Taustar = NULL, *Transmission = NULL;
  // preset arguments
  Real kappaRatio (0.);
  Real q (0.);
  Real u0 (0.5);
  Real beta (1.0); 
  Real h (0.); 
  int isNumerical (0);
  int isAnisotropic (0);
  int isRosseland (0);
  // declare variables
  size_t Tsize = 0;
  size_t Ksize = 0;
  size_t TrSize[2] = {0, 0};
  // get arguments from python
  if (!PyArg_ParseTuple 
      (args, "Odddddiii", &oTAUSTAR, &kappaRatio, &q, &u0, &beta, &h,
       &isNumerical, & isAnisotropic, &isRosseland)) {
    PyErr_SetString (PyExc_ValueError, 
		     "WindAbsorption: Invalid number of parameters.");
    goto _fail;
  }

 // make array object from python numpy array
  Taustar = (PyArrayObject*) PyArray_ContiguousFromAny
    (oTAUSTAR, NPY_FLOAT64, 1, 1);
  if (!Taustar) goto _fail; 
  Tsize = Taustar->dimensions[0];
  // and make an output array object with same size
  Transmission = (PyArrayObject*) PyArray_FromDims
    (1, (int*) &Tsize, NPY_FLOAT64);
  if (!Transmission) goto _fail; 
  // for your protection
  {

    /* The following code is supposed to give us a way to put elements of
     a Real[][] into a 2d array stored as a numpy object*/
    //    int nrows, ncols;
    //    int datatype = NPY_FLOAT64;
    //    PyArray_As2D (op, (char ***) &result, &nrows, &ncols, datatype);

    Real* Taustar_Carray = (npy_float64*) Taustar->data;
    //    Real* kappaRatio_Carray = (npy_float64*) kappaRatio->data;
    Real* Transmission_Carray = (npy_float64*) Transmission->data;

   // calculate answer; protected from goto by braces
    Velocity* V;
    OpticalDepth* Tau;
    OpticalDepth* TauHeII;
    AngleAveragedTransmission* T;
    IntegratedLuminosity* L;
    V = new Velocity (beta, 0.);
    bool isHeII = false; // fixed for now
    Tau = new OpticalDepth (0., h , beta, (bool) isNumerical, (bool) isAnisotropic, (bool) isRosseland, false, false);
    TauHeII = new OpticalDepth (0., h , beta, true, (bool) isAnisotropic, (bool) isRosseland, false, true);
    T = new AngleAveragedTransmission (Tau, TauHeII, 0.); 
    // start with kappaRatio = 0., will change it later
    L = new IntegratedLuminosity (q, u0, 0., V, T);
    T->setEpsRel (1.e-4); // default accuracy
    L->setEpsRel (1.e-4);
    // do calculation:
    L->setTransparentCore (true);  
    Real IntrinsicLuminosity = L->getLuminosity ();
    L->setTransparentCore (false);
    T->setKappaRatio (kappaRatio);

    for (size_t i = 0; i < Tsize; i++) {
      Tau->setParameters (Taustar_Carray[i], h);
      TauHeII->setParameters (Taustar_Carray[i], h);
      Real Luminosity = L->getLuminosity ();
      Transmission_Carray[i] = Luminosity / IntrinsicLuminosity;
    }
    // clean up
    delete L;
    delete T;
    delete Tau;
    delete TauHeII;
    delete V;
  } // end protect braces

  Py_DECREF (Taustar);
  return PyArray_Return (Transmission);
 _fail:
  Py_XDECREF (Taustar);
  Py_XDECREF (Transmission);
  return NULL;

}

 // arguments are: 
 // q, Taustar, beta, h, isNumerical, isAnisotropic, isRosseland, u[]
 // return value is:
 // flux[] (not renormalized)
static PyObject* Py_FractionalWindEmission (PyObject* obj, PyObject* args)
{
  // pyobjects; taustar is input array, transmission is output array
  PyObject *oU;
  PyArrayObject *u = NULL, *Flux = NULL;
  // preset arguments
  Real q (0.);
  Real Taustar (1.);
  Real u0 (0.5);
  Real beta (1.0); 
  Real h (0.); 
  int isNumerical (0);
  int isAnisotropic (0);
  int isRosseland (0);
  // declare variables
  size_t Fsize = 0;
  // get arguments from python
  if (!PyArg_ParseTuple 
      (args, "Oddddiii", &oU, &q, &Taustar, &beta, &h, &isNumerical, 
       &isAnisotropic, &isRosseland)) {
    PyErr_SetString (PyExc_ValueError, 
		     "WindAbsorption: Invalid number of parameters.");
    goto _fail;
  }
 // make array object from python numpy array
  u = (PyArrayObject*) PyArray_ContiguousFromAny (oU, NPY_FLOAT64, 1, 1);
  if (!u) goto _fail; 
  Fsize = u->dimensions[0];
  // and make an output array object with same size
  Flux = (PyArrayObject*) PyArray_FromDims (1, (int*) &Fsize, NPY_FLOAT64);
  if (!Flux) goto _fail; 
  // for your protection
  {
    Real* U_Carray = (npy_float64*) u->data;
    Real* Flux_Carray = (npy_float64*) Flux->data;

   // calculate answer; protected from goto by braces
    Velocity* V;
    OpticalDepth* Tau;
    AngleAveragedTransmission* T;
    IntegratedLuminosity* L;
    V = new Velocity (beta, 0.);
    bool isHeII = false; // fixed for now
    Tau = new OpticalDepth (Taustar, h , beta, (bool) isNumerical, (bool) isAnisotropic, (bool) isRosseland, false, isHeII);
    T = new AngleAveragedTransmission (Tau);
    L = new IntegratedLuminosity (q, u0, 0., V, T);
    T->setEpsRel (1.e-4); // default accuracy
    L->setEpsRel (1.e-4);
    // do calculation:
    L->setTransparentCore (false);
    for (size_t i = 0; i < Fsize; i++) {
      L->setU0 (U_Carray[i]);
      Real Luminosity = L->getLuminosity ();
      Flux_Carray[i] = Luminosity;
    }
    // clean up
    delete L;
    delete T;
    delete Tau;
    delete V;
  } // end protect braces

  Py_DECREF (u);
  return PyArray_Return (Flux);
 _fail:
  Py_XDECREF (u);
  Py_XDECREF (Flux);
  return NULL;

}


// arguments are: 
 // Taustar, beta, h, isNumerical, isAnisotropic, isRosseland, u[]
 // return value is:
 // transmission[]
static PyObject* Py_AngleAveragedTransmission (PyObject* obj, PyObject* args)
{
  // pyobjects; taustar is input array, transmission is output array
  PyObject *oU;
  PyArrayObject *u = NULL, *Transmission = NULL;
  // preset arguments
  Real Taustar (0.);
  Real kappaRatio (0.); // ratio of kappaHeII to kappaMetalIons
  Real beta (1.0); 
  Real h (0.); 
  int isNumerical (0);
  int isAnisotropic (0);
  int isRosseland (0);
  int isHeII (0);
  // declare variables
  size_t Tsize = 0;
  // get arguments from python
  if (!PyArg_ParseTuple 
      (args, "Oddddiiii", &oU, &Taustar, &kappaRatio, &beta, &h, 
       &isNumerical, & isAnisotropic, &isRosseland, &isHeII)) {
    char errorString[] =
      "AngleAveragedTransmission: Invalid number of parameters.";
    PyErr_SetString (PyExc_ValueError, errorString);
    goto _fail;
  }
 // make array object from python numpy array
  u = (PyArrayObject*) PyArray_ContiguousFromAny
    (oU, NPY_FLOAT64, 1, 1);
  if (!u) goto _fail; 
  Tsize = u->dimensions[0];
  // and make an output array object with same size
  Transmission = (PyArrayObject*) PyArray_FromDims
    (1, (int*) &Tsize, NPY_FLOAT64);
  if (!Transmission) goto _fail; 
  // for your protection
  {
    Real* u_Carray = (npy_float64*) u->data;
    Real* Transmission_Carray = (npy_float64*) Transmission->data;

    // calculate answer; protected from goto by braces
    //    Velocity* V;
    OpticalDepth* Tau;
    OpticalDepth* TauHeII;
    AngleAveragedTransmission* T;
    //    V = new Velocity (beta, 0.);
    Tau = new OpticalDepth (Taustar, h , beta, (bool) isNumerical, (bool) isAnisotropic, (bool) isRosseland, false, false);
    if (isHeII) {
      TauHeII = new OpticalDepth (Taustar, h, beta, true, (bool) isAnisotropic, (bool) isRosseland, false, true);
      T = new AngleAveragedTransmission (Tau, TauHeII, kappaRatio);
    }
    else {
      T = new AngleAveragedTransmission (Tau);
    }
    T->setEpsRel (1.e-4); // default accuracy
    // do calculation:
    for (size_t i = 0; i < Tsize; i++) {
      Transmission_Carray[i] = T->getTransmission (u_Carray[i]);
    }
    // clean up
    delete T;
    delete Tau;
    //    delete V;
    if (isHeII) {
      delete TauHeII;
    }
  } // end protect braces

  Py_DECREF (u);
  return PyArray_Return (Transmission);
 _fail:
  Py_XDECREF (u);
  Py_XDECREF (Transmission);
  return NULL;
}

static PyMethodDef windabsorptionMethods[] = {
  {"WindAbsorption", Py_WindAbsorption, METH_VARARGS, 
   "Calculate transmission of whole wind"},
  {"WindAbsorptionHeII", Py_WindAbsorptionHeII, METH_VARARGS, 
   "Calculate transmission of whole wind, including HeII recombination"},
  {"WindAbsorptionHeII_Fixed_Kappa", Py_WindAbsorptionHeII_Fixed_Kappa, METH_VARARGS, 
   "Calculate transmission of whole wind, including HeII recombination - 1D output with fixed kappa"},
  {"FractionalWindEmission", Py_FractionalWindEmission, METH_VARARGS,
   "Calculate unrenormalized emission from fraction of wind"},
  {"AngleAveragedTransmission", Py_AngleAveragedTransmission, METH_VARARGS, 
   "Calculate Transmission as a function of radius"},
  {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC initwindabsorption (void)
{
  import_array () (void) Py_InitModule ("windabsorption", windabsorptionMethods);
}
