/***************************************************************************
    PyWindProfile.cpp   - Provides a python interface for Lx and OpticalDepth.

                             -------------------
    begin				: Winter 2009
    copyright			: (C) 2009 by Maurice Leutenegger
    email				: maurice.a.leutenegger@nasa.gov
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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// disallow deprecated numpy methods

#include "Python.h"
#include "numpy/arrayobject.h"
#include "xsTypes.h"
#include "../../OpticalDepth.h"
#include "../../Lx.h"
#include "../../RAD_OpticalDepth.h"

static PyObject* Py_OpticalDepth (PyObject* obj, PyObject* args)
{
  Real p = 0.;
  Real z = 0.;
  Real Taustar = 0.;
  Real h = 0.;
  Real beta = 1.;
  int numerical = 0;
  int anisotropic = 0;
  int rosseland = 0;
  int  expansion = 0;
  Real tau = 1.e6;
  // THERE IS PROBABLY A BUG HERE
  if (!PyArg_ParseTuple (args, "dddddiiii", &p, &z, &Taustar, &h, &beta, 
			 &numerical, &anisotropic, &rosseland, &expansion)) {
    PyErr_SetString (PyExc_ValueError, 
		     "OpticalDepth: Invalid number of parameters.");
    return NULL;
  }
  
  /* declare OpticalDepth object with parameters */
  OpticalDepth TAU (Taustar, h, beta, (bool) numerical, (bool) anisotropic,
		    (bool) rosseland, (bool) expansion);  
  
  tau = TAU.getOpticalDepth (p, z);
  
  return Py_BuildValue ("d", tau);
  
}

static PyObject* Py_OpticalDepth2d (PyObject* obj, PyObject* args)
{
  PyObject *oP = NULL, *oZ = NULL;
  PyArrayObject *p = NULL, *z = NULL, *tau = NULL;
  Real Taustar = 0.;
  Real h = 0.;
  Real beta = 1.;
  Real tempt = 0.;
  Real tempp = 0.;
  Real tempz = 0.;
  int numerical = 0;
  int anisotropic = 0;
  int rosseland = 0;
  int  expansion = 0;
  int HeII = 0;
  npy_intp tdims[2] = {0, 0};
  npy_intp psize = 0;
  npy_intp zsize = 0;
  // set up descriptor for array allocation
  PyArray_Descr *tau_descriptor = PyArray_DescrFromType(NPY_FLOAT64);
    
  // THERE IS PROABABLY A BUG HERE
  if (!PyArg_ParseTuple (args, "OOdddiiiii", &oP, &oZ, &Taustar, &h, &beta, 
			 &numerical, &anisotropic, &rosseland, &expansion, &HeII)) {
    PyErr_SetString (PyExc_ValueError, 
		     "OpticalDepth: Invalid number of parameters.");
    goto _fail;
  }
  // make array objects from python numpy arrays
  p = (PyArrayObject*) PyArray_ContiguousFromAny
    (oP, NPY_FLOAT64, 1, 1);
  if (!p) {
    PyErr_SetString (PyExc_ValueError,
                     "Py_OpticalDepth2d: Failed to allocate p array.");
    goto _fail;
  }
  z = (PyArrayObject*) PyArray_ContiguousFromAny
    (oZ, NPY_FLOAT64, 1, 1);
  if (!z) {
    PyErr_SetString (PyExc_ValueError,
                     "Py_OpticalDepth2d: Failed to allocate z array.");
    goto _fail;
  }

  // get and set array dimensions
  psize = PyArray_DIM (p, 0);
  zsize = PyArray_DIM (z, 0);
  tdims[0] = psize;
  tdims[1] = zsize;
  // make a 2-dimensional array for tau, initialized to zeros
  tau = (PyArrayObject*) PyArray_Zeros
    (2, tdims, tau_descriptor, 0);
  if (!tau) {
    PyErr_SetString (PyExc_ValueError,
                     "Py_OpticalDepth2d: Failed to allocate tau array.");
    goto _fail;
  }

  { // braces ensure that this stuff never goes out of scope due to goto !
    
    /* declare OpticalDepth object with parameters */
    OpticalDepth TAU (Taustar, h, beta, (bool) numerical, (bool) anisotropic,
		      (bool) rosseland, (bool) expansion, (bool) HeII);  

    /* calculate tau on p,z grid */
    /*
      Here we use GETPTR1 to return a void* pointer to the array at element i;
      this must then be cast to npy_float64* so we know the type;
      and then the whole thing must be dereferenced to get the actual value.
      For the 2D array tau, we use GETPTR2 in the same way.
    */
    for (npy_intp i = 0; i < psize; i++) {
      for (npy_intp j = 0; j < zsize; j++) {
        tempp = *((npy_float64*) PyArray_GETPTR1 (p, i));
        tempz = *((npy_float64*) PyArray_GETPTR1 (z, j));
        tempt = TAU.getOpticalDepth (tempp, tempz);
        *((npy_float64*) PyArray_GETPTR2 (tau, i, j)) = tempt;
      }
    }
  } // end brace protection
  /* clean up */
  Py_XDECREF (z);
  Py_XDECREF (p);
  return PyArray_Return (tau);
 _fail:
  cout << "optical depth failing\n" << flush;
  Py_XDECREF (z);
  Py_XDECREF (p);
  Py_XDECREF (tau);
  return NULL;
}

static PyObject* Py_RAD_OpticalDepth (PyObject* obj, PyObject* args)
{
  Real p = 0.;
  Real z = 0.;
  Real Tau0 = 0.;
  Real beta = 1.;
  Real deltaE = 0.;
  Real gamma = 0.;
  Real vinfty = 0.;
  Real tau = 1.e6;
  Velocity* V = NULL;
  // THERE IS PROBABLY A BUG HERE
  if (!PyArg_ParseTuple (args, "ddddddd", &p, &z, &Tau0, &beta, 
			 &deltaE, &gamma, &vinfty)) {
    PyErr_SetString (PyExc_ValueError, 
		     "RAD_OpticalDepth: Invalid number of parameters.");
    return NULL;
  }

  /* declare velocity object */
  V = new Velocity (beta, 0.);
  /* declare OpticalDepth object with parameters */
  RAD_OpticalDepth TAU (V, deltaE, gamma, Tau0, vinfty);

  tau = TAU.getOpticalDepth (p, z);

  delete V;
  V = NULL;
  
  return Py_BuildValue ("d", tau);
  
}

// edit me
// add a zarray input and return an array
static PyObject* Py_RAD_OpticalDepth_integrand (PyObject* obj, PyObject* args)
{
  PyObject *oZ = NULL;
  PyArrayObject *z = NULL, *integrand = NULL;
  Real p = 0.;
  Real z0 = 0.;
  Real Tau0 = 1.; // dummy
  Real beta = 1.;
  Real deltaE = 0.;
  Real gamma = 0.;
  Real vinfty = 0.;
  //Real tau = 1.e6;
  Velocity* V = NULL;

  npy_intp zsize = 0;
  // set up descriptor for array allocation
  PyArray_Descr *integrand_descriptor = PyArray_DescrFromType(NPY_FLOAT64);
  
  // THERE IS PROBABLY A BUG HERE
  if (!PyArg_ParseTuple (args, "Odddddd", &oZ, &p, &z0, &beta, 
			 &deltaE, &gamma, &vinfty)) {
    PyErr_SetString (PyExc_ValueError, 
		     "RAD_OpticalDepth: Invalid number of parameters.");
    return NULL;
  }

  // make array objects from python numpy arrays
  z = (PyArrayObject*) PyArray_ContiguousFromAny
    (oZ, NPY_FLOAT64, 1, 1);
  if (!z) {
    PyErr_SetString (PyExc_ValueError,
                     "Py_OpticalDepth2d: Failed to allocate z array.");
    goto _fail;
  }

  // get array dimension
  zsize = PyArray_DIM (z, 0);

  // and make an integrand array object with same size as z
  integrand = (PyArrayObject*) PyArray_Zeros
    (1, &zsize, integrand_descriptor, 0);
  if (!integrand) {
    PyErr_SetString (PyExc_ValueError,
                     "RAD_OpticalDepth_integrand: Failed to allocate integrand array.");
    goto _fail;
  }

  { // braces protect objects from goto
    /* declare velocity object */
    V = new Velocity (beta, 0.);
    /* declare OpticalDepth object with parameters */
    RAD_OpticalDepth TAU (V, deltaE, gamma, Tau0, vinfty);

    TAU.initialize (p,z0);

    // get integrand
    /*
      Here we use GETPTR1 to return a void* pointer to the array at element i;
      this must then be cast to npy_float64* so we know the type;
      and then the whole thing must be dereferenced to get the actual value.
    */
    for (npy_intp i = 0; i < zsize; i++) {
      *((npy_float64*) PyArray_GETPTR1 (integrand, i)) =
        TAU.integrand (*((npy_float64*) PyArray_GETPTR1 (z, i)));
    }
  
    delete V;
    V = NULL;

  } // end protect goto
  
  Py_DECREF (z);
  return PyArray_Return (integrand);

   _fail:
  Py_XDECREF (z);
  Py_XDECREF (integrand);
  return NULL;
}

static PyObject* Py_Lx (PyObject* obj, PyObject* args)
{
  // pyobjects
  PyObject *oX;
  PyArrayObject *x = NULL, *flux = NULL;
  // windprofile objects
  Lx* lx;
  Velocity* V;
  HeLikeRatio* He;
  ResonanceScattering* RS;
  OpticalDepth* Tau;
  OpticalDepth* TauHeII;
  // preset parameters
  Real q = 0.;
  Real U0 = 0.5;
  Real Umin = 0.;
  Real beta = 1.;
  Real TauStar = 0.;
  Real h = 0.;
  Real Tau0Star = 0.;
  Real betaSobolev = 0.;
  Real kappaRatio = 0.;
  bool isOpticallyThick = false;
  bool isNumerical = false;
  bool isAnisotropic = false;
  bool isProlate = false;
  bool isRosseland = false;
  bool isExpansion = false;
  bool isHeII = false;
  // somehow the behavior of PyArg_ParseTuple changed,
  // and we need to use ints instead of bools
  int thick = 0;
  int num = 0;
  int aniso = 0;
  int ross = 0;
  int exp = 0;
  int HeII = 0;
  npy_intp xsize = 0;
  PyArray_Descr *flux_descriptor = PyArray_DescrFromType(NPY_FLOAT64);
  // get arguments from python
  if (!PyArg_ParseTuple 
      (args, "Odddddddddiiiiii", &oX, &q, &U0, &Umin, &beta, &TauStar, &h, 
       &Tau0Star, &betaSobolev, &kappaRatio,
       &thick, &num, &aniso, 
       &ross, &exp, &HeII)) {
    PyErr_SetString (PyExc_ValueError, 
                     "Lx: Invalid number of parameters.");
    goto _fail;
  }
  isOpticallyThick = (bool) thick;
  isNumerical = (bool) num;
  if (aniso == 0) {
    isAnisotropic = false;
    isProlate = false;
  } else if (aniso == 1) {
    isAnisotropic = true;
    isProlate = false;
  } else if (aniso == 2) {
    isAnisotropic = true;
    isProlate = true;
  } else {
    cerr << "Lx: got unexpected value for aniso:\t" << aniso << endl;
  }
  isRosseland = (bool) ross;
  isExpansion = (bool) exp;
  isHeII = (bool) HeII;

  // make array object from python numpy array
  x = (PyArrayObject*) PyArray_ContiguousFromAny
    (oX, NPY_FLOAT64, 1, 1);

  if (!x) {
    PyErr_SetString (PyExc_ValueError,
                     "Lx: Failed to allocate x array.");
    goto _fail;
  }
  xsize = PyArray_DIM (x, 0);

  // and make a flux array object with same size as x
  flux = (PyArrayObject*) PyArray_Zeros
    (1, &xsize, flux_descriptor, 0);
  //flux = (PyArrayObject*) PyArray_SimpleNew
  //  (1, (long const *) &xsize, NPY_FLOAT64);
  if (!flux) {
    PyErr_SetString (PyExc_ValueError,
                     "Lx: Failed to allocate flux array.");
    goto _fail;
  }

  // for your protection
  {    // braces protect the objects from goto
    // instantiate Lx and other objects
    V = new Velocity (beta);
    He = new HeLikeRatio (); // default case is harmless
    RS = new ResonanceScattering (Tau0Star, betaSobolev, isOpticallyThick, V);
    Tau = new OpticalDepth 
      (TauStar, h, beta, isNumerical, isAnisotropic, isProlate, isRosseland, isExpansion, false);
    if (isHeII) {
      TauHeII = new OpticalDepth 
        (TauStar, h, beta, true , isAnisotropic, isProlate, isRosseland, isExpansion, true);
      // Note that TauHeII needs to have numerical set to true to work properly.
      lx = new Lx (q, U0, Umin, beta, kappaRatio, wResonance, V, He, RS, Tau, TauHeII);
    } else {
      lx = new Lx (q, U0, Umin, beta, wResonance, V, He, RS, Tau);
    }
    // calculate flux
    /*
      Here we use GETPTR1 to return a void* pointer to the array at element i;
      this must then be cast to npy_float64* so we know the type;
      and then the whole thing must be dereferenced to get the actual value.
    */
    for (npy_intp i = 0; i < xsize; i++) {
      *((npy_float64*) PyArray_GETPTR1 (flux, i)) =
        lx->getLx (*((npy_float64*) PyArray_GETPTR1 (x, i)));
    }
    delete lx;
    delete Tau;
    if (isHeII) {
      delete TauHeII;
    }
    delete RS;
    delete He;
    delete V;
  } // end protect braces
  Py_DECREF (x);
  return PyArray_Return (flux);
 _fail:
  Py_XDECREF (x);
  Py_XDECREF (flux);
  return NULL;
}

static PyMethodDef PyWindProfileMethods[] = {
  {"OpticalDepth", Py_OpticalDepth, METH_VARARGS, "Calculate t(p,z), scalar"},
  {"OpticalDepth2d", Py_OpticalDepth2d, METH_VARARGS, "Calculate t(p,z), 2d"},
  {"RAD_OpticalDepth", Py_RAD_OpticalDepth, METH_VARARGS,
   "Calculate RAD t(p,z)"},
  {"RAD_OpticalDepth_integrand", Py_RAD_OpticalDepth_integrand, METH_VARARGS,
   "Calculate RAD integrand (p,z0) on z array"},
  {"Lx", Py_Lx, METH_VARARGS, "Calculate Lx(x)"},
  {NULL, NULL, 0, NULL} /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
/*
  This structure defines the module, and makes use of the methods structure.
 */
static struct PyModuleDef PyWindProfileDef =
{
    PyModuleDef_HEAD_INIT,
    "PyWindProfile", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    PyWindProfileMethods
};

/*
  This function creates the module using the module structure.
 */
PyMODINIT_FUNC PyInit_PyWindProfile(void)
{
  PyObject *module;
  module = PyModule_Create(&PyWindProfileDef);
  if(module==NULL) return NULL;
  /* IMPORTANT: this must be called */
  import_array();
  if (PyErr_Occurred()) return NULL;
  return module;
}
#else
// This is the py2.7 way
PyMODINIT_FUNC initPyWindProfile (void)
{
  import_array () (void) Py_InitModule ("PyWindProfile",PyWindProfileMethods);
}
#endif
