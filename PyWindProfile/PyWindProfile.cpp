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

#include "Python.h"
#include "numpy/arrayobject.h"
#include "xsTypes.h"
#include "../OpticalDepth.h"
#include "../Lx.h"

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
  PyObject *oP, *oZ;
  PyArrayObject *p = NULL, *z = NULL, *tau = NULL;
  Real Taustar = 0.;
  Real h = 0.;
  Real beta = 1.;
  int numerical = 0;
  int anisotropic = 0;
  int rosseland = 0;
  int  expansion = 0;
  int HeII = 0;
  int tsize[2] = {0, 0};
  size_t psize = 0;
  size_t zsize = 0;
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
  z = (PyArrayObject*) PyArray_ContiguousFromAny
    (oZ, NPY_FLOAT64, 1, 1);
  if (!p || ! z ) goto _fail;

  // make a 2-dimensional array for tau

  psize = p->dimensions[0];
  zsize = z->dimensions[0];
  tsize[0] = psize;
  tsize[1] = zsize;
  tau = (PyArrayObject*) PyArray_FromDims (2, tsize, NPY_FLOAT64);
  if (!tau) goto _fail;

  { // braces ensure that this stuff never goes out of scope due to goto !

    /* The following code is supposed to give us a way to put elements of
     a Real[][] into a 2d array stored as a numpy object*/
    PyObject ** op = (PyObject**) &tau; // pointer to tau
    double **result; 
    int nrows, ncols;
    int datatype = NPY_FLOAT64;
    PyArray_As2D (op, (char ***) &result, &nrows, &ncols, datatype);
    
    /* declare OpticalDepth object with parameters */
    OpticalDepth TAU (Taustar, h, beta, (bool) numerical, (bool) anisotropic,
		      (bool) rosseland, (bool) expansion, (bool) HeII);  

    /* calculate tau on p,z grid */
    Real* parray = (npy_float64*) p->data;
    Real* zarray = (npy_float64*) z->data;
    for (size_t i = 0; i < psize; i++) {
      for (size_t j = 0; j < zsize; j++) {
	Real temp = TAU.getOpticalDepth (parray[i],zarray[j]);
	result[i][j] = temp;
      }
    }
    // should probably free memory, but this function seems to crash
    //    PyArray_Free (*op, (char *) *result); //get rid of the 2d baggage
  } // end brace protection

  /* clean up */
  Py_DECREF (z);
  Py_DECREF (p);
  return PyArray_Return (tau);
 _fail:
  Py_XDECREF (z);
  Py_XDECREF (p);
  Py_XDECREF (tau);
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
  size_t xsize = 0;
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
  isAnisotropic = (bool) aniso;
  isRosseland = (bool) ross;
  isExpansion = (bool) exp;
  isHeII = (bool) HeII;

 // make array object from python numpy array
  x = (PyArrayObject*) PyArray_ContiguousFromAny
    (oX, NPY_FLOAT64, 1, 1);
  if (!x) goto _fail; 
  xsize = x->dimensions[0];

  // and make a flux array object with same size as x
  flux = (PyArrayObject*) PyArray_FromDims
    (1, (int*) &xsize, NPY_FLOAT64);
  if (!flux) goto _fail; 

  // for your protection
  {
  Real* xarray = (npy_float64*) x->data;
  Real* fluxarray = (npy_float64*) flux->data;

   // braces protect the objects from goto
    // instantiate Lx and other objects
    V = new Velocity (beta);
    He = new HeLikeRatio (); // default case is harmless
    RS = new ResonanceScattering (Tau0Star, betaSobolev, isOpticallyThick, V);
    Tau = new OpticalDepth 
      (TauStar, h, beta, isNumerical, isAnisotropic, isRosseland, isExpansion, false);
    if (isHeII) {
      TauHeII = new OpticalDepth 
	(TauStar, h, beta, true , isAnisotropic, isRosseland, isExpansion, true);
      // Note that TauHeII needs to have numerical set to true to work properly.
      lx = new Lx (q, U0, Umin, beta, kappaRatio, wResonance, V, He, RS, Tau, TauHeII);
    } else {
      lx = new Lx (q, U0, Umin, beta, wResonance, V, He, RS, Tau);
    }
    // calculate flux
    for (size_t i = 0; i < xsize; i++) {
      fluxarray[i] = lx->getLx (xarray[i]);
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
  {"Lx", Py_Lx, METH_VARARGS, "Calculate Lx(x)"},
  {NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initPyWindProfile (void)
{
  import_array () (void) Py_InitModule ("PyWindProfile",PyWindProfileMethods);
}
