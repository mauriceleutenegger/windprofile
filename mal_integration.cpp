/***************************************************************************
    mal_integration.cpp  - wrapper around the GSL integration functions

                             -------------------
    begin				: 2006-12-04
    copyright			: (C) 2006 by Maurice Leutenegger
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

#include "mal_integration.h"

using namespace std;

Integral::Integral (size_t limit, double epsrel, double epsabs)
  : itsEpsAbs (epsabs), itsEpsRel (epsrel), itsLimit (limit),
    isAllocated (false), itsStatus (0), itsResult (0.), itsAbsErr (0.), 
    itsNEval (0), itsNCalls (0)
{
  F.function = &integrandGSL;
  F.params = (Integral*) this;
  AllocateWorkspace ();
  return;
}

Integral::~Integral ()
{
  FreeWorkspace ();
  return;
}

// Set the absolute error tolerance, don't use relative.
void Integral::setEpsAbs (double epsabs)
{
  itsEpsAbs = epsabs;
  itsEpsRel = 0.;
  return;
}

// Set the relative error tolerance, don't use absolute.
void Integral::setEpsRel (double epsrel)
{
  itsEpsRel = epsrel;
  itsEpsAbs = 0.;
  return;
}

// Set both kinds of error tolerance.
void Integral::setEps (double epsabs, double epsrel)
{
  itsEpsAbs = epsabs;
  itsEpsRel = epsrel;
  return;
}

// Change the size of the workspace. Deallocate the old space and
// reallocate the new one.
void Integral::setLimit (size_t limit)
{
  itsLimit = limit;
  FreeWorkspace ();
  AllocateWorkspace ();
  return;
}

void Integral::AllocateWorkspace ()
{
  if (isAllocated) {
    cout << "Integral: Found unexpected allocated workspace.\n";
    return;
  }
  itsWorkspace = gsl_integration_workspace_alloc (itsLimit);
  isAllocated = true;
  return;
}

void Integral::FreeWorkspace ()
{
  if (isAllocated) {
    gsl_integration_workspace_free (itsWorkspace);
    itsWorkspace = NULL;
    isAllocated = false;
    return;
  }
  cout << "Integral: Found unexpected unallocated workspace.\n";
  return;  
}

// The next several functions are calls to the GSL integration routines.
// Compared to GSL, the details are mostly hidden from the user. 
// Most of the routines require only the integration limits as arguments.
// The routines all return the result directly, rather than the status.
// The result, status, error estimate, and number of calls to the integrand
// can all be accessed after the integration is performed using the
// other member functions.

// Non-adaptive integration on the interval [a,b]
double Integral::qng (double a, double b)
{
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qng 
    (&F, a, b, itsEpsAbs, itsEpsRel, &itsResult, &itsAbsErr, &itsNEval);
  if (itsStatus) {
    handleError ("qng");
    return 0.;
  } 
  return itsResult;
}

// Adaptive integration on the interval [a,b] with user supplied key.
// Key values range from 1 to 6. Default key is 1.
double Integral::qag (double a, double b, int key)
{
  if ((key > 6) || (key < 1)) {
    cerr << "Integral::qag: key " << key << " not valid.\n";
    cerr << "Setting key to 1.\n";
    key = 1;
  }
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qag 
    (&F, a, b, itsEpsAbs, itsEpsRel, itsLimit, key, itsWorkspace, 
     &itsResult, &itsAbsErr);
  if (itsStatus) {
    handleError ("qag");
    return 0.;
  }
  return itsResult;
}

// Adaptive integration with arbitrary singularities on [a,b]
double Integral::qags (double a, double b)
{
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qags
    (&F, a, b, itsEpsAbs, itsEpsRel, itsLimit, itsWorkspace, 
     &itsResult, &itsAbsErr);
  if (itsStatus) {
    handleError ("qags");
    return 0.;
  }
  return itsResult;
}

// Adaptive integration with known singularities at pts[i].
// pts[0] and pts[npts-1] give the endpoints of the integration.
double Integral::qagp (double* pts, size_t npts)
{
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qagp 
    (&F, pts, npts, itsEpsAbs, itsEpsRel, itsLimit, itsWorkspace, 
     &itsResult, &itsAbsErr);
  if (itsStatus) {
    handleError ("qagp (N pts)");
    return 0.;
  }
  return itsResult;
}

// Adaptive integration with only endpoint singularities on [a,b]
double Integral::qagp (double a, double b)
{
  size_t npts = 2;
  double pts[2];
  pts[0] = a;
  pts[1] = b;
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qagp 
    (&F, pts, npts, itsEpsAbs, itsEpsRel, itsLimit, itsWorkspace, 
     &itsResult, &itsAbsErr);
  if (itsStatus) {
    handleError ("qagp (2 points)");
    cout << "a, b: " << a << ", " << b << endl;
    return 0.;
  }
  return itsResult;
}

// Adaptive integration on [-infinity, infinity]
double Integral::qagi ()
{
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qagi
    (&F, itsEpsAbs, itsEpsRel, itsLimit, itsWorkspace, &itsResult, &itsAbsErr);
  if (itsStatus) {
    handleError ("qagi");
    return 0.;
  }
  return itsResult;
}

// Adaptive integration on [a, infinity]
double Integral::qagiu (double a)
{
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qagiu
    (&F, a, itsEpsAbs, itsEpsRel, itsLimit, itsWorkspace, &itsResult, 
     &itsAbsErr);
  if (itsStatus) {
    handleError ("qagiu");
    return 0.;
  }
  return itsResult;
}

// Adaptive integration on [-infinity, b]
double Integral::qagil (double b)
{
  gsl_set_error_handler_off ();
  itsStatus = gsl_integration_qagil
    (&F, b, itsEpsAbs, itsEpsRel, itsLimit, itsWorkspace, &itsResult, 
     &itsAbsErr);
  if (itsStatus) {
    handleError ("qagil");
    return 0.;
  }
  return itsResult;
}

// Static function to call the real integrand.
// Also counts the number of calls.
double Integral::integrandGSL (double x, void* object)
{
  Integral* ThisIntegral = (Integral*) object;
  ThisIntegral->itsNCalls ++;
  return ThisIntegral->integrand (x);
}

void Integral::handleError (string functionName) {
  cout << "GSL error in " << functionName << endl;
  if (itsStatus == GSL_EDIVERGE) {
    cout << "Error code is GSL_EDIVERGE; integral is not converging." 
	 << endl;
  } else if (itsStatus == GSL_EINVAL) {
    cout << "Error code is GSL_EINVAL. Some input was invalid.";
  } else if (itsStatus == GSL_EMAXITER) {
    cout << "Error code is GSL_EMAXITER. Integration exceeded maximum number of iterations.";
  } else if (itsStatus == GSL_EROUND) {
    cout << "Error code is GSL_EROUND. Integration failed because of roundoff error.";
  } else {
    cout << "Error code number is " << itsStatus << endl;
  }
  cout << "Please send debugging information to Maurice." << endl;
  return;
}
