/***************************************************************************
    mal_integration.h   - wrapper around the GSL integration functions

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

#ifndef MAL_INTEGRATION_H
#define MAL_INTEGRATION_H

#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace std;

/*-------------------------------------------------------------------------*\
  This is an abstract class for integration. To do an integral, declare a 
  derived class with a member function 
  double DerivedClass::integrand (double x).
  Then declare an instance of DerivedClass and call one of the qxxx 
  member functions to do the integral. 
  
  Example:
  
  class MyIntegral : public Integral
  {
   public:
    MyIntegral (...);
    ...
    double integrand (double x);
   private:
    double some_parameter;
    ...
  }

  main ()
  {
    double a = ...;
    double b = ...;
    MyIntegral I (...);
    I.setEpsRel (1.e-5);
    double answer = I.qag (a, b);
    double error = I.getAbsErr ();
    size_t NCalls = I.getNCalls ();
    return;
  }
\*-------------------------------------------------------------------------*/

class Integral
{
 public:
  // Constructor allocates workspace
  Integral (size_t limit = 1000, double epsrel = 1.e-4, double epsabs = 0.);
  // Destructor deallocates memory for workspace.
  virtual ~Integral ();
  // The copy constructor and assignment operator have been defined
  // empty in the private section so they won't be used.
  // The GSL integration routines. See GSL documentation for details.
  double qng (double a, double b); // non-adaptive
  double qag (double a, double b, int key = 1); // adaptive
  double qags (double a, double b); // adaptive with singularities
  double qagp (double* pts, size_t npts); // adaptive with known singularities
  double qagp (double a, double b); // only endpoint singularities
  double qagi (); // adaptive from -infinity to infinity
  double qagiu (double a); // adaptive from a to infinity
  double qagil (double b); // adaptive from -infinity to b
  // If you want to change the accuracy goal of the integration. 
  //   Setting either parameter individually zeroes out the other.
  void setEpsAbs (double epsabs); 
  void setEpsRel (double epsrel); // default 1.e-4
  void setEps (double epsabs, double epsrel);
  // Changes the size of the workspace; deallocates and reallocates memory.
  void setLimit (size_t limit); // default 1000
  size_t getLimit () const {return itsLimit;}
  // This is the fake integrand. Each call increments itsNCalls
  static double integrandGSL (double x, void* object);
  void resetNCalls () {itsNCalls = 0; return;}
  // This is the real integrand. It is a pure virtual function which is 
  // unimplemented and must be overridden in a derived class. 
  virtual double integrand (double x) = 0;
  // If you want to get out information on the integration after the fact.
  int getStatus () const {return itsStatus;}
  double getResult () const {return itsResult;}
  double getAbsErr () const {return itsAbsErr;}
  size_t getNEval () const {return itsNEval;} // for qng only
  size_t getNCalls () const {return itsNCalls;}
 private:
  // settings
  double itsEpsAbs;
  double itsEpsRel;
  size_t itsLimit;
  // gsl integration objects
  gsl_function F;
  gsl_integration_workspace* itsWorkspace;
  bool isAllocated;
  // results
  int itsStatus;
  double itsResult;
  double itsAbsErr;
  size_t itsNEval; // for qng
  size_t itsNCalls; // for others
  void AllocateWorkspace ();
  void FreeWorkspace ();
  void handleError (string functionName);
  Integral (const Integral& I); // no copy constructor
  //  Integral operator = (const Integral& I); //no assignment operator
};

#endif
