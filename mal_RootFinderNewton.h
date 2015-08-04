/***************************************************************************
    RootFinderNewton.h - Wrapper for the gsl root finder using Newton's
                         method. Derived class must have functions 
                         f, df, and fdf.

                             -------------------
    begin				: March 2007
    copyright			: (C) 2007 by Maurice Leutenegger
    email				: maurice@astro.columbia.edu
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

#ifndef MAL_ROOT_FINDER_NEWTON_H
#define MAL_ROOT_FINDER_NEWTON_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

class RootFinderNewton {
 public:
  RootFinderNewton (double guess = 0., double RelativeError = 1.e-3, 
		    double AbsoluteError = 0.);
  virtual ~RootFinderNewton ();
  double findRoot ();
  double getRoot () {return itsRoot;}
  int getIterations () const {return itsIteration;}
  void setGuess (double guess) {itsRoot = guess; return;}
  static double fGSL (double x, void* params);
  static double dfGSL (double x, void* params);
  static void fdfGSL (double x, void* params, double* y, double* dy);
  virtual double f (double x) = 0;
  virtual double df (double x) = 0;
  virtual void fdf (double x, double& y, double& dy) = 0;
 private:
  const static int MAXIMUM_ITERATIONS = 100;
  int itsStatus;
  int itsIteration;
  double itsRelativeError;
  double itsAbsoluteError;
  double itsRoot;
  const gsl_root_fdfsolver_type* T;
  gsl_root_fdfsolver* s;
  gsl_function_fdf FDF;
  void initializeRootFinderNewton ();
  void freeFDFSolver ();
};


#endif//MAL_ROOT_FINDER_NEWTON_H
