/***************************************************************************
    RootFinderNewton.cpp - Wrapper for the gsl root finder using Newton's
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

#include <iostream>
#include "mal_RootFinderNewton.h"

using namespace std;

RootFinderNewton::RootFinderNewton 
(double guess, double RelativeError, double AbsoluteError) 
  : itsStatus (0), itsIteration (0), itsRelativeError (RelativeError),
    itsAbsoluteError (AbsoluteError), itsRoot (guess),
    T (gsl_root_fdfsolver_newton), s (0)
{
  FDF.f = &fGSL;
  FDF.df = &dfGSL;
  FDF.fdf = &fdfGSL;
  FDF.params = (RootFinderNewton*) this;
  initializeRootFinderNewton ();
  return;
}

RootFinderNewton::~RootFinderNewton ()
{
  freeFDFSolver ();
  return;
}

void RootFinderNewton::initializeRootFinderNewton ()
{
  s = gsl_root_fdfsolver_alloc (T);
  return;
}

void RootFinderNewton::freeFDFSolver ()
{
  gsl_root_fdfsolver_free (s);
  return;
}

double RootFinderNewton::fGSL (double x, void* object)
{
  RootFinderNewton* ThisRootFinderNewton = (RootFinderNewton*) object;
  return ThisRootFinderNewton->f (x);
}

double RootFinderNewton::dfGSL (double x, void* object)
{
  RootFinderNewton* ThisRootFinderNewton = (RootFinderNewton*) object;
  return ThisRootFinderNewton->df (x);
}

void RootFinderNewton::fdfGSL (double x, void* object, double *y, double* dy)
{
  double function;
  double derivative;
  RootFinderNewton* ThisRootFinderNewton = (RootFinderNewton*) object;
  ThisRootFinderNewton->fdf (x, function, derivative);
  *y = function;
  *dy = derivative;
  return;
}

double RootFinderNewton::findRoot ()
{
  gsl_root_fdfsolver_set (s, &FDF, itsRoot);
  itsIteration = 0;
  double LastRoot = itsRoot;
  do {
    ++itsIteration;
    itsStatus = gsl_root_fdfsolver_iterate (s);
    LastRoot = itsRoot;
    itsRoot = gsl_root_fdfsolver_root (s);
    itsStatus = gsl_root_test_delta 
      (itsRoot, LastRoot, itsAbsoluteError, itsRelativeError);
  } while (itsStatus == GSL_CONTINUE && itsIteration < MAXIMUM_ITERATIONS);
  if (itsIteration == MAXIMUM_ITERATIONS) {
    cout << "RootFinderNewton::findRoot (): maximum iterations reached.\n";
  }
  return itsRoot;
}
