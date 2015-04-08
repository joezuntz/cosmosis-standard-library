// Module for X-ray cluster fgas data
// Copyright (c) 2013 Adam Mantz, amantz@slac.stanford.edu

// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#ifndef _FWRAPPER_
#define _FWRAPPER_

#include "wrapper.hpp"

namespace Clusters {
  Wrapper *theWrapper;
  double cangluardistance2(const double, const double);
}

extern "C" {

  typedef double fortran_dl;
  typedef float fortran_real;
  typedef long fortran_integer;
  typedef bool fortran_logical1;
  typedef char fortran_character;

  fortran_integer fclwrapperdatasetsize_(fortran_integer *i);
  fortran_dl fclwrappergetredshift_(fortran_integer *i, fortran_integer *j);
  fortran_integer fclwrapperinit_(fortran_character *a_char, fortran_integer a_len);
  fortran_dl fclwrapperlnp_(fortran_dl *datapars);
  fortran_integer fclwrapperloaddata_(fortran_character *a_char, fortran_integer a_len);
  fortran_integer fclwrapperloadparameters_(fortran_dl *modelpars);
  fortran_integer fclwrappernumdatapars_();
  fortran_integer fclwrappernumparameters_();
  fortran_integer fclwrappernumclusters_();
  fortran_integer fclwrappernumconst_();
  void fclwrappersetconst_(fortran_dl *physics);
  void fclwrappersetdataoptions_(fortran_integer *i, fortran_integer *j);
  void fclwrappersetmasspivot_(fortran_dl *M);
  void fclwrappersimulate_(fortran_dl *datapars, fortran_integer *iscatter, fortran_integer *mscatter);

  fortran_dl fangulardistance2_(fortran_dl*, fortran_dl*); // provided by cosmomc
  
}

#endif
