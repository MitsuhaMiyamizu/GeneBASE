////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

/// @file   stats-distributions.h
/// @brief  stats functions taken from perls Statistics::Distributions

//! We want the output to be the same as the perl code, so
//! we have converted the code to C. Refer to the perldocs
//! for more info.

//! These used to live in stats-distributions.cpp, but are now templates
//! so we get the float and double versions.

#ifndef __STATS_DISTRIBUTIONS_H_
#define __STATS_DISTRIBUTIONS_H_

// avoids indenting all the functions.
namespace affxstat {
  template <typename T1> T1 chisqrprob(int n,T1 x);
  template <typename T1> T1 uprob(T1 x);
}

#ifdef WIN32
#define M_PI       3.14159265358979323846
#endif

#include <math.h>

/// @brief     The chi-squared probability
/// @param     n         degrees of freedom
/// @param     x         chi-square value
/// @return    upper probability of the chi-square
template <typename T1>
T1
affxstat::chisqrprob(int n,T1 x)
{
  double p;

  if (x<=0) {
    p=1;
  } 
  else if (n>100) {
    p=uprob((pow((double)(x/n),(double)(1.0/3.0))
             -(1.0-2.0/9.0/n))/sqrt(((2.0/9.0)/n)));
  } 
  else if (x>400) {
    p=0;
  }
  else {  
    double a, i, i1;
    if ((n%2) != 0) {
      p=2.0*uprob(sqrt(x));
      a=sqrt((float)(2/M_PI))*exp(-x/2)/sqrt(x);
      i1=1;
    } 
    else {
      p=a=exp(-x/2.0);
      i1=2;
    }
    for (i=i1;i<=(n-2);i+=2) {
      a*=x/i;
      p+=a;
    }
  }
  return (T1)p;
}

/// @brief     upper probability of the u distribution (u=-0.85)
/// @param     x         ???
/// @return    0.0-1.0
template <typename T1>
T1
affxstat::uprob(T1 x)
{
  double p;
  double absx;

  p=0;
  absx=fabs(x);

  if (absx<1.9) {
    p=pow((1.0 +
           absx*(0.049867347
                 +absx*(0.0211410061
                        +absx*(0.0032776263
                               +absx*(0.0000380036
                                      +absx*(0.0000488906
                                             +absx*0.000005383)))))),-16)/2; // ?
  } 
  else if (absx<=100) {
    for (int i=18; i>=1;i--) {
      p=i/(absx+p);
    }
    p=exp(-0.5*absx*absx)
      / sqrt(2.0*M_PI)/(absx + p);
  }

  if (x<0) {
    p=1-p;
  }
  return (T1)p;
}

//}

#endif // __STATS_DISTRIBUTIONS_H_
