//
// NumericLimits.h
//
// This file contains the NumericLimits template class. 
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London WC1E 7HX, U.K.
// (classlib@mail.cryst.bbk.ac.uk)
// 
// This library is free software; you can redistribute it and/or modify it under the terms of
// the GNU Library General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined(BTL_NUMERICLIMITS_H)
#define BTL_NUMERICLIMITS_H 1

#include <cfloat>
#include <iomanip>

_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description ="This is a template class that can be used to return
     standard floating point constants. N.B. This file will become obsolete when
     the standard C++ file climits becomes more widely available. It contains a
     numeric_limits class which will replace BTL_NumericLimits defined below."]
    [Summary = "template class that returns floating point constants for the
     template parameter"] 
    [Usage="e.g. const float EPSILON = NumericLimits&ltfloat&gt::epsilon(); "]
    [Restrictions="This template class only works for float, double, and long
       double types."]
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/NumericLimits.h>NumericLimits.h</A>"]
    [Friends="None"]
    [Dependencies="None"]
*/

template <class T>
class NumericLimits
{
public:

    	    /**#: [Description="Returns the value of the difference between 1
    	         and the least value greater than 1 that is representable."]*/ 
    static T 
    epsilon() { return DBL_EPSILON; }

    	    /**#: [Description="Returns the minimum finite value."] */ 
    static T 
    min() { return DBL_MIN; }
       
    	    /**#: [Description="Returns the number of base 10 digits which can
    	           be represented without change"] */
    static int 
    digits10() { return DBL_DIG; }
};

class NumericLimits<float>
{
public:
    static float epsilon() { return FLT_EPSILON; }
    static float min() { return FLT_MIN; }
    static int digits10() { return FLT_DIG; }
};

class NumericLimits<double>
{
public:
    static double epsilon() { return DBL_EPSILON; }
    static double min() { return DBL_MIN; }
    static int digits10() { return DBL_DIG; }
};

class NumericLimits<long double>
{
public:
    static long double epsilon() { return LDBL_EPSILON; }
    static long double min() { return LDBL_MIN; }
    static int digits10() { return LDBL_DIG; }
};

_BTL_END_NAMESPACE

#endif
