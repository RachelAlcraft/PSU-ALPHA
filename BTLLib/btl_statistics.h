//
// btl_statistics.h
//
// This file contains generic statistical functions 
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997-1999 Birkbeck College, Malet Street, London, U.K.
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Library General Public License as published 
// by the Free Software Foundation; either version 2 of the License, or 
// (at your option) any later version.  This library is distributed in the
// hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined (BTL_STATISTICS_H)
#define BTL_STATISTICS_H 1

/**#: [Description ="A collection of generic statistical functions."]
    [Summary = "generic statistical algorithms"] 
    [Authors = "M.A.Williams"]
    [Files = "<A HREF=./btl/btl_statistics.h>btl_statistics.h</A>"]
    [Dependencies="None"]
    [Prerequisites="None<P>"
    <P>
*/

   
#include <cmath>
#include <vector>

#include "BTL.h"
#include "btl_numerics.h"

    	
_BTL_BEGIN_NAMESPACE

using namespace std;

//..........................................................................................
// STATISTICS
//..........................................................................................
	    	/**#: [Description="Calculate the mean value of the elements in a 
		       	container."] */

	template<class InputIterator, class T>
    	T mean(InputIterator first, InputIterator last, T mean)
	{
            typename vector<T>::const_iterator i;
	    vector<T> v(first,last);
	    sort(v.begin(), v.end());
	    for (i=v.begin(); i!=v.end(); i++)
	    	mean += *i;
	    return (mean/(last-first));
	}

//..........................................................................................
	    	/**#: [Description="Calculate mean absolute deviation from the mean of the
			elements in a container."] */

	template<class InputIterator, class T>
    	T mean_absolute_deviation(InputIterator first, InputIterator last, T mad, T& mean)
	{
            typename vector<T>::iterator i;
	    vector<T> v(first,last);
	    for (i=v.begin(); i!=v.end(); i++)
	    	*i = fabs(*i - mean);
	    sort(v.begin(), v.end());
	    for (i=v.begin(); i!=v.end(); i++)
		mad += *i;
	    return (mad/(last-first));
	}

//..........................................................................................
	    	/**#: [Description="Calculate the sample variance of the elements
		       in a container."] */

	template<class InputIterator, class T>
    	T variance(InputIterator first, InputIterator last, T var, T& mean)
	{
            typename vector<T>::iterator i;
	    vector<T> v(first,last);
	    for (i=v.begin(); i!=v.end(); i++)
	    	*i = (*i - mean)*(*i -mean);
	    sort(v.begin(), v.end());
	    for (i=v.begin(); i!=v.end(); i++)
		var += *i;
	    return (var/(last-first-1));
	}

//..........................................................................................
	    	/**#: [Description="Calculate the mean, mean absolute deviation from mean,
			sample variance, skew and kurtosis of the distribution
			of the elements in a container."] */

	template<class InputIterator, class T>
    	void normal_statistics(InputIterator first, InputIterator last, 
			       T& mean, T& mad, T& var, T& skew, T& kurtosis)
	{
	    BTL_REAL reciprocal_n = (1.0/(last-first));
	    BTL_REAL reciprocal_nminusone = (1.0/(last-first-1));
            typename vector<T>::iterator i;
	    vector<T> v(first,last);
	    sort(v.begin(), v.end());

	    for (i=v.begin(); i!=v.end(); i++)
	    	mean += *i;
	    mean *= reciprocal_n;

	    for (i=v.begin(); i!=v.end(); i++)
	    	*i = (*i - mean);
	    sort(v.begin(), v.end(), less_absolute<T>());
	    for (i=v.begin(); i!=v.end(); i++)
		mad += fabs(*i);
	    mad *= reciprocal_n;

	    for (i=v.begin(); i!=v.end(); i++)
		var += *i * *i;
	    var *= reciprocal_nminusone;
	    BTL_REAL reciprocal_standard_dev = (1.0 / sqrt(var));

	    for (i=v.begin(); i!=v.end(); i++)
	    	*i = *i * reciprocal_standard_dev;
	    for (i=v.begin(); i!=v.end(); i++)
		skew += *i * *i * *i;
	    skew *= reciprocal_n;

	    for (i=v.begin(); i!=v.end(); i++)
		*i = *i * *i;
	    for (i=v.begin(); i!=v.end(); i++)
	    	kurtosis += *i * *i;
	    kurtosis *= reciprocal_n;
	}

// NUMERICALLY UNSTABLE, ADD TRAPS FOR UNDERFLOW


_BTL_END_NAMESPACE

#endif

