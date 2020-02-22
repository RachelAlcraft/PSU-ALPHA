//
// btl_vector_algorithms.h
//
// This file contains the generic numeric functions for manipulation of vectors 
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997, 1998 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2004, 2005 University College, Gower Street, London, U.K.
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

#if !defined (BTL_VECTORALGORITHMS_H)
#define BTL_VECTORALGORITHMS_H 1

/**#: [Description ="A collection of generic numerical algorithms for
      vector algebra with an emphasis on the manipulation of vectors in 3D."]
    [Summary = "generic vector algorithms"] 
    [Authors = "D.S.Moss, W.R.Pitt, I.Tickle, M.A.Williams"]
    [Files = "<A HREF=./btl/btl_vector_algorithms.h>btl_vector_algorithms.h</A>"]
    [Dependencies="<A HREF=#numeric_vector>btl_numeric_vector.h</A>,
                   <A HREF=#matrix>btl_matrix.h</A>"]
    [Prerequisites="None<P>"
    <P>
*/

#include <vector>
#include <iterator>

#include "BTL.h"
#include "btl_numeric_vector.h"
#include "btl_matrix.h"

_BTL_BEGIN_NAMESPACE

using namespace std;
    	
//..........................................................................................
// VECTOR ALGEBRA
//..........................................................................................
    	    /**#: [Description="Scalar/dot product of two vectors."]*/
 
	template<class InputIterator1, class InputIterator2, class T >
    	T scalar_product(InputIterator1 firsta, InputIterator1 lasta,
                            InputIterator2 firstb, T init)
	{
            for ( ; firsta != lasta; firsta++, firstb++)
	        init += (*firsta * *firstb);
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Scalar/dot product of two triples ."]*/
 
	template<class InputIterator1, class InputIterator2, class T >
    	T scalar_product(InputIterator1 firsta, InputIterator2 firstb, T init)
	{
	    init += (*firsta * *firstb) + (*(firsta+1) * *(firstb+1)) 
                                           + (*(firsta+2) * *(firstb+2));
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Direct product of two matrices or vectors ."]*/

	template<class InputIterator1, class InputIterator2, class OutputIterator>
    	OutputIterator direct_product(InputIterator1 firsta, InputIterator1 lasta, 
                                      InputIterator2 firstb, OutputIterator result)
	{
	    if(firsta == firstb)
		for (; firsta != lasta; firsta++, result++)
		    *result += *firsta * *firsta;
	    else 
		for (; firsta != lasta; firsta++, firstb++, result++)
		    *result += *firsta * *firstb;
	    return ++result;			 
	}
//..........................................................................................
    	    /**#: [Description="Vector/cross product of two triples."]*/
 
	template<class InputIterator1, class InputIterator2, class OutputIterator>
    	OutputIterator vector_product(InputIterator1 firsta, InputIterator2 firstb, 
                                      OutputIterator result)
	{
	    *result	+= (*(firsta+1) * *(firstb+2)) - (*(firsta+2) * *(firstb+1));
	    *(result+1) += (*(firsta+2) * *firstb)     - (*firsta     * *(firstb+2));
	    *(result+2) += (*firsta     * *(firstb+1)) - (*(firsta+1) * *firstb);
	    return (result+3);
	}

//..........................................................................................
    	    /**#: [Description="Vector product of three triples."]*/
 
	template<class InputIterator1, class InputIterator2, class InputIterator3, 
		 class OutputIterator>
    	OutputIterator triple_vector_product(InputIterator1 firsta, InputIterator2 firstb, 
                                             InputIterator3 firstc, OutputIterator result)
	{
            typename iterator_traits<OutputIterator>::value_type adotc = 0.0;
	    adotc = scalar_product(firsta,firstc,adotc);
            typename iterator_traits<OutputIterator>::value_type adotb = 0.0; 
	    adotb = scalar_product(firsta,firstb,adotb);
	    *result     += (*firstb     * adotc) - (*firstc     * adotb);
	    *(result+1) += (*(firstb+1) * adotc) - (*(firstc+1) * adotb);
	    *(result+2) += (*(firstb+2) * adotc) - (*(firstc+2) * adotb);
	    return (result+3);
	}

//..........................................................................................
    	    /**#: [Description="Scalar product of three triples."]*/
 
	template<class InputIterator1, class InputIterator2, class InputIterator3, class T>
    	T triple_scalar_product(InputIterator1 firsta, InputIterator2 firstb, 
                                  InputIterator3 firstc, T init)
	{
	    init += (*firsta * ((*(firstb+1) * *(firstc+2))-(*(firstb+2)* *(firstc+1))))
		  + (*(firsta+1) * ((*(firstb+2) * *firstc)-(*firstb * *(firstc+2))))
		  + (*(firsta+2) * ((*firstb * *(firstc+1))-(*(firstb+1)* *firstc)));
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Calculate the squared distance between the points 
                                 represented by two vectors."]*/
 
	template<class InputIterator1, class InputIterator2, class T>
    	T separation_squared(InputIterator1 firsta, InputIterator1 lasta,
                            InputIterator2 firstb, T init)
	{
            for ( ; firsta != lasta; firsta++, firstb++)
            {
	        init += (*firsta - *firstb) * (*firsta - *firstb);
	    }	
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Calculate the squared distance between points represented 
                                by two triples."]*/
 
	template<class InputIterator1, class InputIterator2, class T>
    	T separation_squared(InputIterator1 firsta, InputIterator2 firstb, T init)
	{
	    init  += ((*firsta - *firstb) * (*firsta - *firstb)) 
	           + ((*(firsta+1) - *(firstb+1)) * (*(firsta+1) - *(firstb+1)))
	           + ((*(firsta+2) - *(firstb+2)) * (*(firsta+2) - *(firstb+2)));
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Calculate the distance between the points represented 
                                by two vectors."]*/
 
	template<class InputIterator1, class InputIterator2, class T>
    	T separation(InputIterator1 firsta, InputIterator1 lasta,
                            InputIterator2 firstb, T init)
	{
	    return T(sqrt(separation_squared(firsta,lasta,firstb,init)));
	}

//..........................................................................................
    	    /**#: [Description="Calculate the distance between the points represented 
                                by two triples."]*/
 
	template<class InputIterator1, class InputIterator2, class T>
    	T separation(InputIterator1 firsta, InputIterator2 firstb, T init)
	{
	    return T(sqrt(separation_squared(firsta,firstb,init)));
	}


//..........................................................................................
    	    /**#: [Description="Returns the sum of the vector elements."]*/
  
	template<class InputIterator, class T>
    	T sum(InputIterator first, InputIterator last, T init)
	{
            for ( ; first != last; first++)
	        init += *first;
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Returns the sum of the elements of a triple."]*/
 
 	template<class InputIterator, class T>
    	T sum(InputIterator firsta, T init)
	{
	    return init += *firsta + *(firsta+1) + *(firsta+2);
	}

//..........................................................................................
    	    /**#: [Description="Returns the sum of the vector elements, 
                   but accumulates the elements in order 
                   of ascending value in order to minimise rounding errors."]*/
 
	template<class InputIterator, class T>
    	T sum_precise(InputIterator first, InputIterator last, T init)
	{
            typename vector<T>::const_iterator i;
	    vector<T> v(first,last);
	    sort(v.begin(), v.end());
	    for (i=v.begin(); i!=v.end(); i++)
	    	init += *i;
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Returns the sum of the squares of the vector elements."]*/
  
	template<class InputIterator, class T>
    	T sum_of_squares(InputIterator first, InputIterator last, T init)
	{
            for ( ; first != last; first++)
	        init += *first * *first;
	    return init;
	}

//..........................................................................................
    	    /**#: [Description="Returns the sum of the squares of elements of a triple."]*/
 
 	template<class InputIterator, class T>
    	T sum_of_squares(InputIterator firsta, T init)
	{
	    return init += (*firsta * *firsta) + (*(firsta+1) * *(firsta+1)) 
			                        + (*(firsta+2) * *(firsta+2));
	}

//..........................................................................................
    	    /**#: [Description="Returns the sum of the squares of the vector elements, 
                   but accumulates the elements in order 
                   of ascending value in order to minimise rounding errors."]*/
 
	template<class InputIterator, class T>
    	T sum_of_squares_precise(InputIterator first, InputIterator last, T init)
	{
            typename vector<T>::const_iterator i;
	    vector<T> v(first,last);
	    sort(v.begin(), v.end());
	    for (i=v.begin(); i!=v.end(); i++)
	    	init += *i * *i;
	    return init;
	}
//..........................................................................................
    	    /**#: [Description="Returns the length of a vector."]*/
  
	template<class InputIterator1, class T>
    	T magnitude(InputIterator1 firsta, InputIterator1 lasta, T init )
	{
	    return sqrt(sum_of_squares(firsta,lasta,init));
	}

//..........................................................................................
    	    /**#: [Description="Returns the precise length of a vector."]*/
  
	template<class InputIterator1, class T>
    	T magnitude_precise(InputIterator1 firsta, InputIterator1 lasta, T init )
	{
	    return sqrt(sum_of_squares_precise(firsta,lasta,init));
	}

//..........................................................................................
    	    /**#: [Description="Returns the length of a triple."]*/
 
 	template<class InputIterator1, class T>
    	T magnitude(InputIterator1 firsta, T init)
	{
	    return sqrt(sum_of_squares(firsta,init));
	}

//} // end of namespace

//.............................................................................

    	    /**#: [Description="Rotates each triple in a container using a given 
	                        rotation matrix and origin"]*/

template <class ForwardIterator, class ConstForwardIterator1, class ConstForwardIterator2>
void rotate(ForwardIterator begin, ForwardIterator end, 
	    ConstForwardIterator1 rotmatrix, ConstForwardIterator2 origin)
{
	for(; begin!=end; begin+=3)
	{
	    *begin     -= *origin;
 	    *(begin+1) -= *(origin+1);
	    *(begin+2) -= *(origin+2);
	    BTL_REAL temp_x = *begin;
           BTL_REAL temp_y = *(begin+1);

           *begin     = (*begin * *rotmatrix)     + (*(begin+1) * *(rotmatrix+1))
                      + (*(begin+2) * *(rotmatrix+2));

           *(begin+1) = (temp_x * *(rotmatrix+3)) + (temp_y * *(rotmatrix+4))
                      + (*(begin+2) * *(rotmatrix+5));

           *(begin+2) = (temp_x * *(rotmatrix+6)) + (temp_y * *(rotmatrix+7))
                      + (*(begin+2) * *(rotmatrix+8));
	    *begin     += *origin;
 	    *(begin+1) += *(origin+1);
	    *(begin+2) += *(origin+2);
	}
        
}                    				       
//.............................................................................
    	    /**#: [Description="Rotates each triple in a container about a given 
	                        axis and origin by a given number of radians."] */

template <class ForwardIterator, class ConstForwardIterator1, class ConstForwardIterator2, class T>
void rotate(ForwardIterator first, ForwardIterator last, 
	    ConstForwardIterator1 axis, ConstForwardIterator2 origin, T angle)
	{
	// Constructor for 3x3 rotation matrix for a rotation about a given axis by a
	// given number of radians. Create the following matrix:

	//  [cosA+a1x.a1x(1-cosA)      -a3xsinA+a1a2(1-cosA)	a2xsinA+a1xa3x(1-cosA) ] 
	//  [a3xsinA+a1xa2(1-cosA)    cosA+a2xa2x(1-cosA)   	-a1xsinA+a2xa3x(1-cosA)] 
	//  [-a2xsinA+a3xa1x(1-cosA)   a1xsinA+(1-cosA)a2xa3x 	cosA+(1-cosA)a3xa3x   ]

	// where (a1x,a2x,a3x) is the unit vector in the direction of the axis and A is the angle
	// of rotation in radians.

            BTL_REAL recipmagofaxis = 1.0/sqrt(((*axis * *axis)+(*(axis+1) * *(axis+1))+(*(axis+2) * *(axis+2)))) ;
	    BTL_REAL a1x     = *axis * recipmagofaxis;
	    BTL_REAL a2x     = *(axis+1) * recipmagofaxis;
	    BTL_REAL a3x     = *(axis+2) * recipmagofaxis;
	    BTL_REAL cosA    = cos(angle);
	    BTL_REAL sinA    = sin(angle);
	    BTL_REAL a1mcosA = a1x * (1.0-cosA);    	
	    BTL_REAL a2mcosA = a2x * (1.0-cosA);    	
	    BTL_REAL a3mcosA = a3x * (1.0-cosA);    	
	    BTL_REAL a1sinA  = a1x * sinA;    	
	    BTL_REAL a2sinA  = a2x * sinA;    	
	    BTL_REAL a3sinA  = a3x * sinA;

            vector<BTL_REAL> rotation;
	    rotation.push_back(( cosA    + a1x * a1mcosA));
            rotation.push_back((-a3sinA  + a1x * a2mcosA));
            rotation.push_back(( a2sinA  + a1x * a3mcosA));
            rotation.push_back(( a3sinA  + a1x * a2mcosA));
            rotation.push_back(( cosA    + a2x * a2mcosA));
            rotation.push_back((-a1sinA  + a2x * a3mcosA));
            rotation.push_back((-a2sinA  + a1x * a3mcosA));
            rotation.push_back(( a1sinA  + a2x * a3mcosA));
            rotation.push_back(( cosA    + a3x * a3mcosA));

            rotate(first,last,rotation.begin(),origin);
	}

//.............................................................................
    	    /**#: [Description="Translates a each triple in a container 
	                        by a given vector"]*/

template <class ForwardIterator, class ConstForwardIterator>
void translate(ForwardIterator begin, ForwardIterator end, ConstForwardIterator trans)
{
	for(; begin!=end; begin+=3)
	{
	    *begin     += *trans;
 	    *(begin+1) += *(trans+1);
	    *(begin+2) += *(trans+2);
	}	
}


_BTL_END_NAMESPACE

//    const double DEGRAD = M_PI/180.0; 	// converts degrees to radians
//    double radAngle = angle * DEGRAD;

#endif

