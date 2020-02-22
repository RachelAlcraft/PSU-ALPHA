//
// btl_least_squares.h
//
// This file contains the LeastSquaresFit template class. 
//
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2004,2005 University College London, Gower Street, London, U.K.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as 
// published by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.  This library is distributed in the
// hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined (BTL_LEASTSQUARESFIT_H)
#define BTL_LEASTSQUARESFIT_H 1

#include <cmath>
#include <limits>

#include "BTL.h"
#include "btl_matrix.h"
#include "btl_numeric_vector.h"
#include "btl_vector_algorithms.h"
#include "btl_matrix_algorithms.h"

_BTL_BEGIN_NAMESPACE

using std::sqrt;
using std::numeric_limits;

/**#: [Description="This file contains algorithms for calculating the applying
       the least squares fit of one set of x,y,z coordinates onto another. The
       method used is that of Simon K. Kearsley (Acta Cryst.'89, A45, 208-10).
       "]
      [Summary = "least squares superposition of two sets of 3D coordinates"] 
      [Usage="All members in this class are static and therefore an object
       need not be create before the members can be used. The methods in this
       class are generic algorithms and can be used with a wide range of
       container types. This flexibility is achieved through the use of
       iterators. Any container used must have an appropriate iterator type
       that can be used to traverse through its elements. Most methods require
       two iterators per container, marking the first and last element to be
       processed. This style of function is designed to be the same as used in
       the Standard Template Library.<P>
       e.g.<P>
            <CODE>
               // Construct two empty vectors<P>
               vector v1,v2;<P>
               // Put some coordinates into these vector<P>
               .... <P>
               // Do least squares fit on these coodinates<P>
               BTL_REAL rmsd = 0.0;
               rmsd = lsqfit(v1.begin(),v1.end(),v2.begin(),v2.end(),rmsd);<P>
            </CODE>
       "]
      [Files="<A HREF=./btl/btl_least_squares.h>btl_least_squares.h</A>"]
      [Dependencies="<A HREF=#matrix<>>btl_matrix.h<></A>"]
      [Authors="D.Moss, W.R.Pitt, M.A.Williams"]*/


//.............................................................................
	    	/**#: [Description="Calculate the root mean squared deviation
		       after a least squares fit of the two sets of coordinates
	    	       using the method of Kearsley(Acta Cryst. '89, A45, 208-10).
		       Both sets of coordinates remain unchanged by this function"]
		      [Restrictions="Each set must have the same number of
		       coordinates. Each coordinate must be orthogonal and in
		       three dimensions (i.e. normal x,y,z coordinates).
		       WARNING: if coordinates with varying dimensions are input
		       then this may well not be detected."]*/

template <class ConstForwardIterator1, class ConstForwardIterator2, class T>
T lsqfit_rmsd(ConstForwardIterator1 begin1, ConstForwardIterator1 end1,
	      ConstForwardIterator2 begin2, ConstForwardIterator2 end2,
	      T init)			        
{    
    // Calculate Kearsley matrix
    matrix<> matfit(4,4);
    _kearsley_matrix(begin1,end1,begin2,end2,matfit.begin());
    
    // Calculate eigenvectors and eigenvalues of matfit
    matrix<> evec(4,4);
    numeric_vector<> eval(4);
    eigen_solution(matfit.begin(),matfit.end(),4,evec.begin(),eval.begin());

    // Calculate RMSD from smallest eigen value. If this value is very small
    // then round rmsd down to zero.

    if (eval(4) < (eval(1) * numeric_limits<BTL_REAL>::epsilon())) 
    	init += 0.0;
    else 
    	init += sqrt(3*eval(4)/(end1-begin1));

    return init;
}

//.............................................................................
	    /**#: [Description="Least squares fit one set of coordinates to
	    	   another using the method of Kearsley. Return the root
	    	   mean squared deviation after the fit. The second set of
	    	   coordinates remains unaltered by this function."]
		   [Restrictions="Each set must have the same number of
		   coordinates. Each coordinate must be orthogonal and in
		   three dimensions (i.e. normal x,y,z coordinates).
		   WARNING: if coordinates with varying dimensions are input
		   then this may well not be detected."]*/


template <class ForwardIterator, class ConstForwardIterator, class T>
T lsqfit(ForwardIterator begin1, ForwardIterator end1,
	 ConstForwardIterator begin2, ConstForwardIterator end2, T init)
{   
    BTL_INT size = check_3d_data(begin1, end1, begin2, end2);

   // Calculate centroids of two sets of coordinates
    numeric_vector<> centroid1,centroid2,translation; 
    centroid(begin1,end1,centroid1.begin());
    centroid(begin2,end2,centroid2.begin());
    translation = centroid2 - centroid1;

    // Calculate Kearsley matrix
    matrix<> matfit(4,4);
    _kearsley_matrix(begin1,end1,begin2,end2,matfit.begin());
    
    // Calculate eigenvectors and eigenvalues of matfit
    matrix<> evec(4,4);
    numeric_vector<> eval(4);
    eigen_solution(matfit.begin(),matfit.end(),4,evec.begin(),eval.begin());
    transpose(evec.begin(), evec.end(), 4, evec.begin());

    // Create the rotation matrix. O(N)
    matrix<> rot(3,3);
    rotation_from_fit(evec.begin(),rot.begin());

    // Fit atoms of inputMatrix using rotation and translations. O(N)
    rotate(begin1,end1,rot.begin(),centroid1.begin());
    translate(begin1,end1,translation.begin());

    // Calculate RMSD from smallest eigen value. 
    // If this value is very small then round rmsd down to zero.

    if (eval(4) < (eval(1) * numeric_limits<BTL_REAL>::epsilon())) 
	init += 0.0;
    else 
	init += sqrt(eval(4)/size);

    return init;
}  

//.............................................................................
// Calculate the centroid of a set of x,y,z coordinates held in an n x 3 matrix.

template <class ConstForwardIterator, class OutputIterator>
OutputIterator centroid(ConstForwardIterator first, ConstForwardIterator last,
                        OutputIterator result)
{
    BTL_REAL reciprocal_n = 3.0 / (last - first);  

    for(; first != last; first+=3)
    {
        *result += *first;
        *(result+1) += *(first+1);
        *(result+2) += *(first+2);
    }

    *result *= reciprocal_n;
    *(result+1) *= reciprocal_n;
    *(result+2) *= reciprocal_n;

    return (result+3);
}

		
//.............................................................................

template <class ForwardIterator, class OutputIterator>
OutputIterator _kearsley_matrix(ForwardIterator begin1, ForwardIterator end1,
		                ForwardIterator begin2, ForwardIterator end2,
			        OutputIterator matfit)
{ 
    // Calculate centroids of two sets of coordinates
    numeric_vector<> centroid1,centroid2; 
    centroid(begin1,end1,centroid1.begin());
    centroid(begin2,end2,centroid2.begin());

    // Set up matrix to be diagonalized

    numeric_vector<> summ,diff;
    BTL_REAL sumcentres1 = centroid1[0] + centroid2[0];
    BTL_REAL sumcentres2 = centroid1[1] + centroid2[1];
    BTL_REAL sumcentres3 = centroid1[2] + centroid2[2];
    BTL_REAL diffcentres1 = centroid1[0] - centroid2[0];
    BTL_REAL diffcentres2 = centroid1[1] - centroid2[1];
    BTL_REAL diffcentres3 = centroid1[2] - centroid2[2];

    for ( ; begin1 != end1; begin1+=3,begin2+=3)
    {	
    	summ[0] = *begin2     +  *begin1     - sumcentres1;
    	summ[1] = *(begin2+1) +  *(begin1+1) - sumcentres2;
    	summ[2] = *(begin2+2) +  *(begin1+2) - sumcentres3;
    	diff[0] = *begin2     -  *begin1     + diffcentres1;
    	diff[1] = *(begin2+1) -  *(begin1+1) + diffcentres2;
    	diff[2] = *(begin2+2) -  *(begin1+2) + diffcentres3;
	
	*matfit      += diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
	*(matfit+1)  += summ[1]*diff[2] - diff[1]*summ[2];
	*(matfit+2)  += diff[0]*summ[2] - summ[0]*diff[2];
	*(matfit+3)  += summ[0]*diff[1] - diff[0]*summ[1];

	*(matfit+5)  += summ[1]*summ[1] + summ[2]*summ[2] + diff[0]*diff[0];
	*(matfit+6)  += diff[0]*diff[1] - summ[0]*summ[1];
	*(matfit+7)  += diff[0]*diff[2] - summ[0]*summ[2];

	*(matfit+10) += summ[0]*summ[0] + summ[2]*summ[2] + diff[1]*diff[1];
	*(matfit+11) += diff[1]*diff[2] - summ[1]*summ[2];

	*(matfit+15) += summ[0]*summ[0] + summ[1]*summ[1] + diff[2]*diff[2];
    }

	*(matfit+4) = *(matfit+1);
	*(matfit+8) = *(matfit+2);
	*(matfit+9) = *(matfit+6);
	*(matfit+12) = *(matfit+3);
	*(matfit+13) = *(matfit+7);
	*(matfit+14) = *(matfit+11);
	    
    return (matfit+16);
}

//.............................................................................
// Create fitting rotation matrix from the eigen vector that corresponds to the
// smallest eigen vector of Kearsley's matrix.

template <class ForwardIterator, class OutputIterator>
OutputIterator rotation_from_fit(ForwardIterator eigenmatrix, OutputIterator rot)
{
    BTL_REAL q1=*(eigenmatrix+3);
    BTL_REAL q2=*(eigenmatrix+7);
    BTL_REAL q3=*(eigenmatrix+11);
    BTL_REAL q4=*(eigenmatrix+15);

    // initialize rotation matrix
    *rot = q1*q1 + q2*q2 - q3*q3 - q4*q4;
    *(rot+1) = 2.0 * (q2*q3 + q1*q4);
    *(rot+2) = 2.0 * (q2*q4 - q1*q3);

    *(rot+3) = 2.0 * (q2*q3 - q1*q4);
    *(rot+4) = q1*q1 + q3*q3 - q2*q2 - q4*q4;
    *(rot+5) = 2.0 * (q3*q4 + q1*q2);

    *(rot+6) = 2.0 * (q2*q4 + q1*q3);
    *(rot+7) = 2.0 * (q3*q4 - q1*q2);
    *(rot+8) = q1*q1 + q4*q4 - q2*q2 - q3*q3;
    
    return (rot+9);
}                                      

    	     	
//.............................................................................
	    	/**#: [Description="Calculate the root mean squared deviation
		       of two sets of coordinates."] 
		      [Restrictions="Each set must have the same number of
		       coordinates. There is no restriction on the number of
		       dimensions which the coordinates represent (except that
		       it must be the same in every case). WARNING: if
		       coordinates with varying dimensions are input then this
		       may well not be detected."]*/

template <class ConstForwardIterator, class T>
T rmsd(ConstForwardIterator begin1,ConstForwardIterator end1,
	      ConstForwardIterator begin2, T init)
{    
    init = separation_squared(begin1,end1,begin2,init);    
    return sqrt((3*init)/(end1-begin1));
}


//.............................................................................
// Check that both input containers have the same number of elements and that
// number is multiple of 3. Output the number of x,y,z positions

template <class ForwardIterator1, class ForwardIterator2>
BTL_INT check_3d_data(ForwardIterator1 begin1, ForwardIterator1 end1,
 			    ForwardIterator2 begin2, ForwardIterator2 end2)
{
    BTL_INT ndata1=0;
    while (begin1 != end1)
    {
    	ndata1++;
    	begin1++;
    }

    if (ndata1 == 0)
	FATAL_ERROR("There are no elements in 1st container");

    if (ndata1%3 != 0)
    	FATAL_ERROR("The number of elements in the 1st container is not a "   
    	    	    "multiple of 3");
    	    	    
    BTL_INT ndata2=0;
    while (begin2 != end2)
    {
    	ndata2++;
    	begin2++;
    } 

    if (ndata1 == 0)
	FATAL_ERROR("There are no elements in 2nd container");

    if (ndata2%3 != 0)
    	FATAL_ERROR("The number of elements in the 2nd container is not a "   
    	    	    "multiple of 3");

    if (ndata1 != ndata2)
    	FATAL_ERROR("There are different numbers of elements in the two input "
    	            "containers");
    
    return ndata1/3;
}

 _BTL_END_NAMESPACE

#endif // btl_least_squares.h



