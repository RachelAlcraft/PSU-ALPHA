//
// btl_matrix_algorithms.h
//
// This file contains the generic numeric functions for manipulation of matrices
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

#if !defined (BTL_MATRIXALGORITHMS_H)
#define BTL_MATRIXALGORITHMS_H 1

/**#: [Description ="A collection of generic numerical algorithms for
       matrix algebra."]
    [Summary = "generic matrix algorithms"] 
    [Authors = "D.S.Moss, W.R.Pitt, I.Tickle, M.A.Williams"]
    [Files = "<A HREF=./btl/btl_matrix_algorithms.h>btl_matrix_algorithms.h</A>"]
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
// MATRIX ALGEBRA
//...............................................................................................
	    	/**#: [Description="Matrix multiplication.
			The output matrix is nrows1*ncols2 in size."]
	    	      [Restrictions="ncols1 must equal nrows2."] */

	template <class InputIterator1, class InputIterator2, class OutputIterator>
	OutputIterator 
	matrix_product(InputIterator1 first1, InputIterator1 last1, BTL_INT nrows1,
		       InputIterator2 first2, InputIterator2 last2, BTL_INT nrows2, 
		       OutputIterator result)
	{
	    BTL_INT ncols1 = (last1-first1)/nrows1;
	    BTL_INT ncols2 = (last2-first2)/nrows2;

	    #if defined(BTL_DEBUG_VERSION)
	        if (ncols1 != nrows2)
	    	   FATAL_ERROR("The number of rows in both matrices must be equal.");
	    #endif

	    BTL_INT i,j;
	    numeric_vector<typename iterator_traits<OutputIterator>::value_type> temp(nrows2);
	    typename iterator_traits<OutputIterator>::value_type tempresult = 0;
	    
	    for (j=0;j<ncols2;j++)
	    {
		copy_column(first2,last2,nrows2,temp.begin(),(j+1));
	    	for (i=0;i<nrows1;i++)
	    	{
                    tempresult=0;
	            *(result+((i*nrows1)+j)) +=
		        scalar_product((first1+(ncols1*i)),(first1+(ncols1*(i+1))),temp.begin(),tempresult);
		}
	    }
	    return (result+(nrows1*ncols2));
	}


//...............................................................................................
	    	/**#: [Description="Matrix multiplication: Transpose(M1)*M2.
			The output matrix is ncols1*ncols2 in size."]
	    	      [Restrictions="Both matrices must have the same number of rows"] */

	template <class InputIterator1, class InputIterator2, class OutputIterator>
	OutputIterator 
	matrixtranspose_matrix_product(InputIterator1 first1, InputIterator1 last1, BTL_INT nrows1,
				       InputIterator2 first2, InputIterator2 last2, BTL_INT nrows2, 
				       OutputIterator result)
	{
	    BTL_INT ncols1 = (last1-first1)/nrows1;
	    BTL_INT ncols2 = (last2-first2)/nrows2;

	    #if defined(BTL_DEBUG_VERSION)
	        if (nrows != nrows2)
	    	   FATAL_ERROR("The number of rows in both matrices must be equal.");
	    #endif
    
	    InputIterator2 j;
	    BTL_INT col1, col2, row = 1; 
	    for(col1=1;first1<last1;first1++,col1++)
	    {
		if (col1 > ncols1) col1 = 1, row = row + 1;
		for (j=(first2+((row-1)*ncols2)), col2=1; j<(first2+(row*ncols2)); j++, col2++)
		    *(result+((col1-1)*ncols1)+(col2-1)) += *first1 * *j;
	    }
	    return (result+(ncols1*ncols2));
	}

//...............................................................................................
	    	/**#: [Description="Matrix multiplication: M1*Transpose(M2)
			The output matrix is nrows1*nrows2 in size."]
	    	      [Restrictions="Both matrices must have the same number of columns"] */

	template <class InputIterator1, class InputIterator2, class OutputIterator>
	OutputIterator
	matrix_matrixtranspose_product(InputIterator1 first1, InputIterator1 last1, BTL_INT nrows1,
				       InputIterator2 first2, InputIterator2 last2, BTL_INT nrows2, 
				       OutputIterator result)
	{
            BTL_INT ncols = (last1-first1)/nrows1;
	    BTL_INT ncols2 = (last2-first2)/nrows2;
            typename iterator_traits<OutputIterator>::value_type tempresult=0;
            
	    #if defined(BTL_DEBUG_VERSION)
	        if ( ncols != ncols2)
	    	   FATAL_ERROR("The number of columns in both matrices must be equal.");
	    #endif

	    BTL_INT i,j;
		for(i=0;i < nrows1; i++)
		{
		    for(j=0;j < nrows2; j++)
		    {
	                tempresult = 0;
	                *(result+((i*nrows2)+j))
				 += scalar_product((first1+(i*ncols)),(first1+((i+1)*ncols)),
						   (first2+(j*ncols)),tempresult);
		    }
		}
	    return (result+(nrows1*nrows2));
	}

//..........................................................................................
	    	/**#: [Description="Matrix transpose. Can be used in-place."] */

	template <class InputIterator, class OutputIterator>
	OutputIterator transpose(InputIterator first, InputIterator last, BTL_INT nrows,
		        	 OutputIterator result)
	{
	    BTL_INT ncols= (last-first)/nrows;
	    BTL_INT i,j;
	    typename iterator_traits<InputIterator>::value_type temp;

	    if (nrows == ncols)
	    {
	    	for (i=0;i<nrows;i++)
	    	{
		    for (j=i;j<ncols;j++)
		    {
		        temp = *(first+(ncols*i)+j);
		        *(result+(ncols*i)+j) = *(first+i+(ncols*j));
		        *(result+i+(ncols*j)) = temp;
		    }
	    	}
	    }
	    else
	    {
		matrix<typename iterator_traits<InputIterator>::value_type> temp_matrix(first,last,nrows,ncols);
	    	for (i=1;i<=nrows;i++)
	    	{
		    for (j=1;j<=ncols;j++)
		    {
		        *(result+(nrows*(j-1))+(i-1)) = temp_matrix(i,j);
		    }
	    	}

	    }

	    return (result+(nrows*ncols));
	}

// AS ALL SWAPS OCCUR WITHIN PAIRS OF SYMMETRICALLY ARRANGED DIAGONALS THERE IS CERTAINLY A 
// MORE MEMORY EFFICIENT WAY OF DOING THIS TRANSPOSE IN THE NON-SQUARE CASE 

//..........................................................................................
	    	/**#: [Description="Calculate the mean of each column of a
	    	       	matrix and store the answers in a container of size >= ncols."] */
	
 	template<class InputIterator, class OutputIterator>
    	OutputIterator column_means(InputIterator first, InputIterator last, 
                                    BTL_INT nrows, OutputIterator result)
	{
	    BTL_INT ncols = (last-first)/nrows;
	    BTL_REAL reciprocal_nrows = 1.0/nrows;
	    OutputIterator begin = result;
	    OutputIterator end = result+ncols;
	    while (first != last)
	    {
	    	for ( result=begin; result!=end; result++, first++)
	    	    *result += *first;
	    }
	    	for ( result=begin; result!=end; result++, first++)
	    	    *result *= reciprocal_nrows;
	    return ++result;
	}

//...........................................................................................
	    	/**#: [Description="Copy a column of this matrix into another container.
				    Columns are numbered from 1. in the FORTRAN manner."] */

 	template<class InputIterator, class OutputIterator>
    	OutputIterator copy_column(InputIterator first, InputIterator last,
                                   BTL_INT nrows, OutputIterator result, BTL_INT column_index)
	{
	    BTL_INT ncols = (last-first)/nrows; 
	    if (column_index > ncols || column_index < 1)
	    	FATAL_ERROR("matrix index out of range");
    
	    for ( first+=(column_index-1) ; (first - last) < 0; result++, first+=ncols)
	    	*result = *first;
 
	    return ++result;
	}
	    	
//..........................................................................................
	    	/**#: [Description="Calculates matrix determinant"]
	    	      [Restrictions="Square matrices only."] */

 	template<class InputIterator, class T>
    	T determinant(InputIterator first, InputIterator last, BTL_INT nrows, T det)
	{
	    if ( (nrows*nrows) != (last-first))
	    	FATAL_ERROR("Cannot calculate determinant of non-square matrix");
    	
	    // For matrices with up to 4 rows use explicit Laplace expansion

	    if (nrows == 1) det += *first;

	    if (nrows == 2)
	    {
		det += (*first * *(first+3)) - (*(first+1) * *(first+2));
    	    }

	    if (nrows == 3)
	    {
		det += *first     * ((*(first+4) * *(first+8)) - (*(first+7) * *(first+5)))
		     - *(first+3) * ((*(first+1) * *(first+8)) - (*(first+7) * *(first+2)))
		     + *(first+6) * ((*(first+1) * *(first+5)) - (*(first+4) * *(first+2)));
    	    }

	    if (nrows == 4)
	    {
		det += (  ((*first      * *(first+5))  - (*(first+4)  * *(first+1))) 
			* ((*(first+10) * *(first+15)) - (*(first+14) * *(first+11))) )
 
		     + (  ((*(first+1)  * *(first+8))  - (*first      * *(first+9))) 
			* ((*(first+6)  * *(first+15)) - (*(first+14) * *(first+7)))  )
 
 		     + (  ((*first      * *(first+13)) - (*(first+1)  * *(first+12))) 
			* ((*(first+6)  * *(first+11)) - (*(first+10) * *(first+7)))  )
 
 		     + (  ((*(first+4)  * *(first+9))  - (*(first+5)  * *(first+8))) 
			* ((*(first+2)  * *(first+15)) - (*(first+14) * *(first+3)))  )
 
 		     + (  ((*(first+5)  * *(first+12)) - (*(first+4)  * *(first+13))) 
			* ((*(first+2)  * *(first+11)) - (*(first+10) * *(first+3)))  )
 
 		     + (  ((*(first+8)  * *(first+13)) - (*(first+9)  * *(first+12))) 
			* ((*(first+2)  * *(first+7))  - (*(first+6)  * *(first+3)))  ); 
    	    }

	    // For matrices with greater than 4 rows use iterative Laplace expansion 

	    if (nrows > 4)
	    {
	        int sign;
		matrix<T> Minor((nrows-1));	// creates appropriate amount of space
	    	typename matrix<T>::iterator end, begin = Minor.begin();
	    	for (BTL_INT i = 1; i <= nrows; i++)
	    	{
		    end = matrix_minor(first,last,nrows,begin,i,1);
		    sign = (2*(i%2)) - 1;
	            det += sign * *(first+(nrows*(i-1))) 
			 * determinant(begin,end,(nrows-1),T(0));
	    	}
	    }

	    return det;
	}

// CAN POSSIBLY BE MORE EFFICIENTLY WRITTEN USING TEMPLATE METAPROGRAMS.
// ALSO AS THE LAPLACE EXPANSION GROWS AS N! - FOR MATRICES WITH MORE THAN 6 OR 7 ROWS 
// WE SHOULD USE LU DECOMPOSITION. 

//..........................................................................................
    	    /**#: [Description="Find the minor of a matrix. 
		        The minor denoted by minormatrix(.,.,,.,i,j) is found by  
			eliminating the row i and column j from the parent matrix."] */

 	template<class InputIterator, class OutputIterator>
    	OutputIterator matrix_minor(InputIterator first, InputIterator last, BTL_INT nrows,
                             OutputIterator result, BTL_INT i, BTL_INT j)
	{
	    BTL_INT ncols= (last-first)/nrows;
	     
	    if ( i<1 || i>nrows || j<1 || j>ncols)
	    	FATAL_ERROR("One or both indexes are out of range");
    	
	    BTL_INT row,col;
	    BTL_INT nless = (ncols-1);
	    i-=1;
	    j-=1;
	    for (row = 0; row < nrows; row++) 
	    	for (col = 0; col < ncols; col++)
	    	{
	    	    if (row < i)
	    	    {
	    	    	if (col < j)
	    	    	    *(result + (nless*row + col)) 
						= *(first + (ncols*row + col));
	    	    	else if (col > j)
	    	    	    *(result + (nless*row + col - 1)) 
						= *(first + (ncols*row + col));
	    	    }    	    	
	    	    else if (row > i)
	    	    {
	    	    	if (col < j)
	    	    	    *(result + (nless*(row-1) + col)) 
						= *(first + (ncols*row + col));
	    	    	else if (col > j)
	    	    	    *(result + (nless*(row-1) + col - 1)) 
						= *(first + (ncols*row + col));
	    	    }
	    	}
	    return (result+((nrows-1)*(ncols-1)));
	}    	    	


//..........................................................................................
	    	/**#: [Description="Matrix adjoint"]
	    	      [Restrictions="Square matrices only"].*/

 	template<class InputIterator, class OutputIterator>
    	OutputIterator adjoint(InputIterator first, InputIterator last,
                             	       BTL_INT nrows, OutputIterator result)
 	{
	    if ((nrows*nrows) != (last-first))
	    	FATAL_ERROR("Cannot calculate adjoint of non-square matrix");
	    	
	    if (nrows < 3)
	    	FATAL_ERROR("Cannot calculate adjoint of matrix with less than 3 rows");
    
	    matrix<BTL_REAL> Minor((nrows-1));
	    typename matrix<BTL_REAL>::iterator end, begin = Minor.begin();
	    BTL_INT col,row;
	    int sign; // must declare as signed integer can't rely on automatic type conversion
	    for (col = 1; col <= nrows; col++)
	    	for (row = 1; row <= nrows; row++)
	    	{    	    	
		    end=matrix_minor(first,last,nrows,begin,row,col);
		    sign = (1-(2*((row+col)%2)));
	    	    *(result+(nrows*(col-1)+row-1)) 
			= sign * determinant(begin,end,(nrows-1),BTL_REAL(0)) ;
	    	}
	    return (result+(nrows*nrows));         
	}

//..........................................................................................
	    	/**#: [Description="Calculate inverse of a small square matrix
	    	       (inverse = adjoint/determinant). The result matrix should 
			be of a floating point type."] */
	    	    	    	    	
 	template<class InputIterator, class OutputIterator>
    	OutputIterator inverse_square(InputIterator first, InputIterator last,
                             	       BTL_INT nrows, OutputIterator result)
 	{
	    // Work out determinant
	    BTL_REAL det = BTL_REAL(0);
	    det = determinant(first,last,nrows,det);
	    BTL_REAL reciprocal_det = 1.0/det;
	    BTL_REAL small = numeric_limits<BTL_REAL>::epsilon();
	    if ( det*det < small*small)
	    	FATAL_ERROR("Matrix singular to machine precision: try SV_decomposition ");
      
	    // Find Adjoint of this matrix.
            OutputIterator end_of_result = adjoint(first,last,nrows,result);

	    // Divide by the determinant if non-zero
	    for (;result!=end_of_result;result++)
	        *result *= reciprocal_det;
    
	    return ++result;
	}

//..........................................................................................
	    	/**#: [Description="Eigenvalues and eigenvectors of a symmetric
	    	       matrix"]
	    	      [Restrictions="the input matrix must be real, square and 
	    	       symmetric, the output vector of eigen values will be nrows in size
		       and the output matrix of eigenvectors will be nrows*nrows in size."]*/

        // routine originally from IBM SSP manual (see p165) Ian Tickle April 1992,
        // (modified by  David Moss February 1993 and Mark Williams November 1998).
	//
	// n - number of rows in input matrix
	// a - an array of size n*(n+1)/2 containing lower triangle of the original 
	//     n*n matrix in the order:
	// 
	// 	      1      2      3    ...
	//     1    a[0]
	//     2    a[1]   a[2]
	//     3    a[3]   a[4]   a[5]   ...
	//
	//     NOTE a is used as working space and is overwritten.
	//     Eigenvalues are written into the diagonal elements of a
	//     i.e.  a[0]  a[2]  a[5]  for a 3*3 matrix.
	// 
	// r - Resultant matrix of eigenvectors stored columnwise in the same
	//     order as eigenvalues, initially set equal to identity matrix.

 	template<class InputIterator, class OutputMatrixIterator, class OutputVectorIterator>
    	void eigen_solution(InputIterator first, InputIterator last, BTL_INT nrows,
                            OutputMatrixIterator eigenvectors, OutputVectorIterator eigenvalues)
	{
	    int n = nrows;
	    BTL_INT size = n*n;
	    if (size != (last-first)) FATAL_ERROR("This matrix must be square");
           

	    BTL_INT column, row;
	    numeric_vector<typename iterator_traits<InputIterator>::value_type> a(nrows*(nrows+1)/2);
	    typename numeric_vector<typename iterator_traits<InputIterator>::value_type>::iterator trng = a.begin(); 
	    InputIterator in;

	    // Copy lower triangle of the input matrix to a numeric vector.

	    for( row = 1 ; row <= nrows; row++)
	        for(in=(first+(nrows*(row-1))); in!=(first+(nrows*(row-1)+row)); in++, trng++)
		    *trng = *in;


	    // The matrix that will hold the results is initially = I.

            for (int i=0; i< (nrows*nrows); i++)
 	        *(eigenvectors+i) = 0.0; 
            for (int i=0; i< (nrows*nrows); i += (nrows+1))
 	        *(eigenvectors+i)+= 1.0; 

	    // Setup variables  

	    BTL_INT i, il, ilq, ilr, im, imq, imr, ind, iq, j, k, km, l, ll, lm, lq, m, mm, mq;
	    BTL_REAL am, anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;
	    BTL_REAL theshold = numeric_limits<double>::epsilon();

  	    // Initial and final norms (anorm & anrmx).
 
	    anorm=0.0;
	    iq=0;
	    for (i=0; i<n; i++) for (j=0; j<=i; j++)
	    {
                if (j!=i) anorm+=a[iq]*a[iq];
	        ++iq;
	    }

	    if (anorm>0.0)
	    {
                anorm=sqrt(2.0*anorm);
	        anrmx=theshold*anorm/n;

		// Compute threshold and initialise flag.
 
		thr=anorm;
		while (thr>anrmx) // Compare threshold with final norm
	    	{ 
		    thr/=n;
		    ind=1;
	      	    while (ind) 
	      	    { 
			ind=0;
			l=0;
			while (l != n-1) // Test for l beyond penultimate column
			{
			    lq=l*(l+1)/2;
		           ll=l+lq;
		  	    m=l+1;
		  	    ilq=n*l;
		  	    while (m != n) // Test for m beyond last column
		  	    {
			    // Compute sin & cos.
				mq=m*(m+1)/2;
				lm=l+mq;
		             // if (fabs(a[lm])>=thr) underflow bug reported by Ethan Merritt, Wash U
		             // if ((a[lm]*a[lm])>=thr) fix suggested by EM but this changes precision R. Grosse-Kunstleve
		                if ((a[lm]*a[lm])>(thr*thr)) // RGK's and MAW's fix
		    		{
				    ind=1;
				    mm=m+mq;
				    x=0.5*(a[ll]-a[mm]);
			   	    y=-a[lm]/sqrt(a[lm]*a[lm]+x*x);
			      	    if (x<0.0) y=-y;
			      	    sinx=y/sqrt(2.0*(1.0+(sqrt(1.0-y*y))));
			      	    sinx2=sinx*sinx;
			      	    cosx=sqrt(1.0-sinx2);
			      	    cosx2=cosx*cosx;
			      	    sincs=sinx*cosx;

				   // Rotate l & m columns.
 
			      	    imq=n*m;
			      	    for (i=0; i<n; i++)
				    { 
					iq=i*(i+1)/2;
					if (i!=l && i!=m)
					{ 
					    if (i<m) im=i+mq;
					        else im=m+iq;
					    if (i<l) il=i+lq;
					        else il=l+iq;
					    x=a[il]*cosx-a[im]*sinx;
					    a[im]=a[il]*sinx+a[im]*cosx;
  		    	    		    a[il]=x;
					}
				    	ilr=ilq+i;
					imr=imq+i;
		 			x = (*(eigenvectors+ilr)*cosx)
					  - (*(eigenvectors+imr)*sinx);
					*(eigenvectors+imr) = (*(eigenvectors+ilr)*sinx)
							    + (*(eigenvectors+imr)*cosx);
					*(eigenvectors+ilr) = x;
				    }

				    x=2.0*a[lm]*sincs;
				    y=a[ll]*cosx2+a[mm]*sinx2-x;
				    x=a[ll]*sinx2+a[mm]*cosx2+x;
		   		    a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2);
				    a[ll]=y;
				    a[mm]=x;
				}
				m++;
			    }
                            l++;
			}
		    }
	        }
	    }
 


	    // Sort eigenvalues & eigenvectors in order of descending eigenvalue.

	    k=0;
	    for (i=0; i<n-1; i++)
	    {
		im=i;
	    	km=k;
	    	am=a[k];
	    	l=0;
	    	for (j=0; j<n; j++)
	    	{
		    if (j>i && a[l]>am)
		    {
			im=j;
			km=l;
			am=a[l];
		    }
	    	    l+=j+2;
		}
	    	if (im!=i)
	    	{
		    a[km]=a[k];
	            a[k]=am;
	            l=n*i;
	            m=n*im;
		    for (j=0; j<n; j++)
		    {
			am=*(eigenvectors+l);
		  	*(eigenvectors+(l++)) = *(eigenvectors+m);
		  	*(eigenvectors+(m++)) = am;
		    }
	    	}
	    	k+=i+2;
	    }


	    // place sorted eigenvalues into the matrix_vector structure


	    for (j=0, k=0; j<nrows; j++)
	    { 
	    	eigenvalues[j]=a[k];
	      	k+=j+2;
	    }
	}
	    	
//..........................................................................................
	    	/**#: [Description="matrix inverse. KNOWN BUGS"] */

    	template <class InputIterator, class OutputIterator>
        OutputIterator inverse(InputIterator first, InputIterator last, BTL_INT nrows,
		                 OutputIterator result)
	{
	    return _inverse(first,last,nrows,result,
	                    typename iterator_traits<OutputIterator>::value_type(result));
	}

    	template <class InputIterator, class OutputIterator, class T>
        OutputIterator _inverse(InputIterator first, InputIterator last, BTL_INT nrows,
		                 OutputIterator result, T*)
	{
	    BTL_INT ncols = (last-first)/nrows;
	    numeric_vector<BTL_REAL> vecL(ncols);
	    matrix<BTL_REAL> matV(ncols,ncols);
    
	    // Calculate SVD
	    sv_decomposition(first,last,nrows,result,vecL.begin(),matV.begin());
    	
	    // Find maximum singular value.
	    BTL_REAL vecLmax = *(max_element(vecL.begin(),vecL.end()));

	    // Set very small values to zero in the matLinv matrix to reduce rounding errors
	    matrix<BTL_REAL> matLinv(ncols,ncols);
	    BTL_REAL vecLmin = vecLmax * numeric_limits<T>::epsilon();
	    for (BTL_INT i=1; i<=ncols; i++)
	    {
	    	if (vecL(i) < vecLmin)
	      	    matLinv(i,i) = 0.0;
	    	else
	      	    matLinv(i,i) = 1.0/vecL(i);
	    }

	    matrix<BTL_REAL> temp(ncols,nrows);
	    matrix_matrixtranspose_product(matLinv.begin(), matLinv.end(), matLinv.rows(),
					   result,(result+(nrows*ncols)),nrows,temp.begin());

	    return         matrix_product(matV.begin(),matV.end(),matV.rows(),
			   	          temp.begin(),temp.end(),temp.rows(),result);
	}

//..........................................................................................	
	    	/**#: [Description="Single value decomposition.<P> 
	    	       A = U * diag(L) * transpose(V) where A is the input
	    	       matrix. KNOWN BUGS"]
	    	      [Restrictions="On input, the size of vector L (representing a diagonal 
			matrix)must be at least equal to the number of columns in the input
		        matrix and matrix V must be a square matrix with dimensions 
			(ncols,ncols), matrix U is the same size as the input matrix."] */

    	template <class InputIterator, class OutputMatrixIterator, class OutputVectorIterator>
        void sv_decomposition(InputIterator first, InputIterator last, BTL_INT nrows,
		 	      OutputMatrixIterator matU, OutputVectorIterator vecL, 
			      OutputMatrixIterator matV)
	{
            return _sv_decomposition(first,last,nrows,matU,vecL,matV,
	                             typename iterator_traits<InputIterator>::value_type(first));
	}

    	template <class InputIterator, class OutputMatrixIterator, class OutputVectorIterator, class T>
        void _sv_decomposition(InputIterator first, InputIterator last, BTL_INT nrows,
		 OutputMatrixIterator matU, OutputVectorIterator vecL, OutputMatrixIterator matV, T*)
	{

	    // It may be shown that any matrix A with nrows >= ncols can be
	    // expressed as the product of three matrices.
	    // A = U * L * t(V), where U(nrows,ncols) and V(ncols,ncols) are 
	    // orthogonal matrices i.e t(U)U=t(V)V=I
	    // and L(ncols,ncols) is a diagonal matrix and t() means transpose.
	    // => t(A).A = t(ULt(V)).(ULt(V)) = Vt(L)t(U).ULt(V) = V t(L)L t(V)
	    // => t(A).A.V = V.t(L)L, which is a set of eigenproblems  

	    BTL_INT ncols = (last-first)/nrows;
	    OutputVectorIterator beginL = vecL;    
	    OutputVectorIterator endL   = vecL; endL += ncols; 

	    matrix<BTL_REAL> ATA(ncols,ncols);
	    matrixtranspose_matrix_product(first,last,nrows,first,last,nrows,ATA.begin());
	    eigen_solution(ATA.begin(),ATA.end(),ATA.rows(),matV,vecL);

	    // Retrieve L from t(L).L
	    for (; vecL!=endL; vecL++)
	    	*vecL = sqrt(*vecL);
    
	    // U = A * V * inv(L)

	    matrix_product(first,last,nrows,matV,(matV+(ncols*ncols)),ncols,matU);

	    numeric_vector<BTL_REAL> reciprocal_L(ncols,1.0);
	    numeric_vector<BTL_REAL>::iterator i;
	    for(i=reciprocal_L.begin(); beginL!=endL; beginL++,i++)
		*i /= *beginL;

	    BTL_INT col, row;
	    for(row=1;row<=nrows;row++)
	        for (col=1;col<=ncols;col++)
	 	    *matU++ *= reciprocal_L(col);
	}
//
// I THINK THIS PRODUCES A DECOMPOSITION OF A ROW_WISE PERMUTATION OF A.
// OK? MAYBE BUT NOT EFFICIENT. ALSO BETTER TO DO V * inv(L) THEM MULTIPLY BY A 
//..........................................................................................
	   /**#: [Description ="
		    Inverts a symmetric positive definite matrix. Employs Cholesky
		    decomposition in three steps:<P> 
		    (1)M=LL(tr) (2)compute L(-1) (3)M(-1)=L(-1)(tr)L(-1) 
		    Input should contain numerical data that can be interpreted as a 
		    symmetric positive definite matrix.
		    The output is the inverse of this matrix in the form of a matrix
		    of the same size as the input and probably containing data of type
		    BTL_REAL."]*/

	template<class InputIterator, class OutputIterator>
	OutputIterator inverse_cholesky(InputIterator first, InputIterator last,
					  BTL_INT nrows, OutputIterator result)
	{
	    if ((last-first) != (nrows*nrows))
	    	FATAL_ERROR("This matrix must be square");

	    BTL_REAL a;
	    BTL_INT i,j,k;
	    matrix<BTL_REAL> m(first,last,nrows,nrows);

	    // Compute L where M=LL(tr) using Cholesky decomposition.
	    for (i=1; i<=nrows; i++)   
		for (j=1; j<=i; j++)
	    	{
		    a=m(i,j);
		    for (k=1; k<j; k++)
		    	a-=m(i,k)*m(j,k);
		    if (i==j)
		    	if (a > numeric_limits<BTL_REAL>::min() * 100.0)
		    	    m(i,i)=sqrt(a);
		    	else
		    	{
		    	    WARNING("A principal minor is singular\n");
		    	    m(i,i)=0.0;
		    	}
		    else
		    	m(i,j)=a/m(j,j);
	    	}

	    // compute L(-1)
	    for (j=1; j<=nrows; j++)
	    {
	    	m(j,j)=1.0/m(j,j);
	    	for (i=j+1; i<=nrows; i++)
	    	{
	    	    a=0.0;
	    	    for (k=j; k<i; k++)
	    	    	a-=m(i,k)*m(k,j);
	    	    m(i,j)=a/m(i,i);
	    	}
	    }

	    // compute M(-1) = L(-1)(tr)L(-1)
	    for (j=1; j<=nrows; j++)
	    	for (i=j; i<=nrows; i++)
	    	{
	    	    a=0.0;
	    	    for (k=i; k<=nrows; k++)
	    	    	a+=m(k,i)*m(k,j);
	    	    m(i,j)=m(j,i)=a;
	    	}

	    for (i=1;i<=nrows;i++)
		for (j=1;j<=nrows;j++)
		    *result++ = m(i,j);

	    return ++result;	
	}

// COULD BE MORE EFFICIENTLY CODED

_BTL_END_NAMESPACE

#endif

