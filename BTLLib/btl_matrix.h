//
// Matrix.h
//
// This file contains declarations for the matrix class.
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997-1998 Birkbeck College, Malet Street, London, U.K.
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

#if !defined (BTL_MATRIX_H)
#define BTL_MATRIX_H 1

#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <functional>
#include <iostream>
#include <limits>

#include "BTL.h"
#include "btl_numeric_vector.h"
#include "btl_vector_algorithms.h"


_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description ="
     This template class represents a numerical Matrix of any dimension. 
     Matrices dimension 1 are more efficiently modelled using the numeric_vector.<P>

     There are two types associated with this class:<P>

     <CODE>value_type</CODE> : this type is the same as the 
     numeric_vector::value_type and defines the type of the elements of the Matrix.<P>
 
     <CODE>size_type</CODE> : this type is the same as the numeric_vector::size_type and
     defines the type of the Matrix indexes. It should always be an unsigned
     integer type."]
 
    [Summary = "a 2-dimensional array of real numbers of any size"] 
    [Authors = "D.S.Moss, W.R.Pitt, I.Tickle, M.A.Williams."]
    [Files = "<A HREF=./btl/btl_matrix.h>btl_matrix.h</A>"]
    [Friends="Friend equivalents to some functions are available and are 
    	      documented with these functions. Also available is: friend
    	      ostream& operator<<(ostream &os, const matrix<T> &m); "]
    [Dependencies="<A HREF=#numeric_vector>btl_numeric_vector.h</A>,
                   <A HREF=#vector_algorithms>btl_vector_algorithms.h</A>"]

*/

template<class T = BTL_REAL>
class matrix
{
public:

// Define equivalent types to those found in the STL container classes 

    typedef vector<T>  					container;
    typedef typename container::value_type  		value_type;
    typedef typename container::size_type  		size_type;    
    typedef typename container::reference 		reference;
    typedef typename container::const_reference 	const_reference;
    typedef typename container::iterator 		iterator;
    typedef typename container::const_iterator		const_iterator;
    typedef typename container::reverse_iterator	reverse_iterator;
    typedef typename container::const_reverse_iterator	const_reverse_iterator;
    typedef typename container::difference_type 	difference_type;

// Define equivalent type to that found in the TNT container classes 

    typedef value_type  	    	        	element_type;
    	    	    	    	    	    	
private:

	container mat;  	    		// The matrix.
	size_type nrows, ncols;   		// The number of rows and columns.

public:    	


//.............................................................................
// CONSTRUCTORS/DESTRUCTOR.
//.............................................................................

	    /**#: [Description="Constructor for default 3x3 matrix and
	    	   		initialise elements to zero."] */

	matrix(): nrows(3), ncols(3)
	{ 
    	    mat.insert(mat.begin(), 9, value_type(0));
	}

//.............................................................................
	    /**#: [Description="Constructs an identity matrix of size p*p."] */

	matrix(BTL_INT p): nrows(p), ncols(p)
	{

	#if defined(BTL_DEBUG_VERSION)
    	    if (p<1) FATAL_ERROR("The dimensions must be greater than 0");
	#endif

    	    mat.reserve(p*p);
    	    mat.push_back(value_type(1));
    	    for (size_type i=1; i<p; i++)
    	    {
       		mat.insert(mat.end(), p, value_type(0)); 
       		mat.push_back(value_type(1));
    	    }
	}	

//.............................................................................
	    /**#: [Description="Constructs a matrix with p rows and q 
	    	   		columns. Initialises each element to a value v.
                   		The default value is zero"] */

	matrix(BTL_INT p, BTL_INT q, const value_type& v = value_type(0))
		  : nrows(p),ncols(q)
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (p<1 || q<1) FATAL_ERROR("Both dimensions must be greater than 0");
	#endif

	    size_type elements = p*q;
	    mat.insert(mat.begin(), elements, v);
	}
	
//.............................................................................
	    /**#: [Description="Constructs a matrix with p rows and q 
	    	       		columns. Elements are initialised to those in given array.
                       		This is obsolete and will be deleted in future releases"] */

	matrix(const value_type* const array,
	           BTL_INT p, BTL_INT q): nrows(p), ncols(q)
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (p<1 || q<1) FATAL_ERROR("Both dimensions must be greater than 0");
	#endif

	    size_type elements = p*q;
	    mat.reserve(elements);
	    for (size_type i=0; i<elements; i++)
	    	mat.push_back(array[i]);
	} 

//.............................................................................
    	    /**#: [Description="Construct a matrix with p rows from a container 
                                given its beginning and end iterators. "] */

        template <class InputIterator>
        matrix(InputIterator first, InputIterator last, BTL_INT p, BTL_INT q) 
                  : nrows(p), ncols(q)
	{
	    mat.reserve((p*q));
	    for(;first!=last;first++) 
	    	mat.push_back(*first);
	}

//.............................................................................
	    	/**#: [Description="Constructor for a 3x3 Matrix with
	    	       initialisation"] */

	matrix(const value_type& e1, const value_type& e2, const value_type& e3,
	       const value_type& e4, const value_type& e5, const value_type& e6,
	       const value_type& e7, const value_type& e8, const value_type& e9)
		: nrows(3), ncols(3)
	{ 
	    mat.reserve(9);
	    mat.push_back(e1); mat.push_back(e2); mat.push_back(e3);
	    mat.push_back(e4); mat.push_back(e5); mat.push_back(e6);
	    mat.push_back(e7); mat.push_back(e8); mat.push_back(e9);
	}


//.............................................................................
	    	/**#: [Description="Copy constructor"] */

	matrix(const matrix<T> &m) 
                   : mat(m.mat), nrows(m.nrows), ncols(m.ncols) {}

//.............................................................................
	    	/**#: [Description="Destructor"] */
	~matrix() {}


//.............................................................................
// OVERLOADED OPERATORS.
//.............................................................................
// ACCESS FUNCTIONS
//.............................................................................
	    	/**#: [Description="Returns an iterator that points to the first
	    	       element in the Matrix. (There is also a const version of
		       this function.)"] */

    	iterator
    	begin() { return mat.begin(); }

    	const_iterator
    	begin() const { return mat.begin(); }

//.............................................................................
	    	/**#: [Description="Returns an iterator that can be used
	    	       in a comparison for ending a traversal through this
	    	       Matrix. (There is also a const version of
		       this function.)"] */

    	iterator
    	end() { return mat.end(); }

    	const_iterator
    	end() const { return mat.end(); }

//.............................................................................
	    	/**#: [Description="Size of Matrix (the number of rows times
	    	       the number of columns)"] */

    	size_type
    	size() const { return mat.size; }

//.............................................................................
	    	/**#: [Description="This function is added to give
		       compatibility with TNT."] */
    	size_type
	lbound() const { return 1; }

//.............................................................................
    	    	/**#: [Description="Read number of rows (for compatability 
		       with TNT)."] */ 
    	size_type 
    	num_rows() const {return nrows;}


//.............................................................................
    	    	/**#: [Description="Read number of columns (for compatability 
		       with TNT)."] */ 
    	size_type 
    	num_cols() const {return ncols;}

//.............................................................................
    	    	/**#: [Description="Read number of rows."] */ 
    	size_type 
    	rows() const {return nrows;}

//.............................................................................
    	    	/**#: [Description="Set number of columns."] */ 
    	void 
    	cols(BTL_INT n) {ncols = n;}

//.............................................................................
    	    	/**#: [Description="Set number of rows."] */ 
    	void 
    	rows(BTL_INT n) {nrows = n;}

//.............................................................................
    	    	/**#: [Description="Read number of columns."] */ 
    	size_type 
    	cols() const {return ncols;}

	    	
//.............................................................................
	    	/**#: [Description="This function is added to give 
		       compatibility with TNT."] */
    	size_type 
	dim(size_type d) const 
    	{
	#if defined(DEBUG_VERSION)
            if ( d < 1 || d > 2) FATAL_ERROR("d must be equal to 1 or 2");
	#endif
            return (d==1) ? nrows : ((d==2) ? ncols : 0); 
    	}

//.............................................................................
	    	/**#: [Description="Returns a matrix element given its indices
	    	       Fortran style e.g. x(i,j) i = row, j = col N.B. (i,j >= 1).
                       (There is also a const (read only) version of this function.)"] */

	value_type& 
	operator()(const size_type& i, const size_type& j)
	{
	
	#if defined(BTL_DEBUG_VERSION)
	    if (i<1 || j<1 || i>nrows || j>ncols) 
     		FATAL_ERROR("Matrix indices out of range");
	#endif

	    return mat[(ncols*(i-1) + j - 1)]; 
	}

	value_type
	operator()(const size_type& i, const size_type& j) const
	{
	
	#if defined(BTL_DEBUG_VERSION)
	    if (i<1 || j<1 || i>nrows || j>ncols) 
     		FATAL_ERROR("Matrix indices out of range");
	#endif

	    return mat[(ncols*(i-1) + j - 1)]; 
	}

//.............................................................................
	    	/**#: [Description="Returns an iterator that points to the first
	    	       element of a given row. N.B. ( 0 <= i < nrows ). (There
		       is also a const (read only) version of this function.)"] */
	iterator
	operator[](const size_type& i)
	{ 

	#if defined(BTL_DEBUG_VERSION)
	    if (i>=nrows) FATAL_ERROR("Matrix index out of range");
	#endif

	    iterator j = mat.begin() + ncols*i;
	    return j;
	}

	const_iterator
	operator[](const size_type& i) const
	{ 

	#if defined(BTL_DEBUG_VERSION)
	    if (i>=nrows) FATAL_ERROR("Matrix index out of range");
	#endif

	    const_iterator j = mat.begin() + ncols*i;
	    return j;
	}

//.............................................................................
	    	/**#: [Description="Returns a matrix element given its indices
	    	       C-style e.g. x[i,j] i = row, j = col N.B. (i,j >= 0) NON STANDARD.
                       (There is also a const (read only) version of this function.)"]

	value_type&
	operator[](const size_type& i, const size_type& j)
	{
	
	#if defined(BTL_DEBUG_VERSION)
	    if (i<0 || j<0 || i>(nrows-1) || j>(ncols-1))
     		FATAL_ERROR("Matrix indices out of range");
	#endif

	    return mat[(ncols*i + j)];
	}

	value_type
	operator[](const size_type& i, const size_type& j) const
	{
	
	#if defined(BTL_DEBUG_VERSION)
	    if (i<0 || j<0 || i>(nrows-1) || j>(ncols-1))
     		FATAL_ERROR("Matrix indices out of range");
	#endif

	    return mat[(ncols*i + j)];
	}
 */
//.............................................................................
// MATRIX ALGEBRA
//.............................................................................

	    	/**#: [Description="Matrix assignment"] */
	matrix<T>& 
	operator=(const matrix<T> &m)
	{
	    if(this != &m)	 // check for self assignment
	    {
	      mat = m.mat;       // memory (de)allocation handled by the STL vector
	      nrows = m.nrows;
	      ncols = m.ncols;
	    }  
	    return *this;
	}

//.............................................................................
	    	/**#: [Description="Matrix multiplication e.g. matrix m1,m2,m3;
	    	       .... m3 = m1 * m2;"]*/
	matrix<T>
	operator*(const matrix<T>& m) const
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (ncols != m.num_rows()) 
	        FATAL_ERROR("The number of columns in the first Matrix must equal the "
	    	                "number of rows in the second");
	#endif
 
	    matrix<T> result(nrows,m.num_cols());
	    const_iterator j = m.begin();
	    iterator k;
	    size_type column;
	    size_type row = 0, element = 0;
	    for (const_iterator i = mat.begin(); i < mat.end(); i++, element++)
	    {
	        if (int(element/ncols) > row)
	        {
	            j = m.begin();
	            row++;
	        }
	        for (column = 0, k = result[row] ; column < m.num_cols(); column++, j++, k++)
	            *k += *i * *j;
	    }
	    return result;
	}

//.............................................................................
	    	/**#: [Description="Postmultiplication of a Matrix by a
	    	       numeric_vector. e.g. matrix m;
	    	       numeric_vector v1,v2; ... v2 = m * v1; "]
         	      [Restrictions="The size of v must equal the number of
         	       columns in this Matrix."] */
	numeric_vector<T>
	operator*(const numeric_vector<T> &v) const
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (ncols != v.size()) 
	    	FATAL_ERROR("The size of the numeric_vector<T> must equal the number of "
    	            "columns in the Matrix");
	#endif

	    numeric_vector<T> result(nrows);
	    const_iterator i;
	    const_iterator j = mat.begin();
	    for (iterator k = result.begin(); k < result.end(); k++)
	        for (i = v.begin(); i < v.end(); i++, j++)
	            *k += *i * *j; 
	    return result;
	}

//.............................................................................
	    	/**#: [Description="Multiple each element by a number."] */
	matrix<T>&
	operator*=(const value_type& v)
	{
	    for (iterator i=mat.begin(); i!=mat.end(); i++)
		*i *= v;
	    return *this;
	}


//.............................................................................
	    	/**#: [Description="Multiple each element by a number"] */
	matrix<T>
	operator*(const value_type& v) const
	{
	    matrix<T> result = *this;
	    for (iterator i=result.begin(); i!=result.end(); i++)
		*i *= v;
	    return result;
	}

//.............................................................................
	    	/**#: [Description="Divide each element by a number"] */
	matrix<T>&
	operator/=(const value_type& v)
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (v == 0.0) WARNING("Divide by zero");
	#endif

	    BTL_REAL y_reciprocal = 1.0 / v;
	    for (iterator i=mat.begin(); i!=mat.end(); i++)
	    	*i *= y_reciprocal;
	    return *this;
	}


//.............................................................................
	    	/**#: [Description="Divide each element by a number"] */
	matrix<T>
	operator/(const value_type& v) const
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (v == 0.0) WARNING("Divide by zero");
	#endif

	    matrix<T> result= *this;
	    BTL_REAL y_reciprocal = 1.0 / v;
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i *= y_reciprocal;
	    return result;
	}


//.............................................................................
	    	/**#: [Description="Subtraction of a number from each element."]
	    	      */
    	matrix<T>&
    	operator-=(const value_type &v)
	{
	    for (iterator i=mat.begin(); i!=mat.end(); i++)
	    	*i -= v;
	    return *this;
	}    	

//.............................................................................
	    	/**#: [Description="Subtraction of a number from each element."]
	    	      */
	matrix<T>
	operator-(const value_type &v) const
	{
	    matrix<T> result = *this;
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i -= v;
	    return result;
	}    	


//.............................................................................
	    	/**#: [Description="Addition of a number to each element."] */
    	matrix<T>&
    	operator+=(const value_type &v)
	{
	    for (iterator i=mat.begin(); i!=mat.end(); i++)
	    	*i += v;
	    return *this;
	}    	


//.............................................................................
	    	/**#: [Description="Addition of a number to each element."] */
	matrix<T>
	operator+(const value_type &v) const
	{
	    matrix<T> result = *this;
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i += v;
	    return result;
	}    	


//.............................................................................
	    	/**#: [Description="Subtraction of a vector from the rows
	    	       of a matrix"]
	    	      [Restrictions="The the size of numeric_vector v must equal the
	    	       number of columns in this matrix."] */
    	matrix<T>&
    	operator-=(const numeric_vector<T> &v)
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != ncols)
	    	FATAL_ERROR("The the size of the input numeric_vector<T> must equal the number of "
	    	            "columns in the Matrix");
	#endif

	    iterator i = mat.begin();
	    const_iterator j;
	    for (int k=0; k<nrows; k++)
	    	for (j=v.begin(); j!=v.end(); j++, i++)
	    	     *i -= *j;
    	    
	    return *this;
	}


//.............................................................................
	    	/**#: [Description="Subtraction of a vector from the rows
	    	       of a matrix"]
	    	      [Restrictions="The the size of numeric_vector v must equal the
	    	       number of columns in this matrix."] */
	matrix<T>
	operator-(const numeric_vector<T> &v) const
	{
	    matrix<T> result = *this;
	    result -= v;
	    return result;
	}


//.............................................................................
	    	/**#: [Description="Addition of a vector to the rows
	    	       of a matrix"]
	    	      [Restrictions="The the size of numeric_vector v must equal the
	    	       number of columns in this matrix."] */
    	matrix<T>&
    	operator+=(const numeric_vector<T> &v)
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != ncols)
	    	FATAL_ERROR("The the size of the input numeric_vector<T> must equal the number of "
	    	            "columns in this Matrix");
	#endif

	    iterator i = mat.begin();
	    const_iterator j;
	    for (int k=0; k<nrows; k++)
	    	for (j=v.begin(); j!=v.end(); j++, i++)
	    	     *i += *j;
	    return *this;
	}


//.............................................................................
	    	/**#: [Description="Addition of a vector to the rows
	    	       of a matrix"]
	    	      [Restrictions="The the size of numeric_vector v must equal the
	    	       number of columns in this matrix."] */
	matrix<T>
	operator+(const numeric_vector<T> &v) const
	{
	    matrix<T> result = *this;
	    result += v;
	    return result;
	}


//.............................................................................
	    	/**#: [Description="matrix subtraction"]
	    	      [Restrictions="The the size of matrix m must equal the
	    	       size of this matrix."] */
    	matrix<T>&
    	operator-=(const matrix<T> &m)
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (nrows != m.num_rows() || ncols != m.num_cols())
	    	FATAL_ERROR("Both matrices must be the same size");
	#endif

	    iterator i;
	    const_iterator j;
	    for (i=mat.begin(), j=m.begin(); i!=mat.end(); i++, j++)
	    	*i -= *j;
	    return *this;
	}


//.............................................................................
	    	/**#: [Description="Matrix subtraction"]
	    	      [Restrictions="The the size of Matrix m must equal the
	    	       size of this Matrix."] */
	matrix<T>
	operator-(const matrix<T> &m) const
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (nrows != m.num_rows() || ncols != m.num_cols())
	    	FATAL_ERROR("Both matrices must be the same size");
	#endif

	    matrix<T>  result = *this;
	    const_iterator i;
	    iterator j;
	    for (i=m.begin(), j=result.begin(); i!=m.end(); i++, j++)
	    	*j -= *i;
	    return result;
	}


//.............................................................................
	    	/**#: [Description="Matrix addition"]
	    	      [Restrictions="The the size of Matrix m must equal the
	    	       size of this Matrix."] */
    	matrix<T>&
    	operator+=(const matrix<T> &m)
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (nrows != m.num_rows() || ncols != m.num_cols())
	   	FATAL_ERROR("Both matrices must be the same size");
	#endif

	    iterator i;
	    const_iterator j;
	    for (i=mat.begin(), j=m.begin(); i!=mat.end(); i++, j++)
	    	*i += *j;
	    return *this;
	}


//.............................................................................
	    	/**#: [Description="Matrix addition"]
	    	      [Restrictions="The the size of Matrix m must equal the
	    	       size of this Matrix."] */
	matrix<T>
	operator+(const matrix<T> &m) const
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (nrows != m.num_rows() || ncols != m.num_cols())
	    	FATAL_ERROR("Both matrices must be the same size");
	#endif

	    matrix<T>  result = *this;
	    const_iterator i;
	    iterator j;
	    for (i=m.begin(), j=result.begin(); i!=m.end(); i++, j++)
	    	*j += *i;
	    return result;
	}


//.............................................................................
	    	/**#: [Description="Matrix multiplication Transpose(m1) * m2
	    	       N.B. this is less efficient than (but not equivalent to)
	    	       operator%"]
	    	      [Restrictions="Both *this and the input Matrix must have
	    	       the same number of rows"] */
	matrix<T>
	operator&(const matrix<T>& m) const
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (m.num_rows() != nrows)
	    	FATAL_ERROR("The input Matrix must have the same number of rows as this Matrix");
	#endif
    
	    matrix<T> result(ncols,m.num_cols());
	    numeric_vector<T> temp(nrows);
	    size_type i,j,k,iarr;
	    for (i=0, iarr=0; i<ncols; i++, iarr=result.ncols*i)
	    {
	    	for (j=0; j<nrows; j++) temp[j] = mat[ncols*j + i];
	    	for (j=0; j<m.ncols; j++)
	    	    for (k=0; k<nrows; k++)    	    
	    	    	result.mat[iarr + j] += temp[k] * m.mat[m.ncols*k + j]; // m1(i,j)+= *this(k,i) * m(k,j)
	    }   	
	    return result;
	}

//.............................................................................
	    	/**#: [Description="Matrix multiplication m1 * Transpose(m2)
	    	       N.B. this is more efficient than (but not equivalent to)
	    	       operator&"]
	    	      [Restrictions="Both *this and the input Matrix must have
	    	       the same number of columns"] */
	matrix<T> 
	operator%(const matrix<T>& m) const
	{

	#if defined(BTL_DEBUG_VERSION)
	    if (m.num_cols() != ncols)
	    	FATAL_ERROR("Both this and the input Matrix must have the same number "
	    	            "of columns");
	#endif
    
	    matrix<T> result(nrows,m.num_rows());
	    for (size_type i=1; i<=nrows; i++)
	    	for (size_type j=1; j<=m.nrows; j++)
	    	    for (size_type k=1; k<=ncols; k++)
	    	    	result(i,j) += (*this)(i,k) * m(j,k);

	    return result;
	}


//.............................................................................
    	    	/**#: [Description="Equality operator"]*/
	bool 
	operator==(const matrix<T>& m) const
	{ 
	    if ( nrows != m.num_rows() || ncols != m.num_cols() )
	    	return false;
    	
	    iterator i;
	    const_iterator j;
	    for (i=mat.begin(), j=m.begin(); i!=mat.end(); i++, j++)
            {
	    	if (*i != *j) return false;
            }
            return true;
	}

//.............................................................................
// I/O FUNCTIONS
//.............................................................................
// ostream operator

	friend ostream&
	operator<<(ostream &os, const matrix<T> &m)
	{
	        typedef typename matrix<T>::size_type size_type;

		os.setf(ios::showpoint);
		os.setf(ios::fixed, ios::floatfield);

		for ( size_type i=1; i<=m.nrows; i++)
		{
			os << '\n';
			for (size_type j=1; j<=m.ncols; j++)
				os << m(i,j) << " ";
		}
		os << '\n';
		return os;
	}

//.............................................................................
// istream operator

    	friend istream&
	operator>>(istream &is, matrix<T> &m)
	{
    
	    typedef typename matrix<T>::size_type size_type;
	    typedef typename matrix<T>::value_type value_type;

	    size_type M, N;

	    is >> M >> N;

	    if (M<1 || N<1) 
	    {
	    	WARNING("Both dimensions must be greater than 0");
		return is;
	    }

//	    if ( !(M == m.nrows && N == m.ncols) )
//	    {
//	        ~m();
//		matrix<T> m(M,N);
//	    	m.nrows = M;
//		m.ncols = N;
//	    }

	    for (size_type i=0; i<M; i++)
	        for (size_type j=0; j<N; j++)
	        {
	            is >>  m(i,j);
	        }

	    return is;
	}

//.............................................................................
// OTHER FUNCTIONS.
//.............................................................................

    	
}; // class matrix

_BTL_END_NAMESPACE

#endif

