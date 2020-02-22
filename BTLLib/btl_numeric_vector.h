//
// Vector.h
//
// This file contains declarations for the BTL_Vector class.
// This code is part of the Bioinformatics Template Library (BTL).
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



#if !defined(BTL_NUMERICVECTOR_H)
#define BTL_NUMERICVECTOR_H 1

#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

#include "BTL.h"

_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description =" this template class represents a numerical vector of any 
    size.<P>

                             V = xi + yj + zk + ...<P>

    However, its primary use will probably in representing the position of atoms
    in three dimensions. The default contructor contructs a 3 element numeric_vector.<P>

    There are several types associated with this class:<P>

    <CODE>value_type</CODE> : by default the same as the type BTL_REAL which is defined in
                              BTL.h<P>
    <CODE>iterator</CODE> : this is a STL type iterator. It can be used to
    randomly access the elements of the numeric_vector e.g.<P><CODE>

    	    	    numeric_vector<> v1;
    	    	     	...
    	    	     for (numeric_vector<>::iterator i = v1.begin(); i != v1.end(); i++)
    	    	    	cout << *i;

    <P>const_iterator</CODE> : used to access a const numeric_vector.<P>
    <CODE>size_type</CODE> : this is an unsigned integer type that represents
    the size of any numeric_vector object."]
    [Summary = "a vector of real numbers of any size"]
    [Authors = "D.S.Moss, W.R.Pitt, M.A.Williams"]
    [Files = "<A HREF=./btl/btl_numeric_vector.h>btl_numeric_vector.h</A>"]
    [Friends="Friend equivalents to some functions are available and are
    	      documented with these functions. Also available is: friend
    	      ostream& operator<<(ostream &os, const numeric_vector<T> &m); "]
    [Dependencies="None"]*/

template<class T = BTL_REAL>
class numeric_vector
{
public:

// Define equivalent types to those found in the STL container classes 

    	typedef vector<T>  					container;
    	typedef typename container::value_type  		value_type;
    	typedef typename container::size_type  			size_type;    
    	typedef typename container::reference 			reference;
    	typedef typename container::const_reference 		const_reference;
    	typedef typename container::iterator 			iterator;
    	typedef typename container::const_iterator		const_iterator;
    	typedef typename container::reverse_iterator		reverse_iterator;
    	typedef typename container::const_reverse_iterator	const_reverse_iterator;
    	typedef typename container::difference_type 		difference_type;

// Define equivalent type to that found in the TNT container classes 

    	typedef value_type  	    	        		element_type;

private:

    	container vec;					// The initially empty vector

public:

//.............................................................................
// Constructors/destructor

	    /**#: [Description="Constructs a 3D vector
	                        and initialises each value to zero."]*/

    	numeric_vector() { vec.insert(vec.begin(), 3, value_type(0)); }

    	    /**#: [Description="Constructs a 3D vector with initialization
    		   and initialises the values with the numbers given."]*/

    	numeric_vector(const value_type& p, const value_type& q, const value_type& r)
    	      	{ vec.push_back(p); vec.push_back(q); vec.push_back(r); }

    	    /**#: [Description="Constructor for a vector of n elements.
    				All n values initiated to a value v. 
                                The default value is zero ."]*/

    	numeric_vector(const size_type& n, const value_type& v = value_type(0)) 
		{ vec.insert(vec.begin(), n, v);}

    	    /**#: [Description="Construct from an array given a value_type*. 
                                This constructor is obsolete and will be deleted 
                                in future releases."] */

	numeric_vector(const value_type* const array, const size_type& p)
    	{
    	    for (size_type i=0; i<p; i++)
    	   	    vec.push_back(array[i]);
    	}

    	    /**#: [Description="Construct from a container given its beginning
                                and end iterators. "] */

   	template <class InputIterator>
    	numeric_vector(InputIterator first, InputIterator last) 
        {
	    vec.reserve((last-first));
	    for(;first!=last;first++) 
	    	vec.push_back(*first);
	}

    	    /**#: [Description="Copy constructor."]*/

    	numeric_vector(const numeric_vector<T> &v) : vec(v.vec) {}

    	    /**#: [Description="Destructor."]*/

    	~numeric_vector() {}

//.............................................................................
// Access functions
//.............................................................................

    	    /**#: [Description="Returns iterator that can be used to
    		   begin traversing through all elements of the numeric_vector.
		   (There is also a const version of this function.)"]*/
    	iterator
    	begin() { return vec.begin(); }

    	const_iterator
    	begin() const { return vec.begin(); }

//.............................................................................
    	    /**#: [Description="Returns an iterator can be used in
    		   comparison for ending traversal through the numeric_vector.
		   (There is also a const version of this function.)"]*/
    	iterator
    	end() { return vec.end(); }

    	const_iterator
    	end() const { return vec.end(); }

//.............................................................................
    	    /**#: [Description="Returns the element, i positions
    		   from the beginning of the numeric_vector, in constant time.
    		   (C-style i >= 0). (There is also a const version of this
		   function.)"]*/
 
    	value_type&
    	operator[](const size_type& i)
	{ 
	    #if defined(BTL_DEBUG_VERSION)
	        if (i>=size()) FATAL_ERROR("Index i is out of range");
	    #endif

	    return vec[i];
	}

    	value_type
    	operator[](const size_type& i) const
	{ 
	    #if defined(BTL_DEBUG_VERSION)
	        if (i>=size()) FATAL_ERROR("Index i is out of range");
	    #endif

	    return vec[i];
	}

//.............................................................................
    	    /**#: [Description="Returns the i-th element in the numeric_vector
    		   (Fortran style i >= 1). (There is  also a const version 
                   of this function)."]*/

    	value_type&
    	operator()(const size_type& i)
	{
	    #if defined(BTL_DEBUG_VERSION)
	        if (i>size() || i<1) FATAL_ERROR("Index i is out of range");
	    #endif

	    return vec[i-1];
	}

    	value_type
    	operator()(const size_type& i) const
	{
	    #if defined(BTL_DEBUG_VERSION)
	        if (i>size() || i<1) FATAL_ERROR("Index i is out of range");
	    #endif

	    return vec[i-1];
	}

//.............................................................................
    	    /**#: [Description="Returns the number of elements in the
    		   numeric_vector."]*/

    	size_type
    	size() const { return vec.size(); }

	    /**#: [Description="This function is added to give
		       compatibility with TNT."] */
    	size_type
    	dim() const {return vec.size();}

//.............................................................................
    	    /**#: [Description="Vector stream output."]*/

    	friend ostream&
    	operator<<(ostream &os, const numeric_vector<T> &v)
	{
	    if (v.vec.empty()) return os;  
	    os.setf(ios::showpoint);
	    os.setf(ios::fixed, ios::floatfield);
	    os  << "( " << v[0];
	    typename numeric_vector<T>::const_iterator i=v.begin(); 
	    for (i++; i!=v.end(); i++)
	    	    os << ", " << *i;
	    os << " )\n";
	    return os;
	}


//.............................................................................
// VECTOR ALGEBRA
//.............................................................................

    	    /**#: [Description="Vector assignment."]*/

    	numeric_vector<T>&
    	operator=(const numeric_vector<T> &v) 
	{ 
	    if(this != &v) vec = v.vec;return *this;
	}


//.............................................................................
    	    /**#: [Description="Vector addition."]*/

    	numeric_vector<T>
    	operator+(const numeric_vector<T> &v) const
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != size())
	    	WARNING("The numeric_vectors must be the same size for addition");
	#endif
    	
	    numeric_vector<T>  result(*this);
	    const_iterator i;
	    iterator j;
	    for (i=v.begin(), j=result.begin(); i!=v.end(); i++, j++)
	    	*j += *i;
	    return result;
	}

//.............................................................................
    	    /**#: [Description="Vector addition and assignment."]*/

    	numeric_vector<T>&
    	operator+=(const numeric_vector<T> &v)
	{ 
	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != size())
	    	WARNING("The numeric_vectors must be the same size for increment");
	#endif

	    iterator i;
	    const_iterator j;
	    for (i=vec.begin(), j=v.begin(); i!=vec.end(); i++, j++)
	    	*i += *j;
	    return *this;
	}

//.............................................................................
    	    /**#: [Description="Vector subtraction."]*/

    	numeric_vector<T>
    	operator-(const numeric_vector<T> &v) const
	{ 
	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != size())
	    	WARNING("The numeric_vectors must be the same size for subtraction");
	#endif

	    numeric_vector<T>  result = *this;
	    const_iterator i;
	    iterator j;
	    for (i=v.begin(), j=result.begin(); i!=v.end(); i++, j++)
	    	*j -= *i;
	    return result;
	}


//.............................................................................
    	    /**#: [Description="Vector negation (unary minus operator)."]*/

    	numeric_vector<T>
    	operator-() const
	{
	    numeric_vector<T>  negthis = *this;
	    for (iterator i=negthis.begin(); i!=negthis.end(); i++)
	    	*i = -(*i);
	    return negthis;
	}



//.............................................................................
    	    /**#: [Description="Vector subtraction and assignment."]*/

    	numeric_vector<T>&
    	operator-=(const numeric_vector<T> &v)
	{ 
	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != size())
	    	WARNING("The numeric_vectors must be the same size for decrement");
	#endif

	    iterator i;
	    const_iterator j;
	    for (i=vec.begin(), j=v.begin(); i!=vec.end(); i++, j++)
	    	*i -= *j;
	    return *this;
	}

//.............................................................................
    	    /**#: [Description="Add a scalar to each element in
    		   the numeric_vector."]*/

    	numeric_vector<T>
    	operator+(const value_type &s) const
	{
	    numeric_vector<T> result = *this;
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i += s;
	    return result;
	}    	

//.............................................................................
    	    /**#: [Description="Add and assign a scalar to each element in
    		   the numeric_vector."]*/
 
    	numeric_vector<T>&
    	operator+=(const value_type &s)
	{
	    for (iterator i=vec.begin(); i!=vec.end(); i++)
	    	*i += s;
	    return *this;
	}    	

//.............................................................................
    	    /**#: [Description="Subtract a scalar from each element in
    		   the numeric_vector."]*/
 
    	numeric_vector<T>
    	operator-(const value_type &s) const
	{
	    numeric_vector<T> result = *this;
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i -= s;
	    return result;
	}    	

//.............................................................................
    	    /**#: [Description="Subtract and assign a scalar from each element in
    		   the numeric_vector."]*/
 
    	numeric_vector<T>&
    	operator-=(const value_type &s)
	{
	    for (iterator i=vec.begin(); i!=vec.end(); i++)
	    	*i -= s;
	    return *this;
	}    	

//.............................................................................
    	    /**#: [Description="Vector scalar/dot product."]*/
 
    	value_type
    	operator*(const numeric_vector<T> &v) const
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (v.size() != size())
	    	WARNING("The input numeric_vector<T> must be the same size as this numeric_vector");
	#endif

	    value_type product = value_type(0);
	    const_iterator i,j;
	    for (i=vec.begin(), j=v.begin(); i!=vec.end(); i++, j++)
	    	product += *i * *j;

	    return product;
	}


//.............................................................................
    	    /**#: [Description="Vector cross product."]*/
 
    	numeric_vector<T>
    	operator%(const numeric_vector<T> &v) const
	{
	    // This works out the cross product of two vectors. These vectors must
	    // be of size 3.

	#if defined(BTL_DEBUG_VERSION)
	    if (size() != 3 || v.size() != 3)
	    	WARNING("Both numeric_vectors must be of size 3 for a valid cross product");
	#endif

	    numeric_vector<T> result;
	    result[0] = vec[1]*v[2] - vec[2]*v[1];
	    result[1] = vec[2]*v[0] - vec[0]*v[2];
	    result[2] = vec[0]*v[1] - vec[1]*v[0];
	    return result;
	}


//.............................................................................
    	    /**#: [Description="Vector multiplication by a scalar."]*/
 
    	numeric_vector<T>
    	operator*(const value_type& y) const
	{ 
	    numeric_vector<T> result = *this;
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i *= y;
	    return result;
	}

//.............................................................................
    	    /**#: [Description="Vector multiplication by a scalar and
    		  assignment."]*/
 
    	numeric_vector<T>&
    	operator*=(const value_type& y)
	{
	    for (iterator i=vec.begin(); i!=vec.end(); i++)
	    	*i *= y;
	    return *this;
	}

//.............................................................................
    	    /**#: [Description="Vector division by a scalar."]*/
 
    	numeric_vector<T>
    	operator/(const value_type& y) const
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (y == 0.0) WARNING("Divide by zero");
	#endif
    
	    BTL_REAL y_reciprocal = 1.0 / y;
	    numeric_vector<T> result = *this; 
	    for (iterator i=result.begin(); i!=result.end(); i++)
	    	*i *= y_reciprocal;
	    return result;
	}

//.............................................................................
    	    /**#: [Description="Vector division by a scalar and
    		   assignment."]*/
 
    	numeric_vector<T>&
    	operator/=(const value_type& y)
	{
	#if defined(BTL_DEBUG_VERSION)
	    if (y == 0.0) WARNING("Divide by zero");
	#endif
    
	    BTL_REAL y_reciprocal = 1.0 / y;
	    for (iterator i=vec.begin(); i!=vec.end(); i++)
	    	*i *= y_reciprocal;
	    return *this;
	}

//.............................................................................
    	    	/**#: [Description="Equality operator"]*/
	bool 
	operator==(const numeric_vector<T>& v) const
	{ 
    	
	    iterator i;
	    const_iterator j;
	    for (i=vec.begin(), j=v.begin(); i!=vec.end(); i++, j++)
            {
	    	if (*i != *j) return false;
            }
            return true;
	}

//.............................................................................
	    /**#: [Description="Postmultiplication of a Vector by a
	    	  BTL_Matrix. e.g. matrix m;
	    	  numeric_vector v1,v2; ... v2 = v1 * m; "]
         	  [Restrictions="The size of v must equal the number of
         	  rows in the Matrix."]*/

//    	friend numeric_vector<T>
//    	operator*(matrix<T>& m) // currently matrix not recognised as a type
//	{
//	#if defined(BTL_DEBUG_VERSION)
//	    if (vec.size() != m.num_rows())
//	     	FATAL_ERROR("The size of the input numeric_vector<T> must equal the number of rows"
//	     	            " in this Matrix");
//	#endif

//	    numeric_vector<T> result(vec.size());
//	    const_iterator i;
//	    iterator k;    
//	    const_iterator j = m.begin();
//	    for (i = vec.begin(), k = result.begin(); i < vec.end(); i++, k++)
//	        for (size_type l = 0; l < m.num_cols(); l++, j++)
//	            *k += *i * *j;
//	    return result;
//	}

};

_BTL_END_NAMESPACE

#endif // numeric_vector.h
