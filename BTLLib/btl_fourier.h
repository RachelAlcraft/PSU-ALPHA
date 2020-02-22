//
// btl_fourier.h
//
// This file contains the header file for the Fast Fourier Transform Routine. 
// The routine allows one, two or three-dimensional transforms, in either 
// the 'forward' or 'reverse' directions.
//
// Copyright (C) 1997, 1998 Birkbeck College, Malet Street, London, U.K
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
//
// 

#if !defined(BTL_FOURIER_H)
#define BTL_FOURIER_H 1

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>

#include "BTL.h" 

 _BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description =" 
    This routine may be used to carry out a One-dimensional, Two-dimensional,
    or Three-dimensional Fast Fourier Transform. The basic algorithm has been
    adapted from freely-distributed code issued by the HP Reference Library
    It is based on the Cooley-Tukey algorihm (decimation in time), and
    it caters for sizes of the input array which are not a power of two.
    <P>
    The input parameters are;
    <P>
    First Input Parameter: A one, two or three-dimensional object, from
    the Fastft class.
    Second Input Parameter: A boolean variable. When set to 'true', this
    gives the 'forward' transform.
    <P>
    Note: in the comments following, the term 'array' is used to refer to the
    object created using the Fastft class.
    <P>
    Operation:  The routine carries out the standard Fast Fourier Transform. It
    checks to see if the size of the array is divisible by any of the first 12
    primes. If it is it continues using the usual algorithm. It copies z1 to z2
    as it proceeds, and the next time around uses z2 as the input array.
    Each iteration the variable 'arg' is set to 1, and 'omega' is set to the
    appropriate value. Arg is then multiplied by omega to give the relevant
    roots of unity. (For example, if n=8, in the first iteration -1 and 1 is
    used, then -i, 1, i, and -1, then in the last iteration the 8 '8th roots'
    are used).
    The routine should still work if the array size contains a factor which
    is not one of the first 12 primes. However it may be preferable to
    increase the number of primes in this case (in the 'prime' array), or
    alternatively pad the array with extra data points to ensure it is
    divisible by one of the 12 primes, or preferably a power of a prime.
    <P>
    Warnings:
    1. This routine does divide by n when doing the 'forward' transform. This
    is useful if the routine is used to transform forward and then back again.
    Most routines seem to leave this 'division by n' to the calling routine,
    however it is more sensible to let the calling routine 'multiply by n' if
    there is ever a need to do so.
    (E.g. (1, 0, 0, 0) will 'forward' transform to (.25, .25, .25, .25), and
    will then transform back to (1, 0, 0, 0).)
    <P>
    2. The calculation will leave small 'epsilon' values in the array
    where the value is calculated as zero. The calling routine should
    take whatever action is required.
    <P>
    Acknowledgements: This routine was originally adapted from code which
    has been circulated freely from the standard Hewlett-Packard Reference
    Library. It is essentially the same routine with the addition of the
    'esign' variable, which allows 'reverse' transforms as well as
    'forward' transforms; and with 'division by n' during the forward transform.
    The HP code was initially transcribed from Fortran as presented in
    'FFT as Nested Multiplication, with a Twist' by Carl de Boor in
    SIAM Sci. Stat. Comput., Vol 1 No 1, March 1980."]
    [Summary = "a class containing generic algorithms to calculate 1,2, and 3D
     fast Fourier transforms (FFT)."] 
    [Authors = "B.Sweeney, W.R.Pitt"]
    [Files = "<A HREF=./btl/btl_fourier.h>btl_fourier.h</A>"]
    [Dependencies="none"]
*/

static const double BTL_PI2 = 8.0 * atan(1.0);

// SGI CC version 7.1 does not implement the complex class as a template type.
#if !defined(SGI_CC)
    typedef complex<double> complex_type;
#else
    typedef complex complex_type;
#endif

//..................................................................................
// private function

template<class RandomAccessIterator>
void _Fftstp(const RandomAccessIterator zin,
    	      const unsigned int after,
    	      const unsigned int now,
    	      const unsigned int before,
    	      const RandomAccessIterator zout,
             const int esign)
{
    unsigned int now_x_after = now * after; // now*after is 'n'
     
    double angle = BTL_PI2/(now_x_after);	
    
    complex_type omega = complex_type(cos(angle), esign*sin(angle));
    
    complex_type arg = complex_type(1, 0);
    
    complex_type value;
    
    unsigned int j,ia,ib, before_x_after = before * after, offset;
    int in;
    
    for (j=0; j<now; j++)	    	// The 'j' gives us different values
    {				    	// for the zout index.
    	for (ia=0; ia<after; ia++)  	// The 'ia' and 'ib' indices
    	{			    	// give use the correct indices.
    	    for (ib=0; ib<before; ib++)
    	    {
    	    	offset = ia + ib*after;
     	    	value = zin[offset + (now-1) * before_x_after];
    		for (in=(now-2); 0<=in; in--)
    		{
    		    value *= arg;
    		    value += zin[offset + in * before_x_after];
    		}
    		zout[ia + j*after + ib * now_x_after] = value;
    	    }
    	    arg *= omega;   // gives use the appropriate root
       }			 //  of unity each loop.
    }
}       

//............................................................................
    	/**#: [Description="One dimensional FFT."] */
//
// N.B. prerequisites for type RandomAccessIterator are:
//
// 1. an operator[](int)
// 2. must contain type complex type
//
// N.B. start and result must point to separate areas of memory.

template <class RandomAccessIterator>
void fourier_transform(const RandomAccessIterator start, 
    	        	  const RandomAccessIterator result,
    	        	  const unsigned int size, 
    	        	  bool forward)
{
    if (start == result)
    {
    	cerr << "!!Error in Transform\n"
    	     << "start and result must point to separate areas of memory."
    	     << endl;
    	exit(1);
    }
    const unsigned int NEXTMX = 12;
    const int prime[NEXTMX] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
  
    // Here we initialise the 'sign', and the other initial variables...
    int esign = 1;
    if (forward) esign = -1;
    int inzee = 1;
    unsigned int before = size, after = 1,
         	 next = 0,
         	 now;

    do 
    {  	                      	    	// This do loop picks out the relevant
    	int np = prime[next]; 	    	// prime on each iteration.
 
    	if ( (before/np) * np < before )   // Does np divide before?
    	{
    	    if (++next < NEXTMX) continue;
    	    now = before;
    	    before = 1;
    	}
    	else
    	{
    	    now = np;
    	    before /= np;
    	}
    	if (inzee == 1)     	// Alternate between using z1 or z2 as input,
    	   _Fftstp(start, after, now, before, result, esign);
    	else
    	   _Fftstp(result, after, now, before, start, esign);

    	inzee = 3 - inzee;
    	after *= now;
    } while (1 < before);

    // Make sure the result contains the correct data.
    if (inzee == 1)
    {
    	// Copy data from start container to results container
    	RandomAccessIterator ir = result, is = start;
    	for (unsigned int i=0; i<size; i++) *ir++ = *is++; 
    }

    // Here we divide by the number of elements, if we are carrying
    // out the forward transformation.
    if (forward)
    	for (unsigned int i=0; i<size; i++)
    	    result[i] /= complex_type(size);
}

//............................................................................
	    	/**#: [Description="Two dimensional FFT."] */

template <class RandomAccessIterator>
void fourier_transform(const RandomAccessIterator start,
    	        	  const RandomAccessIterator result,
                	  const unsigned int nrows,
                	  const unsigned int ncols,
                	  const bool forward)
{
    if (start == result)
    {
    	cerr << "!!Error in Transform\n"
    	     << "start and result must point to separate areas of memory."
    	     << endl;
    	exit(1);
    }
    unsigned int size = nrows*ncols;
    unsigned int row, col;
    RandomAccessIterator is, ir;
    
    // We carry out the 1D transform on the rows and put answer in result
    // container.
    //
    for (row=0, is=start, ir=result; row<nrows; row++, is+=ncols, ir+=ncols)
    	fourier_transform(is, ir, ncols, forward);

    // We now copy the result to the start container, transposing the rows
    // and columns.
    //
    for (row=0, ir=result; row<nrows; row++) 
    	for (col=0; col<ncols; col++) 
    	    start[nrows*col + row] = *ir++;
   
    // We now transform it again (column wise with respect to original input).
    //
    for (col=0, is=start, ir=result; col<ncols; col++, is+=nrows, ir+=nrows)
    	fourier_transform(is, ir, nrows, forward);
   
    // We now copy back to start container, reversing the rows and columns.
    for (row=0, is=start; row<nrows; row++) 
    	for (col=0; col<ncols; col++) 
    	    *is++ = result[nrows*col + row];
    
    // Copy results for start to result container.
    //
    ir = result;
    is = start;
    for (unsigned int i=0; i<size; i++) *ir++ = *is++; 
}

//............................................................................
    	/**#: [Description="Three dimensional FFT."] */

template <class RandomAccessIterator>
void fourier_transform(const RandomAccessIterator start,
    	      	         const RandomAccessIterator result,
              	  const unsigned int nrows,
              	  const unsigned int ncols,
              	  const unsigned int nlayers,
              	  const bool forward)

{
    if (start == result)
    {
    	cerr << "!!Error in Transform\n"
    	     << "start and result must point to separate areas of memory."
    	     << endl;
    	exit(1);
    }
    
    unsigned int row, col, layer;
    RandomAccessIterator is, ir;
    unsigned int layerSize=nrows*ncols;
    
    // We carry out the 2D transform on the start container putting the answer
    // in the result container.
    //
    for (layer=0, is=start, ir=result; layer<nlayers; layer++,is+=layerSize,ir+=layerSize)
    	fourier_transform(is, ir, nrows, ncols, forward);

    // Perform 1D transform on each row, column combination
    //
    vector<complex_type> tempStart(nlayers), tempResult(nlayers);
    unsigned int row_offset, in_layer_offset;
    for (row=0; row<nrows; row++)
    {
    	row_offset = ncols * row;
    	for (col=0; col<ncols; col++)
    	{
    	    // Copy row,col combination to single continous section of memory
    	    //
    	    in_layer_offset = row_offset + col;
    	    for (layer=0; layer<nlayers; layer++)
    	    	tempStart[layer] = result[layerSize*layer + in_layer_offset];
    	    	
    	    // Carry out 1D FFT
    	    //
    	    fourier_transform(tempStart.begin(), tempResult.begin(), nlayers, forward);
    	    
    	    // Copy result back to original locations
    	    //
    	    for (layer=0; layer<nlayers; layer++)
    	    	result[layerSize*layer + in_layer_offset] = tempResult[layer];
    	}
    }
}

_BTL_END_NAMESPACE

#endif
