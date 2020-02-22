//
// btl_numerics.h
//
// This file contains generic statistical functions 
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

#if !defined (BTL_NUMERICS_H)
#define BTL_NUMERICS_H 1

/**#: [Description ="A collection of generic numerical functions."]
    [Summary = "generic numerical algorithms"] 
    [Authors = "M.A.Williams"]
    [Files = "<A HREF=./btl/btl_numerics.h>btl_numerics.h</A>"]
    [Dependencies="None"]
    [Prerequisites="None<P>"
    <P>
*/

   
#include "BTL.h"

#include <cmath>
#include <functional>

    	
_BTL_BEGIN_NAMESPACE

using namespace std;

//..........................................................................................
// NUMERIC COMPARISON
//..........................................................................................

    	    /**#: [Description="Returns the maximum difference between
    		   equivalent elements in two containers in the sense of (a-b)"]*/

        template <class ForwardIterator1, class ForwardIterator2, class T>
        T max_mismatch(ForwardIterator1 firsta, ForwardIterator1 lasta, 
                                ForwardIterator2 firstb, T init )
        {
    	    for( ; firsta < lasta; firsta++, firstb++)
            {
               init = max(init,(*firsta-*firstb));
            }
            return init; 
        }
//..........................................................................................

    	    /**#: [Description="Returns the largest absolute difference between
    		   equivalent elements in two containers"]*/

        template <class ForwardIterator1, class ForwardIterator2, class T>
        T max_absolute_mismatch(ForwardIterator1 firsta, ForwardIterator1 lasta, 
                                ForwardIterator2 firstb, T init )
        {
    	    for( ; firsta < lasta; firsta++, firstb++)
            {
               init = max(init,fabs(*firsta-*firstb));
            }
            return init; 
        }

//..........................................................................................

    	    /**#: [Description="Returns the minimum difference between
    		   equivalent elements in two containers in the sense of (a-b)"]*/

        template <class ForwardIterator1, class ForwardIterator2, class T>
        T min_mismatch(ForwardIterator1 firsta, ForwardIterator1 lasta, 
                                ForwardIterator2 firstb, T init )
        {
    	    for( ; firsta < lasta; firsta++, firstb++)
            {
               init = min(init,(*firsta-*firstb));
            }
            return init; 
        }

//..........................................................................................

    	    /**#: [Description="Returns the smallest absolute difference between
    		   equivalent elements in two containers"]*/

        template <class ForwardIterator1, class ForwardIterator2, class T>
        T min_absolute_mismatch(ForwardIterator1 firsta, ForwardIterator1 lasta, 
                                ForwardIterator2 firstb, T init )
        {
    	    for( ; firsta < lasta; firsta++, firstb++)
            {
               init = min(init,fabs(*firsta-*firstb));
            }
            return init; 
        }

//..........................................................................................


    	    /**#: [Description="Returns true if every element of a is less than 
                                or equal to the equivalent element of b"]*/

        template <class ForwardIterator1, class ForwardIterator2>
        bool every_less_equal(ForwardIterator1 firsta, ForwardIterator1 lasta, 
                              ForwardIterator2 firstb)
        {
            bool result=true;
    	    for( ; firsta < lasta; firsta++, firstb++)
            {
               if(*firstb < *firsta)
                  return result=false;
            }
            return result; 
        }

//..........................................................................................


    	    /**#: [Description="Returns true if every element of a is less than 
                                the equivalent element of b"]*/

        template <class ForwardIterator1, class ForwardIterator2>
        bool every_less(ForwardIterator1 firsta, ForwardIterator1 lasta, 
                              ForwardIterator2 firstb)
        {
            bool result=true;
    	    for( ; firsta < lasta; firsta++, firstb++)
            {
               if(*firstb <= *firsta)
                  return result=false;
            }
            return result; 
        }

//..........................................................................................
	    	/**#: [Description="Compare the absolute values of two arguements
			 	    true if |x| < |y|"] */

        template<class T> 
        struct less_absolute : public binary_function<T,T,bool>
        {
            bool operator() (const T& x, const T& y) const 
            {
	    return (x*x) < (y*y);
            }
	};

//............................................................................................

	    	/**#: [Description="Portable random number generator that produces
                                    float,double, int and long random numbers.
                                    In the case of float & double the numbers are uniformly 
                                    distributed in the range 0 to 1. In the case of long they are
                                    uniformly distributed in the range 0 to 2**31 and
                                    in the int case from 0 to 2**15.
                                    It is an implementation of Donald Knuth's floating point 
                                    ranf_array generator, a lagged additive generator
                                    x(i) = (x(i-KK) + x(i-LL))mod 1.0 for float and double numbers
                                    and Knuth's lagged subtractive generator 
                                    x(i) = (x(i-KK) - x(i-LL))mod 2**MM where MM = 15
                                    for int and MM=30 for long.
                                    An object containing random numbers is constructed and used by
                                    
                                    random_uniform<T> numbers(seed);
                                    T temp = numbers();
                                    
                                    The seed integer should be < 1073741822 for float, double and long
                                    numbers and < 32765 in the case of of int data. 
                                    N.B. it is best not to generate a sequence of more than 
                                    2**70 numbers using the same seed. This generator fails 
                                    the birthday spacings test i.e some differences between pairs
                                    of numbers occur more often than desired. 
                                    See Knuth's Art of Computer Programming 3rd Edition for solutions"] */

// This implementation reproduces the results of the rudimentary test given in Section 3.6
// of the Art of Computer Programming 3rd Edition and differs from the code given there in using
// the C++ standard math functions ldexp and modf, and no bitwise operations.


template <class T = BTL_REAL>
class random_uniform
{
private:
	static const int KK = 100;       // the longer lag
	static const int LL = 37;        // the short lag
                                         // lag pairs (127, 30) and (258, 83) can also be used
	static const int TT = 70;        // guarantees many seperate streams of numbers 
					 // ~ 2^TT streams when TT  =< KK - 1 ?
	static const int DK = 200;       // (KK + KK)
	static const int KL = 63;        // (KK - LL)
	double ran_u[KK];
        int current_position;
        static const int BUFFER_SIZE = 1009;
        double buffer[BUFFER_SIZE];

void fill_buffer(int n)		// fill an array->vector size n
{
    	register int i,j;
    	double a; double* pd = &a;
    	for (j=0;j<KK;j++)     buffer[j] = ran_u[j];
    	for (;j<n;j++)         buffer[j] = modf((buffer[j-KK] + buffer[j-LL]),pd); 
    	for (i=0;i<LL;i++,j++)  ran_u[i] = modf((buffer[j-KK] + buffer[j-LL]),pd); 
    	for (;i<KK;i++,j++)     ran_u[i] = modf((buffer[j-KK] + ran_u[i-LL]),pd);
        current_position = 0;
}

protected:

T next_number() 
{
        if(current_position == BUFFER_SIZE) fill_buffer(BUFFER_SIZE);
        return buffer[current_position++];
}


public:

random_uniform(long seed)
{ 
    	register int s,t,j;
    	double a; double* pd = &a;
    	double  u[DK-1]; 
    	double ul[DK-1];
    	double ulp = ldexp(1.0,-52);  			// 2 to the -52
    	double ss = 2.0*ulp*((seed&0x3fffffff)+2);      // the & construction zeros the two highest bits of the seed
    	for (j=0;j<KK;j++)
    	{
            u[j] = ss;
            ul[j] = 0.0;                         	// bootstrap the buffer
            ss += ss;
            if(ss>=1.0) ss -= (1.0 - 2*ulp);            // cyclic shift of 51 bits
    	}
    	for(;j<DK-1;j++) u[j] = ul[j] = 0.0;
    	u[1] += ulp;
    	ul[1] = ulp;                                    // make u[1] and only u[1] odd
    	s = seed&0x3fffffff;				// s and t are counters
    	t = TT - 1;

    	while (t)
    	{ 
            for (j=KK-1;j>0;j--)      ul[j+j] = ul[j],    u[j+j] = u[j];
            for (j=DK-2;j>KL;j-=2)    ul[DK-1-j] = 0.0,   u[DK-1-j] = u[j]-ul[j];
            for (j=DK-2;j>=KK;j--)
            {
                if (ul[j])
                {
                    ul[j-KL] = ulp - ul[j-KL], u[j-KL] = modf((u[j-KL] + u[j]),pd);
                    ul[j-KK] = ulp - ul[j-KK], u[j-KK] = modf((u[j-KK] + u[j]),pd);
                }
	    }
            if ((s%2)==1)                               // if s is odd
            {
                for (j=KK;j>0;j--) ul[j] = ul[j-1], u[j] = u[j-1];
                ul[0] = ul[KK], u[0] = u[KK];
                if (ul[KK]) ul[LL] = ulp - ul[LL], u[LL] = modf((u[LL] + u[KK]),pd);
            }

            if (s) s /= 2;          // successive integer divides of seed value by 2 
	    else t--;		    // when 0 is reached countdown TT further cycles
        }
 
        for(j=0;j<LL;j++) ran_u[j+KL] = u[j];
        for(;j<KK;j++)    ran_u[j-LL] = u[j];
        fill_buffer(BUFFER_SIZE);
}

T operator() () 
{
        if(current_position == BUFFER_SIZE) fill_buffer(BUFFER_SIZE);
        return buffer[current_position++];
}

}; // end of rng for floating point numbers


//............................................................................................

template <>
class random_uniform<long>
{
private:

// The lagged subtractive generator x(i) = (x(i-KK) - x(i-LL))mod 2**MM 

	static const long MM = 1L<<30;   // modulus = 2**30
	static const int KK = 100;       // the long lag
	static const int LL = 37;        // the short lag
                                         // lag pairs (127, 30) and (258, 83) can also be used
	static const int TT = 70;        // guarantees many seperate streams of numbers 
					 // ~ 2**TT streams when TT + MM = KK
	static const int DK = 200;       // KK+ KK
	static const int KL = 63;        // KK - LL

	long ran_x[KK];		         // current state of random number generator

        int current_position;
        static const int BUFFER_SIZE = 1009;
        long buffer[BUFFER_SIZE];

void fill_buffer(int n)
{
	register int i,j;
	for (j=0;j<KK;j++)     buffer[j] = ran_x[j];
	for (;j<n;j++)         buffer[j] = ((buffer[j-KK] - buffer[j-LL])&(MM-1));  // sum mod MM-1 
	for (i=0;i<LL;i++,j++)  ran_x[i] = ((buffer[j-KK] - buffer[j-LL])&(MM-1));
	for (;i<KK;i++,j++)     ran_x[i] = ((buffer[j-KK] - ran_x[i-LL])&(MM-1));
        current_position = 0;
}

public:

random_uniform(long seed)
{ 
// add debug for n > 1073741821 i.e. (MM-3)
    register int t,j;
    long  x[DK-1]; 
    register long ss = ((seed+2)&(MM-2));

    for (j=0;j<KK;j++)
    {
        x[j] = ss;
        ss *= 2;
        if(ss>=MM) ss -= (MM-2);                
    }
    for(;j<DK-1;j++) x[j] = 0;
    x[1]++;
    ss = seed&(MM-1);                  // ss forced < MM, both s and t are counters
    t = TT - 1;

    while (t)
    { 
        for (j=KK-1;j>0;j--)  x[j+j] = x[j];
        for (j=DK-2;j>KK-LL;j-=2) x[DK-1-j] = (x[j]&(MM-2));
        for (j=DK-2;j>=KK;j--)
        {
            if ((x[j]%2) == 1)
            {
                x[j-KL] = ((x[j-KL] - x[j])&(MM-1));
                x[j-KK] = ((x[j-KK] - x[j])&(MM-1));
            }
	}
        if ((ss%2) == 1)                               // if s is odd
        {
            for (j=KK;j>0;j--) x[j] = x[j-1];
            x[0] = x[KK];
            if ((x[KK]%2) == 1) x[LL] = ((x[LL] - x[KK])&(MM-1)) ; //???
        }

         if (ss) ss /= 2;             // successive integer divides of seed value by 2 
	 else t--;		      // when 0 is reached countdown TT further cycles
     }
 
     for(j=0;j<LL;j++) ran_x[j+KL] = x[j];
     for(;j<KK;j++)    ran_x[j-LL] = x[j];
     fill_buffer(BUFFER_SIZE);
}


long operator() () 
{
        if(current_position == BUFFER_SIZE) fill_buffer(BUFFER_SIZE);
        return buffer[current_position++];
}

}; // end of rng for long integers


//............................................................................................

template <>
class random_uniform<int>
{
private:

// The lagged subtractive generator x(i) = (x(i-KK) - x(i-LL))mod 2**MM 

	static const unsigned int MM = 1U<<15;   // modulus = 2**15
	static const int KK = 100;       // the long lag
	static const int LL = 37;        // the short lag
                                         // lag pairs (127, 30) and (258, 83) can also be used
	static const int TT = 70;        // guarantees many seperate streams of numbers 
					 // ~ 2**TT streams when TT + MM = KK
	static const int DK = 200;       // KK + KK
	static const int KL = 63;        // KK - LL

	int ran_x[KK];		         // current state of random number generator

        int current_position;
        static const int BUFFER_SIZE = 1009;
        int buffer[BUFFER_SIZE];

void fill_buffer(int n)
{
	register int i,j;
	for (j=0;j<KK;j++)     buffer[j] = ran_x[j];
	for (;j<n;j++)         buffer[j] = ((buffer[j-KK] - buffer[j-LL])&(MM-1));  // sum mod MM-1 
	for (i=0;i<LL;i++,j++)  ran_x[i] = ((buffer[j-KK] - buffer[j-LL])&(MM-1));
	for (;i<KK;i++,j++)     ran_x[i] = ((buffer[j-KK] - ran_x[i-LL])&(MM-1));
        current_position = 0;
}

public:

random_uniform(int seed)
{ 
    register int t,j;
    int x[DK-1]; 
    register unsigned int ss = ((seed+2)&(MM-2));

    for (j=0;j<KK;j++)
    {
        x[j] = ss;
        ss *= 2;
        if(ss>=MM) ss -= (MM-2);                
    }
    for(;j<DK-1;j++) x[j] = 0;
    x[1]++;
    ss = seed&(MM-1);                  // ss forced < MM, both s and t are counters
    t = TT - 1;

    while (t)
    { 
        for (j=KK-1;j>0;j--)  x[j+j] = x[j];
        for (j=DK-2;j>KK-LL;j-=2) x[DK-1-j] = (x[j]&(MM-2));
        for (j=DK-2;j>=KK;j--)
        {
            if ((x[j]%2) == 1)
            {
                x[j-KL] = ((x[j-KL] - x[j])&(MM-1));
                x[j-KK] = ((x[j-KK] - x[j])&(MM-1));
            }
	}
        if ((ss%2) == 1)                               // if s is odd
        {
            for (j=KK;j>0;j--) x[j] = x[j-1];
            x[0] = x[KK];
            if ((x[KK]%2) == 1) x[LL] = ((x[LL] - x[KK])&(MM-1)) ; //???
        }

         if (ss) ss /= 2;             // successive integer divides of seed value by 2 
	 else t--;		      // when 0 is reached countdown TT further cycles
     }
 
     for(j=0;j<LL;j++) ran_x[j+KL] = x[j];
     for(;j<KK;j++)    ran_x[j-LL] = x[j];
     fill_buffer(BUFFER_SIZE);
}


int operator() () 
{
        if(current_position == BUFFER_SIZE) fill_buffer(BUFFER_SIZE);
        return buffer[current_position++];
}

}; // end of rng for integers

//...................................................................................

	    	/**#: [Description="Portable random number generator that produces float
                                    or double precision random numbers from a normal distibution
                                    using the ratio method of Kinderman & Monahan (see Knuth TAOCP)."] */

template <class T = BTL_REAL>
class random_normal : public random_uniform<T>
{
private:

        double m, sd, mult1, mult2, temp;
	T u,v;

public:

random_normal(long seed) : random_uniform<T>(seed)
{
        m = 0.0; sd = 1.0;
        mult1 = -0.5*exp(1.0); mult2 = sqrt(8.0*exp(-1.0))*sd;
}
random_normal(long seed, T mean, T standard_deviation) : random_uniform<T>(seed) 
{
        m = mean; sd = standard_deviation;
        mult1 = -0.5*exp(1.0); mult2 = sqrt(8.0*exp(-1.0)) * sd;
}

T operator() () 
{
        u = random_uniform<T>::next_number();
        v = random_uniform<T>::next_number() - 0.5;
        while (v*v > mult1*u*u*log(u))
        {
            u = random_uniform<T>::next_number();
            v = random_uniform<T>::next_number() - 0.5;
        }
        return m + mult2*v/u;
}

}; // end of rng for normally distributed floating point numbers
 

//........................................................................................


	    	/**#: [Description="Portable random number generator that produces float
                                    or double precision random numbers from an exponential distibution."] */

template <class T = BTL_REAL>
class random_exponential : public random_uniform<T>
{
private:

	double m;

public:

random_exponential(long seed) : random_uniform<T>(seed)
{
	m = 1.0; 
}
random_exponential(long seed, T mean) : random_uniform<T>(seed)
{
	m = mean;
}

T operator() () 
{
	return m * log(random_uniform<T>::next_number()); 
}

}; // end of rng for exponentially distributed floating point numbers
 
//.........................................................................................

_BTL_END_NAMESPACE

#endif

