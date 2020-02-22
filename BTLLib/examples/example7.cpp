//
// This piece of code is a test program that uses the random number generators
// in btl_numerics.h.
//
// Copyright (C) 1998,1999 Birkbeck College, Malet Street, London WC1E 7HX, U.K.
// d.moss@mail.cryst.bbk.ac.uk or m.williams@biochemistry.ucl.ac.uk
// 
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////////
//
// Author: Mark Williams                
// 
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 
// Brief Description of Code:
// 
// Implements the rudimentary tests given in Section 3.6 (and associated exercises)
// of D.E. Knuth's the Art of Computer Programming 3rd Edition. 
//
///////////////////////////////////////////////////////////////////////////////
// 
 
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
using namespace std;

#include "btl_numerics.h"
using namespace btl;

int main()
{
	register int m,j;

        // Test the generation of floats (doubles would use same code code internally) 

        double a[2009]; 
	clock_t initialize = clock();
	clock_t generated = clock();
	clock_t start = clock();
	double speed = CLOCKS_PER_SEC;
	random_uniform<float> tempran(310952);
	initialize = clock();
        for(m=0;m<=2009;m++){
            for(j=0;j<1009;j++)  a[j] = tempran();   
        }
	generated = clock();
        cout << "Testing the BTL random number generator \n\n";
	cout << "Initialisation  " << double(initialize - start)/speed << " seconds \n";
	cout << "Generation of 2 million seperate numbers " << double(generated - start)/speed << " seconds \n";
        cout << " Floating point test result should equal 0.274526 \n";
	cout << "Actual result " << a[0] << " \n\n";

        // Test the generation of 32 bit integers 

	long la[2009];
	start = clock();
	random_uniform<long> templran(310952);
	initialize = clock();
        for(m=0;m<=2009;m++){
            for(j=0;j<1009;j++)  la[j] = templran();   
        }
	generated = clock();
	cout << "Initialisation  " << double(initialize - start)/speed << " seconds \n";
	cout << "Generation of 2 million seperate numbers " << double(generated - start)/speed << " seconds \n";
        cout << "Long integer test result should be 461390032 \n";
	cout << "Actual result " << la[0] << " \n\n";

        // Test the generation of 16 bit integers

	int ia[2009];
	start = clock();
	random_uniform<int> tempiran(12509);
	initialize = clock();
        for(m=0;m<=2009;m++){
        for(j=0;j<1009;j++)  ia[j] = tempiran();   
        }
	generated = clock();
	cout << "Initialisation  " << double(initialize - start)/speed << " seconds \n";
	cout << "Generation of 2 million seperate numbers " << double(generated - start)/speed << " seconds \n";
        cout << "Integer test result should be 9387 \n";
	cout << "Actual result " << ia[0] << " \n\n";

        // The normal distribution
	double na[2009];
	start = clock();
	random_normal<double> tempnran(310952);
	initialize = clock();
        for(m=0;m<=2009;m++){
        for(j=0;j<1009;j++)  na[j] = tempnran();   
        }
	generated = clock();
	cout << "Initialisation  " << double(initialize - start)/speed << " seconds \n";
	cout << "Generation of 2 million normally distributed numbers " << double(generated - start)/speed << " seconds \n\n";

        // The exponential distribution
	double ea[2009];
	start = clock();
	random_exponential<double> temperan(310952);
	initialize = clock();
        for(m=0;m<=2009;m++){
        for(j=0;j<1009;j++)  ea[j] = temperan();   
        }
	generated = clock();
	cout << "Initialisation  " << double(initialize - start)/speed << " seconds \n";
	cout << "Generation of 2 million exponentially distributed numbers " << double(generated - start)/speed << " seconds \n\n";


	return 0;
}

