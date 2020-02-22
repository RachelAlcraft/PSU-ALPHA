//
// This file contains some example code, which may be used to call 
// the Fast Fourier Transform Class.
// This code is not meant to be usefull in itself, but is provided as 
// an example of how the FFT class may be used.
//
// Copyright (C) 1997-1999 Software Engineering Group, Crystallography Department,
// Birkbeck College, Malet Street, London WC1E 7HX, U.K.
// (d.moss@mail.cryst.bbk.ac.uk or m.williams@biochemistry.ucl.ac.uk)
// 
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Library General Public License as published by 
// the Free Software Foundation; either version 2 of the License, or (at your
// Handle) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////
//
// Author: Will Pitt, Mary Steven & Breen Sweeney, Mark Williams 
// 
///////////////////////////////////////////////////////////////////////////
// 
// Brief Description of Code:
// 
// This routine is a program which uses the Fast Fourier Transform class 
// (btl_fourier class) to carry out a One-dimensional, Two-dimensional, 
// and a Three-dimensional Fast Fourier Transform.
// This program is not meant to be used, and carries out no useful function
// itself. It is provided as an example of how the routines can be called
// from other C++ programs.
//
// For further details of the fourier_transform algorithms see the 
// documentation for the individual classes.
///////////////////////////////////////////////////////////////////////////

#include <string>
#include <complex>
#include <fstream>
using namespace std;

#include "btl_fourier.h"
#include "btl_numerics.h"
using namespace btl;


void TestResults(const string& testName, const complex_type& error)
{
    cout << " Ideally difference should be zero "
         << ">>>>> Difference = "
         << error
         << " : "
         << testName
         << " <<<<<\n"
         << endl;
}

//.............................................................................

template<class RandomAccessIterator>
void FillComplexContainer(RandomAccessIterator start,
                          const unsigned int size)
{
    // Set different random number seed each time this function is called.
    static unsigned int seed = 0;
    seed++;
    srand(seed);

    for (unsigned int i=0; i<size; i++)
    	start[i] = complex_type((rand()/10000.0),(rand()/10000.0));
}

//.............................................................................

template<class RandomAccessIterator >
complex<double> max_mismatch(RandomAccessIterator firsta,
                             RandomAccessIterator lasta,
                             RandomAccessIterator firstb,
                             complex<double> result)
{
    double maxError = 0.0, error = 0.0;
    for (; firsta<lasta; firsta++, firstb++) 
    {
         error = norm(*firsta - *firstb);
         if (error > maxError){
             result = *firsta - *firstb;
             maxError = error;
         }
    }           
    return result;
}

//..............................................................................

template <class RandomAccessIterator>
complex_type OneDimForwardBackwardTest(RandomAccessIterator testarr,
    	    	         	       RandomAccessIterator intermediate,
    	    	         	       RandomAccessIterator copy,
    	    	         	       const unsigned int size)
{    	    	     	 

    complex_type ZERO = 0.0;
    	// Transform in forward direction
    	//
    fourier_transform(testarr, intermediate, size, true);
   
    	// Transform back - overwriting the origin data
    	//
    fourier_transform(intermediate, testarr, size, false);
    
    	// Compare results with original
    	//
    return( max_mismatch(testarr,(testarr+size),copy,ZERO) );
}

//..............................................................................

template <class RandomAccessIterator>
complex_type TwoDimForwardBackwardTest(RandomAccessIterator testarr,
    	    	         	       RandomAccessIterator intermediate,
    	    	         	       RandomAccessIterator copy,
    	    	         	       const unsigned int nrows,
    	    	         	       const unsigned int ncols)
{    	    	     	 

    complex_type ZERO = 0.0;    
    	// Transform in forward direction
    	//
    fourier_transform(testarr, intermediate, nrows, ncols, true);
   
    	// Transform back - overwriting the origin data
    	//
    fourier_transform(intermediate, testarr, nrows, ncols, false);
    
    	// Compare results with original
    	//
    unsigned int size = nrows*ncols;    	
    return( max_mismatch(testarr,(testarr+size),copy,ZERO) );
}

//..............................................................................

template <class RandomAccessIterator>
complex_type ThreeDimForwardBackwardTest(RandomAccessIterator testarr,
    	    	         	         RandomAccessIterator intermediate,
    	    	         	         RandomAccessIterator copy,
    	    	         	         const unsigned int nrows,
    	    	         	         const unsigned int ncols,
    	    	         	         const unsigned int nlayers)
{    	    	     	 
    
    complex_type ZERO = 0.0;    
    	// Transform in forward direction
    	//
    fourier_transform(testarr, intermediate, nrows, ncols, nlayers, true);
   
    	// Transform back - overwriting the origin data
    	//
    fourier_transform(intermediate, testarr, nrows, ncols, nlayers, false);
    
    	// Compare results with original
    	//
    unsigned int size = nrows*ncols*nlayers;    	
    return( max_mismatch(testarr,(testarr+size),copy,ZERO)  );
}

//..............................................................................

void Check1DVectorInput(const unsigned int size)
{
   
    	// Create test data in array form
    	//
    vector<complex_type> testarr(size);
    FillComplexContainer(testarr.begin(),size);

    	// Copy original for comparison purposes
    	//
    vector<complex_type> copy = testarr;

    	// Create array to store intermediate results;
    	//
    vector<complex_type> intermediate(size);

    	// Do test
    	//
    complex_type error = OneDimForwardBackwardTest(testarr.begin(),
                                             	   intermediate.begin(),
                                             	   copy.begin(),
                                             	   size);
                                             
    TestResults("Results for vector of complex",error);
}

//..............................................................................

void Check2DVectorInput(const unsigned int nrows,
                            const unsigned int ncols)
{    
    	// Create test data in array form
    	//
    unsigned int size = nrows*ncols;
    vector<complex_type> testarr(size);
    FillComplexContainer(testarr.begin(),size);

    	// Copy original for comparison purposes
    	//
    vector<complex_type> copy = testarr;

    	// Create array to store intermediate results;
    	//
    vector<complex_type> intermediate(size);

    	// Do test
    	//
    complex_type error = TwoDimForwardBackwardTest(testarr.begin(),
                                             	   intermediate.begin(),
                                             	   copy.begin(),
                                             	   nrows,
                                             	   ncols);
                                             
    TestResults("Results for vector containing 2D data.",error);
}

//..............................................................................

void Check3DVectorInput(const unsigned int nrows,
                        const unsigned int ncols,
                        const unsigned int nlayers)
{
    	// Create test data in array form
    	//
    unsigned int size = nrows*ncols*nlayers;
    vector<complex_type> testarr(size);
    FillComplexContainer(testarr.begin(),size);

    	// Copy original for comparison purposes
    	//
    vector<complex_type> copy = testarr;

    	// Create array to store intermediate results;
    	//
    vector<complex_type> intermediate(size);

    	// Do test
    	//
    complex_type error = ThreeDimForwardBackwardTest(testarr.begin(),
                                             	     intermediate.begin(),
                                             	     copy.begin(),
                                             	     nrows,
                                             	     ncols,
                                             	     nlayers);
                                             
    TestResults("Results for vector containing 3D data.",error);
}


//..............................................................................

void Check1DArrayInput(const unsigned int size)
{
    	// Create test data in array form
    	//

    complex_type *testarr = new complex_type[size];
    FillComplexContainer(testarr,size);

    	// Copy original for comparison purposes
    	//
    complex_type *copy = new complex_type[size];
    for (unsigned int i=0; i<size; i++) copy[i] = testarr[i];

    	// Create array to store intermediate results;
    	//
    complex_type *intermediate = new complex_type[size];

    	// Do test
    	//
    complex_type error = OneDimForwardBackwardTest(testarr,
                                                   intermediate,
                                                   copy,
                                                   size);
                                                   
    TestResults("Results for array of complex",error);
     
    	// delete arrays
    	//
    delete[] testarr;
    delete[] copy;
    delete[] intermediate; 
}

//..............................................................................

void Check2DArrayInput(const unsigned int nrows, const unsigned int ncols)
{
    	// Create test data in array form
    	//
    unsigned int size = nrows*ncols;
    complex_type *testarr = new complex_type[size];
    FillComplexContainer(testarr,size);

    	// Copy original for comparison purposes
    	//
    complex_type *copy = new complex_type[size];
    for (unsigned int i=0; i<size; i++) copy[i] = testarr[i];

    	// Create array to store intermediate results;
    	//
    complex_type *intermediate = new complex_type[size];

    	// Do test
    	//
    complex_type error = TwoDimForwardBackwardTest(testarr,
                                                   intermediate,
                                                   copy,
                                                   nrows,
                                                   ncols);
                                                   
    TestResults("Results for array containing 2D data.",error);
     
    	// delete arrays
    	//
    delete[] testarr;
    delete[] copy;
    delete[] intermediate; 
}

//..............................................................................

void Check3DArrayInput(const unsigned int nrows, 
    	    	       const unsigned int ncols,
    	    	       const unsigned int nlayers)
{
    	// Create test data in array form
    	//
    unsigned int size = nrows*ncols*nlayers;
    complex_type *testarr = new complex_type[size];
    FillComplexContainer(testarr,size);

    	// Copy original for comparison purposes
    	//
    complex_type *copy = new complex_type[size];
    for (unsigned int i=0; i<size; i++) copy[i] = testarr[i];

    	// Create array to store intermediate results;
    	//
    complex_type *intermediate = new complex_type[size];

    	// Do test
    	//
    complex_type error = ThreeDimForwardBackwardTest(testarr,
                                                     intermediate,
                                                     copy,
                                                     nrows,
                                                     ncols,
                                                     nlayers);
                                                   
    TestResults("Results for array containing 3D data.",error);
     
    	// delete arrays
    	//
    delete[] testarr;
    delete[] copy;
    delete[] intermediate; 
}

//..............................................................................

void Series_1Dtest()
{

    typedef vector<complex_type> Container;
    unsigned int dataSize = 1024;
    Container sin_freq_10(dataSize),
    	      sin_freq_20(dataSize),
    	      sin_freq_30(dataSize),
    	      sin_freq_sum(dataSize);

    double step=0.0, stepSize = 2.0 * M_PI / dataSize;
    for (unsigned int i=0; i<dataSize; i++, step += stepSize)
    //for (int i=0; i<dataSize; i++, step += stepSize)
    {
    	sin_freq_10[i] = complex_type(sin(10.0*step));
    	sin_freq_20[i] = complex_type(sin(20.0*step));
    	sin_freq_30[i] = complex_type(sin(30.0*step));
    	sin_freq_sum[i] = sin_freq_10[i] + sin_freq_20[i] + sin_freq_30[i];
    }
    
    ofstream fin;
    fin.open("fft1D.in");

    // check that file can be opened
    if (!fin)
    {
        cerr << "\n\nError creating fft1D.in" << endl;
        return;
    }
    
    
    // Please note the ii,iii etc variables are for the benefit of the SGI CC
    // compiler 
    
    for (unsigned int ii=0; ii<dataSize; ii++)          // DEBUG was i<datasize
   
    //for (int ii=0; i<dataSize; ii++)
    	
    fin << real(sin_freq_sum[ii]) << '\n';
        
    Container output_10(dataSize),
    	      output_20(dataSize),
    	      output_30(dataSize),
    	      output_sum(dataSize);
    	      
    fourier_transform(sin_freq_10.begin(),output_10.begin(), dataSize, true);
    fourier_transform(sin_freq_20.begin(),output_20.begin(), dataSize, true);
    fourier_transform(sin_freq_30.begin(),output_30.begin(), dataSize, true);
    fourier_transform(sin_freq_sum.begin(),output_sum.begin(), dataSize, true);
   
    vector<double> powerSpectrum_sum(dataSize);
    Container FFT_sum(dataSize);
    double theoretical_result, diff, maxdiff = 0.0;

    for (unsigned int iii=0; iii<dataSize; iii++)
    {
    	FFT_sum[iii] = output_10[iii]+output_20[iii]+output_30[iii];
    	theoretical_result  = norm(FFT_sum[iii]);
    	powerSpectrum_sum[iii] = norm(output_sum[iii]);
    	diff = fabs(theoretical_result - powerSpectrum_sum[iii]);
    	if (diff > maxdiff) maxdiff = diff;
    }

    complex_type ZERO = 0.0;
    cout << "1D Sin series test: max difference from ideal result = " << maxdiff
         << "\n";
    
    ofstream fout;
    fout.open("fft1D.out");
    // check that file can be opened
    if (!fout.good())
    {
        cerr << "\n\nError creating fft1D.out" << endl;
        return;
    }
    
    // Perform a shift so that the origin lies in the centre.
    unsigned int j;
    for (unsigned int iiii=0; iiii<dataSize; iiii++)
    {
    	if (iiii < dataSize/2)
    	    j = dataSize/2 + iiii;
    	else
    	    j = iiii - dataSize/2;
    	fout << powerSpectrum_sum[j] << '\n';
    }    
    
    fout.close();       
}

//..............................................................................

void Series_2Dtest()
{

    // Create containers to contain test data
    //
    
    typedef vector< complex_type > Container;
    unsigned int dim = 64, nrows = dim, ncols = dim, dataSize = nrows*ncols;
    Container sin_freq_1(dataSize),
    	      sin_freq_2(dataSize),
    	      sin_freq_3(dataSize),
    	      sin_freq_sum(dataSize);
 
    // Put sin series with 4 different frequencies into containers
    //
    double rowstep, colstep, stepSize = 2.0 * M_PI / dim, start=-M_PI;
    unsigned int row,col,iarr;
    for (row=0, iarr=0, rowstep=start; row<nrows; row++, rowstep+=stepSize)
    	for (col=0, colstep=start; col<ncols; col++, colstep+=stepSize, iarr++)
    	{
    	      sin_freq_1[iarr] = complex_type(sin(1.0*rowstep)*sin(1.0*colstep));
    	      sin_freq_2[iarr] = complex_type(sin(2.0*rowstep)*sin(2.0*colstep));
    	      sin_freq_3[iarr] = complex_type(sin(3.0*rowstep)*sin(3.0*colstep));
    	    sin_freq_sum[iarr] = sin_freq_1[iarr] + 
    	                         sin_freq_2[iarr] +
    	                         sin_freq_3[iarr];
    	}
    
    // Output input data to a file
    //    
    ofstream fin;
    fin.open("fft2D.in");

    // check that file can be opened
    if (!fin)
    {
        cerr << "\n\nError creating fft2D.in" << endl;
        return;
    }
    for (row=0, iarr=0; row<nrows; row++)
    {
    	for (col=0; col<ncols; col++, iarr++)
    	    fin << real(sin_freq_sum[iarr]) << '\n';
    	fin << '\n';
    }
    
    // Create containers to contain Fourier transforms of the input data
    //
    Container output_1(dataSize),
    	      output_2(dataSize),
    	      output_3(dataSize),
    	      output_sum(dataSize);

    // Do the FFT calculations
    //    	      

    fourier_transform(sin_freq_1.begin(),output_1.begin(), nrows, ncols);
    fourier_transform(sin_freq_2.begin(),output_2.begin(), nrows, ncols);
    fourier_transform(sin_freq_3.begin(),output_3.begin(), nrows, ncols);
    fourier_transform(sin_freq_sum.begin(),output_sum.begin(), nrows, ncols);
   

    // Calculate the power spectrum of the resulting complex numbers to 
    // remove the imaginary part. Compare the theoretical result with the
    // calculated result
    //
    vector<double> powerSpectrum_sum(dataSize);
    double theoretical_result, diff, maxdiff = 0.0;

    for (unsigned int i=0; i<dataSize; i++)
    {
    	theoretical_result  = norm(output_1[i]+output_2[i]+output_3[i]);
    	                      
    	powerSpectrum_sum[i] = norm(output_sum[i]);
    	diff = fabs(theoretical_result - powerSpectrum_sum[i]);
    	if (diff > maxdiff) maxdiff = diff;
    }
    
    // Output the result to a file
    //
    cout << "2D Sin series test: max difference from ideal result = " << maxdiff << endl;
    
    ofstream fout;
    fout.open("fft2D.out");
    // check that file can be opened
    if (!fout.good())
    {
        cerr << "\n\nError creating fft2D.out" << endl;
        return;
    }
    
    // Perform a shift so that the origin lies in the centre.
    unsigned int x,y;
    for (row=0; row<nrows; row++)
    {
    	if (row < nrows/2)
    	    x = nrows/2 + row;
    	else
    	    x = row - nrows/2;
    	for (col=0; col<ncols; col++)
    	{
    	    if (col < ncols/2) 
    	    	y = ncols/2 + col;
    	    else
    	    	y = col - ncols/2;
    	    iarr = ncols*x + y; 
    	    fout << powerSpectrum_sum[iarr] << '\n';
    	}
    	fout << '\n';
    }    
}

//..............................................................................

template <class InputIterator>
void ArrayToFile3D(char *fileName, InputIterator begin, InputIterator end, 
                  unsigned int nrows,
                  unsigned int ncols,
                  unsigned int nslabs,
                  bool takeLogs=false)
{
    ofstream fout;
    fout.open(fileName);
    // check that file can be opened
    if (!fout.good())
    {
        cerr << "\n\nError creating " << fileName << endl;
        return;
    }
    
    unsigned int dataSize = nrows*ncols*nslabs;
    
    fout << "# vtk DataFile Version 1.0\n"
        << "title\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << nrows << " " << ncols << " "  << nslabs << '\n'
        << "ORIGIN 0.0 0.0 0.0\n"
        << "ASPECT_RATIO 1.0 1.0 1.0\n\n"
    	<< "POINT_DATA " << dataSize << '\n'
    	<< "SCALARS scalars float\n"
    	<< "LOOKUP_TABLE default\n";
    
    if (takeLogs)
    	for (unsigned int i=0; begin != end; begin++, i++) 
    	{
    	    fout << log(1+*begin) << " ";
    	    if (i==nslabs) 
    	    {
    	    	i=0;
    	    	fout << '\n';
    	    }
    	}
    else
    	for (unsigned int i=0; begin != end; begin++, i++) 
    	{
    	    fout << *begin << " ";
    	    if (i==nslabs) 
    	    {
    	    	i=0;
    	    	fout << '\n';
    	    }
    	}
}    
    
// Perform a shift so that the origin lies in the centre. 
//
template <class Iterator>
void Shift(Iterator begin, 
           Iterator tempBegin,
    	   unsigned int nrows,
           unsigned int ncols,
	   unsigned int nslabs)
{
    unsigned int dataSize = nrows*ncols*nslabs;

    // Copy data to temp container
    //
    Iterator i=begin, j = tempBegin;   
    for (unsigned int k=0; k<dataSize; k++)
    	*j++ = *i++;

    
    unsigned int row,col,slab,x,y,z,iarr,jarr;
    for (row=0; row<nrows; row++)
    {
    	if (row < nrows/2)
    	    x = nrows/2 + row;
    	else
    	    x = row - nrows/2;
    	    
    	for (col=0; col<ncols; col++)
    	{
    	    if (col < ncols/2)
    		y = ncols/2 + col;
    	    else
    		y = col - ncols/2;
    		
    	    for (slab=0; slab<nslabs; slab++)
    	    {
    	    	if (slab < nslabs/2)
    	    	    z = nslabs/2 + slab;
    	    	else
    	    	    z = slab - nslabs/2;

    	    	iarr = ncols*nslabs*row + nslabs*col + slab;
    	    	jarr = ncols*nslabs*x + nslabs*y + z;

    	    	begin[iarr] = tempBegin[jarr];
    	    }
    	}
    }   	   
}

void Series_3Dtest()
{

    // Create containers to contain test data
    //
    
    typedef vector< complex_type > Container;
    unsigned int dim = 32, nrows = dim, ncols = dim, nslabs = dim,
                 dataSize = nrows*ncols*nslabs;
    Container sin_freq_1(dataSize),
    	      sin_freq_2(dataSize),
    	      sin_freq_3(dataSize),
    	      sin_freq_sum(dataSize);
 
    // Put sin series with 4 different frequencies into containers
    //
    double rowstep, colstep, slabstep, stepSize = 2.0 * M_PI / dim, start=-M_PI;
    unsigned int row,col,slab,iarr;
    for (row=0, iarr=0, rowstep=start; row<nrows; row++, rowstep+=stepSize)
    	for (col=0, colstep=start; col<ncols; col++, colstep+=stepSize)
    	    for (slab=0, slabstep=start; slab<nslabs; slab++,
    	         slabstep+=stepSize, iarr++)
    	    {
    	      	sin_freq_1[iarr] =complex_type(sin(1*rowstep)*
                                              sin(1*colstep)*sin(1*slabstep));
    	      	sin_freq_2[iarr] =complex_type(sin(2*rowstep)*
                                              sin(2*colstep)*sin(2*slabstep));
    	      	sin_freq_3[iarr] =complex_type(sin(3*rowstep)*
                                              sin(3*colstep)*sin(3*slabstep));
    	    	sin_freq_sum[iarr] = sin_freq_1[iarr] + 
    	                             sin_freq_2[iarr] +
    	                             sin_freq_3[iarr];
    	    }
    	    
    // Output input data to a file
    //    
    ofstream fin;
    fin.open("fft3D_in.vtk");

    // check that file can be opened
    if (!fin)
    {
        cerr << "\n\nError creating fft3D_in.vtk" << endl;
        return;
    }
    fin << "# vtk DataFile Version 1.0\n"
        << "Sum of three sin functions of different frequencies "
        << "- test input to FFTTest.cpp\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << nrows << " " << ncols << " "  << nslabs << '\n'
        << "ORIGIN 0.0 0.0 0.0\n"
        << "ASPECT_RATIO 1.0 1.0 1.0\n\n"
    	<< "POINT_DATA " << dataSize << '\n'
    	<< "SCALARS scalars float\n"
    	<< "LOOKUP_TABLE default\n";
    for (row=0, iarr=0; row<nrows; row++)
    	for (col=0; col<ncols; col++)
    	{
    	    for (slab=0; slab<nslabs; slab++, iarr++)
    	    	fin << real(sin_freq_sum[iarr]) << " ";
    	    fin << '\n';
    	}
    
    
    
    // Create containers to contain Fourier transforms of the input data
    //
    Container output_1(dataSize),
    	      output_2(dataSize),
    	      output_3(dataSize),
    	      output_sum(dataSize);

    // Do the FFT calculations
    //    	      

    fourier_transform(sin_freq_1.begin(),output_1.begin(), nrows, ncols, nslabs);
    fourier_transform(sin_freq_2.begin(),output_2.begin(), nrows, ncols, nslabs);
    fourier_transform(sin_freq_3.begin(),output_3.begin(), nrows, ncols, nslabs);
    fourier_transform(sin_freq_sum.begin(),output_sum.begin(), nrows,ncols,nslabs);
   

    // Calculate the power spectrum of the resulting complex numbers to 
    // remove the imaginary part. Compare the theoretical result with the
    // calculated result
    //
    vector<double> powerSpectrum_sum(dataSize);
    double theoretical_result, diff, maxdiff = 0.0;

    for (unsigned int i=0; i<dataSize; i++)
    {
    	theoretical_result  = norm(output_1[i]+output_2[i]+output_3[i]);
    	powerSpectrum_sum[i] = norm(output_sum[i]);
    	diff = fabs(theoretical_result - powerSpectrum_sum[i]);
    	if (diff > maxdiff) maxdiff = diff;
    }    
    // Output the error to the screen
    //
    cout << "3D Sin series test: max difference from ideal result = " << maxdiff << endl;
    
    vector<double> temp(dataSize);
    Shift(powerSpectrum_sum.begin(), temp.begin(), nrows, ncols, nslabs);
    ArrayToFile3D("fft3D_out.vtk",
                  powerSpectrum_sum.begin(),
                  powerSpectrum_sum.end(),
                  nrows, ncols, nslabs, true);
    
}


//..............................................................................

void SinglePoint_3Dtest()
{

    // Create containers to contain test data
    //
    
    typedef vector< complex_type > Container;
    unsigned int dim = 32, nrows = dim, ncols = dim, nslabs = dim,
                 dataSize = nrows*ncols*nslabs;

    // Create container of all zero's
    //
    Container testData(dataSize, (complex_type) 0.0);
    Container testResults(dataSize);
    
    // Make first element equal to 1
    //
    testData[0] = (complex_type) 1;

    // Do the FFT calculations
    //    	      
    fourier_transform(testData.begin(),testResults.begin(), nrows, ncols, nslabs);
    
    // The resulting transform should be evenly spread out.
    //

    
    vector<double> powerSpec(dataSize);
    double normal = norm(testResults[0]), max=normal, min=normal;
    for (unsigned int i=0; i<dataSize; i++)
    {
    	normal  = norm(testResults[i]);
    	if (normal < min) min = normal;
    	if (normal > max) max = normal;
    	powerSpec[i] = normal;
    }

    vector<double> temp(dataSize);
    Shift(powerSpec.begin(), temp.begin(), nrows, ncols, nslabs);
    ArrayToFile3D("fft3D_point.vtk",
                  powerSpec.begin(),
                  powerSpec.end(),
                  nrows, ncols, nslabs, true);
    
    complex_type error = max-min;
    
    TestResults("3D Single point transform",error);
}

//..............................................................................

int main(int argc, char* argv[])
{

    SinglePoint_3Dtest();
    Series_1Dtest();
    Series_2Dtest();
    Series_3Dtest();
    
    unsigned int dim1 = 30, dim2 = 20, dim3 = 7;
    if (argc > 1) dim1 = atoi(argv[1]);
    if (argc > 2) dim2 = atoi(argv[2]);
    if (argc > 3) dim3 = atoi(argv[3]);
    
    cout << "\nSize of 1D containers: " << dim1 << "\n\n";
    
    Check1DVectorInput(dim1);
    Check1DArrayInput(dim1);

    cout << "\nSize of 2D containers: " << dim1 << " x " << dim2 << "\n\n";
    
    Check2DVectorInput(dim1,dim2);
    Check2DArrayInput(dim1,dim2);

    cout << "\nSize of 3D containers: " 
    	 << dim1 << " x " 
    	 << dim2 << " x " 
    	 << dim3 << "\n\n";
    
    Check3DVectorInput(dim1,dim2,dim3);
    Check3DArrayInput(dim1,dim2,dim3);

    return 0;
}
