//
// This file contains some example code, which may be used to call 
// the 3d coordinate fitting algorithms.
//
// This code is not necessarily meant to be useful in itself, but is provided as 
// an example of how the class may be used.
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
/////////////////////////////////////////////////////////////////////////////////////
//
// Author: Will Pitt, Breen Sweeney & Mary Steven
// 
/////////////////////////////////////////////////////////////////////////////////////
// 
// Brief Description of Code:
// 
// Fit a molecule to itself and to another of identical size.
//
// For further details see the documentation for the matrix, vector classes 
// and the algorithms header files.
/////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <complex>
#include <iomanip>
using namespace std;

#include "btl_numeric_vector.h"
#include "btl_matrix.h"
#include "btl_vector_algorithms.h"
#include "btl_least_squares.h"
using namespace btl;

//.............................................................................

void Check_lsqfit_rmsd(matrix<BTL_REAL>& m)
{
  cout<<"\n#####################################\n"
      << "### Check lsqfit_rmsd()          ####\n"
      << "###                              ####\n"
      << "### Fitting a matrix to itself   ####\n"
      << "#####################################\n";


  // Copy input to matrix a.
  matrix<BTL_REAL> a(m);

  // Define types to be used
  typedef matrix<BTL_REAL>::iterator iterator;
  typedef matrix<BTL_REAL>::const_iterator const_iterator;
  
  // (1)
  //
  // Calculate root mean squared deviation after a least squares fit of a to m
  // and output the result
    
  BTL_REAL ZERO = 0.0;
  BTL_REAL RMSD = lsqfit_rmsd(a.begin(), a.end(), m.begin(), m.end(), ZERO);
  
  cout << "\nRMSD's (should be close to zero)\n\n";
  cout << setprecision(16) << "After fitting to itself : " << RMSD << '\n';

  //    // (2)
  //    //
  //    // Translate input matrix then fit to the original.
  numeric_vector<BTL_REAL> trans(43.1,8.0,-9.0);
  a = m + trans;
    
  RMSD = lsqfit_rmsd(a.begin(), a.end(), m.begin(), m.end(), ZERO);
                                                           
  cout << "After translation, then fitting : " << RMSD << '\n';

  // (3)
  //
  // Rotate and translate input matrix then fit to the original.
  numeric_vector<BTL_REAL> axis(45.0,56.9,-9999.9);
  BTL_REAL angle = 1.0;  
  rotate(a.begin(),a.end(),axis.begin(),trans.begin(),angle);

  RMSD = lsqfit_rmsd(a.begin(), a.end(), m.begin(), m.end(), ZERO);

  cout << "After rotation and translation, then fitting : " << RMSD 
       << "\n\n";
         
}
//.............................................................................

void Check_lsqfit(matrix<BTL_REAL>& m)
{
  cout<<"\n####################################\n"
      << "### Check lsqfit()              ####\n"
      << "###                             ####\n"
      << "### Fitting a matrix to itself  ####\n"
      << "####################################\n";

  // Copy input to matrix a.
  matrix<BTL_REAL> a(m);

  // Define types to be used
  typedef matrix<BTL_REAL>::iterator iterator;
  typedef matrix<BTL_REAL>::const_iterator const_iterator;
    
  // (1)
  //
  // Translate input matrix then fit to the original.
  numeric_vector<BTL_REAL> trans(0.0,8.0,-9.0);
  a = m + trans;

  BTL_REAL RMSD = 0.0;
  RMSD = lsqfit(a.begin(), a.end(), m.begin(), m.end(), RMSD);

  cout << setprecision(16) << "\nRMSD's (should be close to zero)\n";
  cout << "\nRMSD after translation then fitting (Kearsley's method) : " 
       << RMSD;

  // Calculate RMSD using normal method

  RMSD = 0.0;
  RMSD = rmsd(a.begin(), a.end(), m.begin(), RMSD);

  cout << "\nRMSD after translation then fitting (explicit method)   : " 
       << RMSD;

  // (2)
  //
  // Rotate and translate input matrix then fit to the original.
  numeric_vector<BTL_REAL> axis(45.0,56.9,-99.9);
  BTL_REAL angle = 1.0;  
  rotate(a.begin(),a.end(),axis.begin(),trans.begin(),angle);
 
  RMSD = 0.0;
  RMSD = lsqfit(a.begin(), a.end(), m.begin(), m.end(), RMSD);
    
  cout << "\nRMSD after rotation, translation then fitting "
       << "(Kearsley's method) : " << RMSD;
  
  RMSD = 0.0;
  RMSD = rmsd(a.begin(), a.end(), m.begin(), RMSD);

  cout << "\nRMSD after rotation, translation then fitting "
       << "(explicit method)   : " << RMSD << "\n\n";
}

//.............................................................................

template<class iterator, class const_iterator>
void Check_lsqfit(iterator begin1, iterator end1, 
		  const_iterator begin2, const_iterator end2)
{
  cout<<"\n###############################\n"
      << "### Check lsqfit()         ####\n"
      << "### Fitting insulin (4ins) ####\n"
      << "### chain A to chain C     ####\n"
      << "###############################\n";

  // Do fit
  BTL_REAL RMSD = 0.0;
  RMSD = lsqfit(begin1, end1, begin2, end2, RMSD);

  cout << "\nRMSD (Kearsley's method) : "  << RMSD
       << " (Quanta and Profit 1.314)";

  // Check RMSD using normal method of calculation 
    
  RMSD = 0.0;
  RMSD = rmsd(begin1, end1, begin2, RMSD);

  cout << "\nRMSD (explicit method)   : "  << RMSD << "\n\n";
}
//.............................................................................

template<class iterator, class const_iterator>
void Check_lsqfit_rmsd(iterator begin1, iterator end1, 
		       const_iterator begin2, const_iterator end2)
{
  cout<<"\n###############################\n"
      << "### Check lsqfit_rmsd()    ####\n"
      << "### Fitting insulin (4ins) ####\n"
      << "### chain A to chain C     ####\n"
      << "###############################\n";

  // Do fit
  BTL_REAL RMSD = 0.0;
  RMSD = lsqfit_rmsd(begin1, end1, begin2, end2, RMSD);

  cout << "\nRMSD (Kearsley's method) : "  << RMSD 
       << " (Quanta and Profit 1.314)\n\n";

}

//.............................................................................

template<class iterator, class const_iterator>
void CheckRMSD(iterator begin1, iterator end1, 
	       const_iterator begin2, const_iterator end2)
{
  cout<<"\n###############################\n"
      << "### Check rmsd()           ####\n"
      << "### Using insulin (4ins)   ####\n"
      << "### chain A to chain C     ####\n"
      << "###############################\n";

  // Do calculation
  BTL_REAL RMSD = 0.0;
  RMSD = rmsd(begin1, end1, begin2, RMSD);

  cout << "\nRMSD (explicit method) : "  << RMSD 
       << " (Quanta 29.44 for 161 atoms)\n\n";
}

//.............................................................................

int main()
{

  // Set up data in the form of two 1D arrays - with pointers to the 
  // beginning and end of each array
  //
#include "test_p4ins.dat"
  BTL_REAL *beginA = p4ins_A, *endA = beginA + 163*3,
    *beginC = p4ins_C, *endC = beginC + 163*3;

  // (1)
  //
  // Test fitting routines: fitting Matrix to itself after rigid body
  // movements
  //
  matrix<BTL_REAL> matA(p4ins_A,163,3);
  matrix<BTL_REAL> matC(p4ins_C,163,3);

  Check_lsqfit_rmsd(matA);
  Check_lsqfit(matA);

  // (2)
  //
  // Check fitting of two sets of coordinates: in Matrix form
  //
  CheckRMSD(matA.begin(),matA.end(),matC.begin(),matC.end());
  Check_lsqfit_rmsd(matA.begin(),matA.end(),matC.begin(),matC.end());
  Check_lsqfit(matA.begin(),matA.end(),matC.begin(),matC.end());

  // (3)
  //
  // Test with STL vectors
  //
  vector<BTL_REAL> vecA(beginA,endA), vecC(beginC,endC);

  cout << "\n********* Repeat last 3 tests using STL vectors "
       << "*********\n\n";

  CheckRMSD(vecA.begin(),vecA.end(),vecC.begin(),vecC.end());
  Check_lsqfit_rmsd(vecA.begin(),vecA.end(),vecC.begin(),vecC.end());
  Check_lsqfit(vecA.begin(),vecA.end(),vecC.begin(),vecC.end());

  // (4)
  //
  // Test with 1D arrays instead of Matrix's
  //

  cout << "\n********* Repeat last 3 tests using arrays instead of Matrix's "
       << "*********\n\n";
  
  CheckRMSD(beginA,endA,beginC,endC);
  Check_lsqfit_rmsd(beginA,endA,beginC,endC);
  Check_lsqfit(beginA,endA,beginC,endC);

  return 0;
}

















