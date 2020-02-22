//
// This file contains some example code, which may be used to call 
// the matrix & vector algorithms amd the least squares classes.
//
// This code is not necessarily meant to be useful in itself, but is provided as 
// an example of how the class may be used.
//
// Copyright (C) 1999 Software Engineering Group, Crystallography Department,
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
// Author: Mark Williams 
// 
/////////////////////////////////////////////////////////////////////////////////////
// 
// Brief Description of Code:
// 
// Superimpose two molecular structures in PDB format.
//
// For further details of the btl_least_squares algorithms, 
// see the documentation for the classes.
/////////////////////////////////////////////////////////////////////////////////////

// Standard header files
#include <vector>
#include <iostream>
using namespace std;

// BTL header files
#include "btl_biomolecular_data.h"
#include "btl_least_squares.h"
#include "btl_matrix.h"
#include "btl_numeric_vector.h"
#include "btl_matrix_algorithms.h"
using namespace btl;

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Usage: program_name firstPDBFile secondPDBFile" << endl;
        exit(1);
    }

    // Create objects to represent each structure using one of the file processor classes from the BTL
    // Read information from PDB files (reading only chains M and N, and the B atoms when alternatives are given)
    ATOM_processor A; A.ReadFile(argv[1],"MN ",'B');
    ATOM_processor B; B.ReadFile(argv[2],"MN ",'B');

    // The Coords() member function of ATOM_processor returns an STL vector containing the coordinates.
    // Consequently, the number of atoms in each file can be retrieved using the standard size() member function. 
    if (A.Coords().size() != B.Coords().size() ) {
        cerr << "Number of atoms unequal" << endl;
        exit(1);
    }
   
    bool long_way=false;
    if(long_way){ 
    // Do the superposition the long way in order to demonstrate the vector and matrix algorithms 
    // The geometric centre of each structure is declared as a BTL numeric_vector with 3 elements of
    // BTL_REAL(0.0) (the default). The coordinates of the centres are calculated using the generic 
    // BTL centroid algorithm is in this case operating on both STL and BTL vectors.
    numeric_vector<> centreA, centreB; 
    centroid(A.Coords().begin(), A.Coords().end(), centreA.begin());
    centroid(B.Coords().begin(), B.Coords().end(), centreB.begin());

    // Move protein A such that the protein centres are superimposed using the generic BTL algorithm `translate'
     numeric_vector<> translation = centreB - centreA;
     translate(A.Coords().begin(), A.Coords().end(), translation.begin());

    // Determine and perform the rotation necessary to superimpose structures
    // First calculate the Kearsley matrix and determine its eigenvalues and eigenvectors
    matrix<> matfit(4,4), evector(4,4); numeric_vector<> evalue(4);
    _kearsley_matrix(A.Coords().begin(), A.Coords().end(), B.Coords().begin(), B.Coords().end(), matfit.begin());
    
    eigen_solution(matfit.begin(), matfit.end(), 4 ,evector.begin(), evalue.begin());
    transpose(evector.begin(), evector.end(), 4, evector.begin());
    
    // Then rotate A about its centre in order to effect the superposition
    matrix<> rotation(3,3);
    rotation_from_fit(evector.begin(),rotation.begin());
    rotate(A.Coords().begin(), A.Coords().end(), rotation.begin(), centreB.begin());
    }

   else { 
   // Alternatively, and much shorter, the above steps are incorporated in a single algorithm in which the 
    // first protein's  coordinates are overwritten.  Here again we apply a BTL algorithm to the coordinate
    // data held in STL vectors. 
    double rmsd = 0.0;
    rmsd = lsqfit(A.Coords().begin(), A.Coords().end(), B.Coords().begin(), B.Coords().end(), rmsd);   
    cout << "Root mean square distance : " << rmsd << "\n";      
    }

    // The outstream operator << is overloaded to write the contents of an ATOM_processor object in PDB format.
    cout << A;		
    return 0;
}


