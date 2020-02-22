//
// This file contains some example code, which may be used to call 
// the sequence algorithms class.
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
// Find the number of elements that that are identical between two aligned sequences.
//
// For further details of the btl_sequence_algorithms, see the documentation for
// the class.
/////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
using namespace std;

#include "btl_sequence_algorithms.h"
using namespace btl;

int main() 
{

   // The sequences need to be held in containers that provide appropriate iterators 
   // for the BTL alignment algorithms. In this case, the sequences are entered from the keyboard and
   // the standard STL member function push_back() is used to place the data in STL vectors.

   vector<char> sequence, target; char temp[80];

   cout << "Please enter the first sequence (maximum 80 characters); \n";
   cin >> temp;
   int i = 0;
   while (temp[i] != '\0'){
       sequence.push_back(temp[i]);
        i++;
    }
    cout << "Please enter the second sequence (maximum 80 characters); \n";
    cin >> temp;
    i = 0;
    while (temp[i] != '\0'){
         target.push_back(temp[i]);
         i++;
    }

    // Calculate the similarity score for the two sequences using the default identity matrix scoring scheme in 
    // the BTL pairwise alignment algorithm. The iterators pointing to the beginning and end of each sequence 
    // are found using the standard begin() and end() STL container member functions.

    float score = 0.0;
    score = needleman_wunsch_similarity(sequence.begin(),sequence.end(),target.begin(),target.end(),1,0,0,0,score);
    cout << "\nThe number of identities between the aligned sequences is " << score << "\n";	

    return 0;
}

