//
// This piece of code is a test program that uses the algorithms in 
// btl_sequence_algorithms.h, currently only a Pairwise Sequence Alignment. 
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
///////////////////////////////////////////////////////////////////////////
//
// Author: Breen Sweeney, Will Pitt, Mark Williams                
// 
///////////////////////////////////////////////////////////////////////////
// Modification History:
// Original version: B.S. (assistance from W.P)             Date: 30/08/98
// Last modified by: M.A.W.                                 Date: 31/08/99
// Comments:   
///////////////////////////////////////////////////////////////////////////
// 
// Brief Description of Code:
// 
// Carry out alignment of two sequences entered from keyboard or from a PIR 
// format file. User selectable gap penalty. With the supplied test.pir file
// this reproduces the results in fig. 2 of the original Needleman & Wunsch 
// (1970) J. Mol. Biol. 48, 443-453 publication
//
// For further details see the documentation for the matrix, vector classes 
// and the algorithms header files.
///////////////////////////////////////////////////////////////////////////////
// 

#include <vector>
#include <algorithm>
using namespace std;

#include "btl_biomolecular_data.h"
#include "btl_sequence_algorithms.h"
using namespace btl;


int main()
{
  	cout << "Welcome to the Pairwise Sequence Alignment program.\n";
  	cout << "We use the dynamic programming algorithm, as specified by\n";
  	cout << "Needleman and Wunsch.\n";

	// Declare some of the variables.
  	vector<char> sequence_s;
  	vector<char> sequence_t;

	// Here we decide whether to read from a file, or not.
  	cout << "\n\nEnter 1 if you want to enter your own short sequences \n";
	cout << "enter 2 to read the test PIR file : ";
  	int decide;
  	cin >> decide;
  	switch(decide) 
	{
  	    case 1: 
  	    {
  		 char seq_t[80], seq_s[80];
  		 cout << "Please enter the first sequence (S); \n";
  		 cin >> seq_s;
                 int i = 0;
 	         while (seq_s[i] != '\0')
                 {
	             sequence_s.push_back(seq_s[i]);
                     i++;
                 }
  		 cout << "Please enter the second sequence (T); \n";
  		 cin >> seq_t;
                 i = 0;
 	         while (seq_t[i] != '\0')
                 {
	             sequence_t.push_back(seq_t[i]);
                     i++;
                 }
  		 cout << "The sequences are now stored in a vector.\n";
  		 break;
	    }

  	    case 2: 
  	    {
    	     	 // Read the test file instead, and put into sequence_s, sequence_t.
    		 PIR_processor p;
    		 p.ReadFile("test.pir");
    		 if (p.empty())
    	        {
     		     cout << "Error, file empty.\n";
     		     exit(1);
    	    	 }
        	
    		 string seq1 = p.Seq(0);
    		 string seq2 = p.Seq(1);
    		 if (seq1.empty())
    		 {
     		     cout << "Error, sequence 1 empty.\n";
                   exit(1);
    		 }
    		 if (seq2.empty())
    		 {
     		     cout << "Error, sequence 2 empty.\n";
       	     exit(1);
    		 }

    		 // We now want to put the two strings into two vector sequences. 
    		 string::const_iterator i1 = seq1.begin();
    		 sequence_s.reserve(seq1.size());
    		 while (i1 != seq1.end()) sequence_s.push_back(*i1++); 
    		 string::const_iterator i2 = seq2.begin();
    		 sequence_t.reserve(seq2.size());
    		 while (i2 != seq2.end()) sequence_t.push_back(*i2++); 
  		 cout << "The sequences are now read from test.pir and stored in a vector.\n";
   		 break;
 	    }

  	    default:
               cout << "ERROR - input has gone wrong.......\n\n";
  	}

  	
  	// Now continue with program.
	cout << "This code implements a score of \n";
        cout << " (match_score or mismatch_score) - (gap_penalty + gap_length*gap_multiplier)\n";
        cout << "Needleman + Wunsch originally used a match score of 1 with all other terms zero.\n";

	cout << "\nPlease enter the gap_penalty : ";
  	double gap_penalty;
  	cin >> gap_penalty;
	cout << "\nPlease enter the gap length multiplier gap_length_multiplier r: ";
  	int gap_length_multiplier;
  	cin >> gap_length_multiplier;
	cout << "\nPlease enter the score for a match : ";
  	int match_score;
  	cin >> match_score;
	cout << "\nPlease enter the score for a mismatch : ";
  	int mismatch_score;
  	cin >> mismatch_score;

        // We have our standard template vectors, now get their sizes.
        int length_s = sequence_s.size();
  	int length_t = sequence_t.size();

        double score = 0.0;
        score = needleman_wunsch_similarity(sequence_s.begin(),sequence_s.end(),
                                    sequence_t.begin(),sequence_t.end(),
                                    match_score, mismatch_score, gap_penalty, 
                                    gap_length_multiplier,score);

	// Perform the alignment
        vector<char> alignment((2*max(length_s,length_t)));
        vector<char>::iterator alignend;
        alignend = needleman_wunsch_alignment(sequence_s.begin(),sequence_s.end(),
                               sequence_t.begin(),sequence_t.end(),
                               match_score, mismatch_score, gap_penalty, 
                               gap_length_multiplier, alignment.begin());

        cout << "\nThe similarity score is " << score << "\n";

	
  	// Display the aligned sequences.
  	int counter = 0;
        int maxcount = 0;
  	int display_length = 80;
        int aligned_length = alignend-alignment.begin() ;
  	cout << "\nThe alignment:\n";
        while(counter < aligned_length)
        {
            maxcount=min(aligned_length,counter+2*display_length);
 	    cout << "Sequence S:" << " ";
  	    for (int lj=counter; lj<maxcount; lj+=2) 
  	        cout << alignment[lj]; 
  	    cout << "\n";
  	    cout << "Sequence T:" << " ";
  	    for (int lj=(counter+1); lj<maxcount; lj+=2)
  	        cout << alignment[lj]; 
  	    cout << "\n\n";
            counter += 2*display_length;
         }

  	cout << "\n\n End program.\n";

        return 0;


} // End of main 
