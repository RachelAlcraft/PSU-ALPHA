//
// pdb_sort.h
//
// This file contains declarations for the pdb_sort class.
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London, U.K.
//
// This library is free software; you can
// redistribute it and/or modify it under the terms of
// the GNU Library General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined(BTL_PDBSORT_H)
#define BTL_PDBSORT_H 1


#include <string>
#include <vector>

#include "BTL.h"

_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description ="
     This class defines a functor or function object which can be used by
     a generic sort routine to sort atoms into the order used in 
     Brookhaven pdb format files. Order is defined by strings, 2 or 3 
     characters long."]
      [Summary = "function object (functor) used for sorting objects into
      Brookhaven Protein Data Bank (pdb) format atom order"] 
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/pdb_sort.h>pdb_sort.h</A>"]
    [Dependencies="none"]

*/
    static const short  ATOMMAX = 68;

    static const string ATOMNAME[68]= { 

    "N",     //  0
    "CA",    //  1
    "C",     //  2
    "O",     //  3
    "CB",    //  4
    "SG",    //  5
    "OG",    //  6
    "CG",    //  7
    "OG1",   //  8
    "CG1",   //  9
    "CG2",   // 10
    "CD",    // 11
    "OD",    // 12
    "SD",    // 13
    "CD1",   // 14
    "OD1",   // 15
    "ND1",   // 16
    "CD2",   // 17
    "OD2",   // 18
    "ND2",   // 19
    "CE",    // 20
    "NE",    // 21
    "CE1",   // 22
    "OE1",   // 23
    "NE1",   // 24
    "CE2",   // 25
    "OE2",   // 26
    "NE2",   // 27
    "CE3",   // 28
    "CZ",    // 29
    "NZ",    // 30
    "CZ2",   // 31
    "CZ3",   // 32
    "OH",    // 33
    "NH1",   // 34
    "NH2",   // 35
    "CH2",   // 36
    "OXT",   // 37
    "P",     // 38
    "O1P",   // 39
    "O2P",   // 40
    "O5*",   // 41
    "C5*",   // 42
    "C4*",   // 43
    "O4*",   // 44
    "C3*",   // 45
    "O3*",   // 46
    "C2*",   // 47
    "O2*",   // 48
    "C1*",   // 49
    "N9",    // 50
    "C8",    // 51
    "N7",    // 52
    "C5",    // 53
    "C6",    // 54
    "O6",    // 55
    "N6",    // 56
    "N1",    // 57
    "C2",    // 58
    "O2",    // 59
    "N2",    // 60
    "N3",    // 61
    "C4",    // 62
    "O4",    // 63
    "N4",    // 64
    "C5",    // 65
    "C5M",   // 66
    "C6" };  // 67

class pdb_sort
{

private:

	    	/**#: [Hidden] */
    bool    
    EqualsAtomName(const string atomId, const short nameIdx) const
    {
        if (atomId == ATOMNAME[nameIdx] ) return true;
        if (nameIdx > 40 && nameIdx < 50 && 
	atomId.size() == 3 &&
	atomId[0] == ATOMNAME[nameIdx][0] &&
	atomId[1]  == ATOMNAME[nameIdx][0]) return true;
    return false;
    }
    	    
	    	/**#: [Hidden] */
    short   
    FindIndx(const string& atomId1) const
    {
    // Remove all blanks
    string s1 = "";
        for(unsigned short i = 0; i < atomId1.size(); i++)
        {
    	    if (atomId1[i] != ' ')
    	        s1 += atomId1[i];
        }

    short i = 0;
    bool found = false;
        while (i != ATOMMAX && !found)
        {
    	    found = EqualsAtomName(s1,i);
    	    i++;
        }
    
    if (found)
    	return i;
    	
    return ATOMMAX;
    }    

public:

	    	/**#: [Description="Returns true is 1st string comes before
	    	       the 2nd in pdb order. Returns false otherwise"] */
    bool    operator()(const string& atomId1, const string& atomId2) const
    {

    // Find index number of first atom name.
    short idx1 = FindIndx(atomId1);
    
    // Find index number of second atom name.
    short idx2 = FindIndx(atomId2);
    	
    if (idx1 != ATOMMAX && idx2 != ATOMMAX)
	return idx1 < idx2;
    else if (idx1 != ATOMMAX && idx2 == ATOMMAX)
	return true;
    else if (idx1 == ATOMMAX && idx2 != ATOMMAX)
	return false;
    else 
	return atomId1 < atomId2;
    }


	    	/**#: [Description="The c style string version of the above
	    	       function"] */
    bool    operator()(const char* atomId1, const char* atomId2) const
    { 
    string s1 = atomId1;
    string s2 = atomId2;
    return operator()(s1,s2); 
    }    	    	
};

_BTL_END_NAMESPACE
  	    	        	    	        	    	 
#endif  // pdb_sort.h
