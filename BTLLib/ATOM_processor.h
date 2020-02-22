//
// ATOM_processor.h
//
// This file contains declarations and definitions for the ATOMProcessor class.
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1998 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2005 University College, Gower Street, London, U.K.
//
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Library General Public License as published 
// by the Free Software Foundation; either version 2 of the License, or 
// (at your option) any later version. This library is distributed in the 
// hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined(BTL_ATOMPROCESSOR_H)
#define BTL_ATOMPROCESSOR_H 1

#include "BTL.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

_BTL_BEGIN_NAMESPACE

using namespace std;


/**#: [Description ="This is simple class for reading ATOM records from 
	clean PDB format files."]
    [Summary = "reads PDB format files and extracts ATOM record data into vectors."]
    [Friends = "an output operator"]
    [Authors = "M.A.Williams"]
    [Files = "<A HREF=./btl/ATOMProcessor.h>ATOMProcessor.h</A>"]
    [Dependencies="None"]
*/

class ATOM_processor
{
public:

    typedef vector<int>        IdStore;
    typedef vector<BTL_REAL>   RealStore;
    typedef vector<char>       NameStore;

private:

    IdStore	    atom_no;
    NameStore       atom_name;
    NameStore       residue_name;
    NameStore       chain_name;
    IdStore	    residue_no;
    NameStore       alternate_name;
    RealStore	    coords;
    RealStore	    occupancy;
    RealStore	    b_factor;

	    	/**#: [Hidden] */
    void
    OpenFile(const char* fileName, ifstream& fileStream);

    void
    CloseFile(const char* fileName, ifstream& fileStream);

    NameStore&
    Alternate() { return alternate_name; }


public:

	    	/**#: [Description="Default constructor (does nothing)"] */

    ATOM_processor() {}

	    	/**#: [Description="Default destructor (does nothing)"] */

    ~ATOM_processor() {}

	    	/**#: [Description="Returns true if no coordinates have been read"] */
    bool
    empty() const { return coords.empty(); }

	    	/**#: [Description="Read a PDB format file with the given name,
                selecting only those chains whose id letter occurs in a given string
                and using the alternate set of coordinates identified by the 
                user-supplied character. The default reads records with any 
                chain id and either no alternate id or the first encounterd
                alternate id."] */
    void
    ReadFile(const char *fileName, const char *chains = " ", const char alternate = ' ');

	    	/**#: [Description="Return a reference to the atom number in a
                STL vector of int."] */
    IdStore&
    AtomNo() { return atom_no; }

	    	/**#: [Description="Return a reference to the atom n in a
                STL vector of char. Each name is stored as 4 chars."] */
    NameStore&
    AtomName() { return atom_name; }

	    	/**#: [Description="Return a reference to the atom names in a
                STL vector of char. Each name is stored as 3 chars."] */
    NameStore&
    ResidueName() { return residue_name; }

	    	/**#: [Description="Return a reference to the residue number in a
                STL vector of int."] */
    IdStore&
    ResidueNo() { return residue_no; }

	    	/**#: [Description="Return a reference to the chain names in a
                STL vector of char. Each name is stored as a single char."] */
    NameStore&
    ChainName() { return chain_name; }

	    	/**#: [Description="Return a reference to the read coordinates in a
                STL vector of float. Coordinates are stored in the following
                order: x,y,z,x,y,z,x,..."] */
    RealStore&
    Coords() { return coords; }

	    	/**#: [Description="Return a reference to the occupancy in a
                STL vector of float. "] */
    RealStore&
    Occupancy() { return occupancy; }

	    	/**#: [Description="Return a reference to the b_factor in a
                STL vector of float. "] */
    RealStore&
    BFactor() { return b_factor; }

    friend ostream&
    operator<<(ostream& os, const ATOM_processor& x);
};

ostream& operator<<(ostream& os, const ATOM_processor& x)
{

    // Pass information to the output stream in PDB format
 
    cout << showpoint << fixed << right;

    unsigned int nOfAtoms = x.atom_no.size();
    for (unsigned int i=0; i<nOfAtoms; i++)
    {
    	os << "ATOM  " 
           << setw(5) << x.atom_no[i] 
           << " "
    	   << x.atom_name[4*i]
    	   << x.atom_name[4*i+1]
    	   << x.atom_name[4*i+2]
    	   << x.atom_name[4*i+3]
           << " "
    	   << x.residue_name[3*i]
    	   << x.residue_name[3*i+1]
    	   << x.residue_name[3*i+2]
           << " "
           << x.chain_name[i]
           << setw(4) << x.residue_no[i] 
           << x.alternate_name[i]
           << "   "
 	   << setw(8) << setprecision(3) << x.coords[3*i] 
    	   << setw(8) << setprecision(3) << x.coords[3*i+1]
    	   << setw(8) << setprecision(3) << x.coords[3*i+2]
    	   << setw(6) << setprecision(2) << x.occupancy[i]
    	   << setw(6) << setprecision(2) << x.b_factor[i]
    	   << '\n';
    }
    os << '\n';
    return os;
}

void ATOM_processor::OpenFile(const char* fileName, ifstream& fileStream)
{
    fileStream.open(fileName);

    // check that file can be opened
    if (!fileStream.is_open())
    {
        cerr << "\n\nFile " << fileName << " cannot be opened" << endl;
        fileStream.clear();
        exit(1);
    }
}

void ATOM_processor::CloseFile(const char* fileName, ifstream& fileStream)
{
    fileStream.clear();
    fileStream.close();
}

void ATOM_processor::ReadFile(const char *fileName, const char *chains, const char alternate)
{ 
    // Open a file as an input stream object

    ifstream pdbin;
    OpenFile(fileName,pdbin); 

    string line;

    // To ensure efficient loading of the vectors we reserve sufficient
    //  contiguous memory for each vector in advance.
    // First read file to determine total number of ATOM records
    // N.B. It remains to be tested whether the additional read offsets 
    //      any memory allocation advantage

    unsigned int nOfAtoms = 0;
    while(getline(pdbin,line))
    {
        if (line.substr(0,4) == "ATOM") nOfAtoms++;
    }

    // Close the input file

    CloseFile(fileName,pdbin); 

    // Erase current contents of each vector if neccessary
    if (!empty())
    {
        atom_no.erase(atom_no.begin(),atom_no.end());
        atom_name.erase(atom_name.begin(),atom_name.end());
        residue_name.erase(residue_name.begin(),residue_name.end());
        chain_name.erase(chain_name.begin(),chain_name.end());
        residue_no.erase(residue_no.begin(),residue_no.end());
        alternate_name.erase(alternate_name.begin(),alternate_name.end());
        coords.erase(coords.begin(),coords.end());
        occupancy.erase(occupancy.begin(),occupancy.end());
        b_factor.erase(b_factor.begin(),b_factor.end());
    }

    // Reserve space for each new vector
    atom_no.reserve(nOfAtoms);
    atom_name.reserve(nOfAtoms*4);
    residue_name.reserve(nOfAtoms*3);
    chain_name.reserve(nOfAtoms);
    residue_no.reserve(nOfAtoms);
    alternate_name.reserve(nOfAtoms);
    coords.reserve(nOfAtoms*3);
    occupancy.reserve(nOfAtoms);
    b_factor.reserve(nOfAtoms);

    // Read in the data
    // PDB format is 6x,i5,x,a4,x,a3,x,a,i4,a,3x,3f8.3,2f6.2 fortran style

    // Open the file again at the beginning 
    OpenFile(fileName,pdbin); 

    // Read data into the appropriate vectors
    char first_alternate = alternate;
    string s; s = chains;     // separate initialization and assignment for string 

    while(getline(pdbin,line))
    {
        if (line.substr(0,4) == "ATOM")
	{
        //  if default behaviour or chain id matches a character in the user defined string
            if ((chains == " ") || (s.find(line[21]) < s.size()) )
	    {             
            //  if default behaviour or alternate id matches user defined character
                if(first_alternate == ' ') first_alternate = line[26];
                if(line[26] == ' '  || first_alternate == line[26])
                {
    	            atom_name.push_back(line[12]);
    	            atom_name.push_back(line[13]);
    	            atom_name.push_back(line[14]);
    	            atom_name.push_back(line[15]);
                    residue_name.push_back(line[17]);
                    residue_name.push_back(line[18]);
                    residue_name.push_back(line[19]);
                    chain_name.push_back(line[21]);
                    alternate_name.push_back(line[26]);

    	            atom_no.push_back(atoi(line.substr(6, 5).c_str()));
    	            residue_no.push_back(atoi(line.substr(22, 4).c_str()));

    	            for (int j=0; j<3; j++)
    	            {
    	                coords.push_back((BTL_REAL) atof(line.substr(30 + j*8,8).c_str()));
    	            }    
    	            occupancy.push_back((BTL_REAL) atof(line.substr(54, 6).c_str()));
    	            b_factor.push_back((BTL_REAL) atof(line.substr(60, 6).c_str()));
                }
	    }
        }
    }

    CloseFile(fileName,pdbin);

}

_BTL_END_NAMESPACE

#endif
