//
// XYZ_processor.h
//
// This file contains declarations for the XYZProcessor class.
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London WC1E 7HX, U.K.
// (classlib@mail.cryst.bbk.ac.uk)
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

#if !defined(BTL_XYZPROCESSOR_H)
#define BTL_XYZPROCESSOR_H 1

#include "BTL.h"

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description ="This is simple class for reading xyz format files.
     It is really only included in the library so that real data can be used
     in the test programs."]
    [Summary = "reads xyz format coordinate files and puts the data in vectors"]
    [Friends = "an output operator"]
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/XYZProcessor.h>XYZProcessor.h</A>"]
    [Dependencies="None"]
*/

class XYZ_processor
{
    typedef vector<float>   CoordStore;
    typedef vector<char>    ElementStore;

private:

    CoordStore	    coords;
    ElementStore    elements;
    string  	    title;

	    	/**#: [Hidden] */
    void
    OpenFile(const char* fileName, ifstream& fileStream);


public:

	    	/**#: [Description="Default constructor (does nothing)"] */

    XYZ_processor() {}

	    	/**#: [Description="Returns true if no coordinates have been read"] */
    bool
    empty() const { return coords.empty(); }

	    	/**#: [Description="Read an xyz format file with the given name"] */
    void
    ReadFile(const char *fileName);

	    	/**#: [Description="Return a reference to the read coordinates in a
                STL vector of float. Coordinates are stored in the following
                order: x,y,z,x,y,z,x,..."] */
    CoordStore&
    Coords() { return coords; }

	    	/**#: [Description="Return a reference to the read atomic elements in a
                STL vector of char. Each element is stored as 2 chars."] */
    ElementStore&
    Elements() { return elements; }

    friend ostream&
    operator<<(ostream& os, const XYZ_processor& x);
};

ostream& operator<<(ostream& os, const XYZ_processor& x)
{
    os << x.title << "\n\n";

    unsigned int nAtoms = x.coords.size()/3;
    for (unsigned int i=0, j=0, k=0; i<nAtoms; i++)
    {
    	os << setw(3) << i << " : " 
    	   << x.elements[j++] 
    	   << x.elements[j++]    	   
    	   << " " << setw(10) << x.coords[k++] 
    	   << " " << setw(10) << x.coords[k++]
    	   << " " << setw(10) << x.coords[k++]
    	   << '\n';
    }
    os << '\n';
    return os;
}

void XYZ_processor::OpenFile(const char* fileName, ifstream& fileStream)
{
    fileStream.open(fileName);

    // check that file can be opened
    if (!fileStream.good())
    {
        cerr << "\n\nFile " << fileName << " cannot be opened" << endl;
        exit(1);
    }
    
}
void XYZ_processor::ReadFile(const char *fileName)
{ 
    ifstream fin;
    OpenFile(fileName,fin);

    string line;

    // Number of coords
    //
    if (!getline(fin,line)) return;
    unsigned int nAtoms = atoi(line.c_str());
    
    elements.erase(elements.begin(),elements.end());
    elements.reserve(nAtoms*2);
    coords.erase(coords.begin(),coords.end());
    coords.reserve(nAtoms*3);
    
    // Title
    //
    if (!getline(fin,line)) return;
    title = line;
    
    // Read coords
    //
    for (unsigned int i=0; i<nAtoms; i++)
    {    	
    	getline(fin,line,' ');
    	elements.push_back(line[0]);
    	elements.push_back(line[1]);
        getline(fin,line,' ');
    	coords.push_back(atof(line.c_str()));
        getline(fin,line,' ');
    	coords.push_back(atof(line.c_str()));
        getline(fin,line);
    	coords.push_back(atof(line.c_str()));
    }
} 

_BTL_END_NAMESPACE

#endif
