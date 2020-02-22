//
// PIR_processor.h
//
// This file contains declarations for the PIRProcessor class.
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997, 1998 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2004, 2005 University College, Gower Street, London, U.K.
//
// This library is free software; you can redistribute it and/or modify 
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; either version 2 of the License, or 
// (at your option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined(BTL_PIRPROCESSOR_H)
#define BTL_PIRPROCESSOR_H 1

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

_BTL_BEGIN_NAMESPACE

using namespace std;


/**#: [Description ="This class is a simple PIR sequence file processor. It is
    not complete and is only included for use by the demo programs that come
    with this library. It can parse a PIR format containing any number
    of sequences and provide the data as a vector<string> or as individual
    strings."]
    [Summary = "simple class for parsing PIR format sequence files"] 
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/PIRProcessor.h>PIRProcessor.h</A>,
              <A HREF=./btl/PIRProcessor.cpp>PIRProcessor.cpp</A>"]
    [Dependencies="none"]

*/

class PIR_processor
{
private:

    vector<string> seqs;
    
	    	/**#: [Hidden] */
    void 
    OpenFile(const char* fileName, ifstream& fileStream);
    
	    	/**#: [Hidden] */
    int 
    HowManySeqs(const char* fileName, ifstream& fileStream);

public:

	    	/**#: [Description="Default constructor (does nothing)"] */
    PIR_processor() {}
    
	    	/**#: [Description="Return true if no sequences have been
	    	       read in, false otherwise."] */
    bool
    empty() const { return seqs.empty(); } 
    
	    	/**#: [Description="Read a file with a given name. If this 
	    	       function is called more than once, the new sequences
	    	       are added to those already there."] */
    void
    ReadFile(const char *fileName);
    
	    	/**#: [Description="Return any sequences that been read in."] */
    vector<string>&
    Seqs() { return seqs; }
    
	    	/**#: [Description="Return a single sequence, given a valid
	    	       index number"] */
    string&
    Seq(const unsigned int iseq) 
    {   
    	if (iseq>=seqs.size()) 
    	{
    	    cerr << "!!Invalid sequence no." << endl;
    	    exit(1);
    	}
    	return seqs[iseq];
    }
};

void PIR_processor::OpenFile(const char* fileName, ifstream& fileStream)
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

int PIR_processor::HowManySeqs(const char* fileName, ifstream& fileStream)
{
    int nseqs=0;
    fileStream.open(fileName);
    string line;
    while(getline(fileStream,line))
     	if( line[0] == '>' ) nseqs++;
    fileStream.clear();
    fileStream.close();

    cout << nseqs << " sequences found in " << fileName << "\n\n";
    return nseqs;
}    

void PIR_processor::ReadFile(const char *fileName)
{ 
    ifstream fin;

    int nseqs = HowManySeqs(fileName,fin);
    
    OpenFile(fileName,fin);
    //
    // Extract sequences from file 
    //
    int iseq = 0;
    string line;

    // Read each sequence in turn
    //
    while(iseq<nseqs)
    {
        char c;

        // Read until the beginning of a sequence header
        //
        do
        {
            if (!getline(fin,line)) break;
            c = line[0];
        }
        while(c != '>');
            
        if (c!='>') cerr << "\n!!!Error reading " << fileName 
                         << ": initial > expected.\n" << endl;

    	// Skip comment field
    	//
        if (!getline(fin,line)) break;
    
    	// Read until the end of the sequence
    	//
    	string seq;

    	do
        {
            if (!getline(fin,line)) break;
            for (int i=0; i<line.size(); i++)
            {
                c=line[i];
                if(c == '\n' || c == '\0' || c == '*')
                    break;
                if (c != ' ')
                    seq += c;
            }
    	}
        while(c!='*');
        
        if (c!='*') cerr << "\n!!!Error reading " << fileName 
                         << ": terminal * expected.\n" << endl;

    	if (!seq.empty()) seqs.push_back(seq);

        iseq++;    
    }
    fin.clear();
    fin.close();
}

_BTL_END_NAMESPACE

#endif
