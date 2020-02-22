//
// PPM_processor.h
//
// This file contains the PPM_processor class.
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

#if !defined (BTL_PPMPROCESSOR_H)
#define BTL_PPMPROCESSOR_H 1

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>

_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description ="This class is a simple ascii PPM image file processor. It 
    is not complete and is only included for use by the demo programs that come
    with this library."] 
    [Summary = "simple class for parsing ascii PPM format image files"] 
    [Authors = "W.R.Pitt"]
    [Friends = "ostream operator"]
    [Files = "<A HREF=./btl/PPMProcessor.h>PPMProcessor.h</A>,
              <A HREF=./btl/PPMProcessor.cpp>PPMProcessor.cpp</A>"]
    [Dependencies="none"]

*/


class PPM_processor
{
private:

    unsigned int width, height;
    vector<double> imageData;
    
	    	/**#: [Hidden] */
    void 
    OpenFile(const char* fileName, ifstream& fileStream);

	    	/**#: [Hidden] */
    void 
    GetDimensions(ifstream& fileStream);

	    	/**#: [Hidden] */
    void 
    GetImageData(unsigned int size, ifstream& fileStream);
    
public:
    
	    	/**#: [Description="Read in a file with a given name."] */
    void
    ReadFile(const char *fileName);
    
	    	/**#: [Description="Get image data."] */
    vector<double>& 
    ReadImageData() { return imageData; }

	    	/**#: [Description="Write image data."] */
    void
    WriteImageData(const vector<double>& idata, 
                   const unsigned int iwidth,
                   const unsigned int iheight); 

	    	/**#: [Description="Read width of image."] */
    unsigned int
    ReadWidth() const { return width; }
     
	    	/**#: [Description="Read image height."] */
    unsigned int
    ReadHeight() const { return height; }
     
    friend ostream&
    operator<<(ostream &os, const PPM_processor &v);
    
};

void PPM_processor::OpenFile(const char* fileName, ifstream& fileStream)
{
    fileStream.open(fileName);

    // check that file can be opened
    if (!fileStream.good())
    {
        cerr << "\n\nFile " << fileName << " cannot be opened" << endl;
        exit(1);
    }
    
}

void PPM_processor::GetDimensions(ifstream& fileStream)
{
    string line;
    
    // Read and ignore 1st line and any comments
    if (!getline(fileStream,line))
    {
    	cerr << "!! Error reading image dimensions." << endl;
    	exit(1);
    }
    
    do getline(fileStream.line); while (line[0] == '#');
    
    getline(fileStream,line,' ');
    width = atoi(line.c_str());
    getline(fileStream,line);
    height = atoi(word);
}

void PPM_processor::GetImageData(unsigned int size, ifstream& fileStream)
{
    unsigned int i=0;	
    char line[80];
    while (i<size)
    {	
    	if (!fileStream.getline(line,80))
    	{
    	    cerr << "!! Error reading image data." << endl;
    	    exit(1);
    	}
    	char intensity[4];
    	intensity[3] = 0;
    	unsigned int j=0;
    	while (i<size && j<15)
    	{
    	    for (int k=0; k<3; k++)
    	    {
    	    	strncpy(intensity, line + j++ * 4, 3);
    	    	imageData.push_back((double) atoi(intensity));
    	    }
    	    i++;
    	}
    }
}

void PPM_processor::WriteImageData(const vector<double>& idata, 
                   		  const unsigned int iwidth,
                    		  const unsigned int iheight)
{
    // Check dimensions make sense
    if (idata.size() != iwidth * iheight * 3)
    {
    	cerr << "!! In consistency in input image data." << endl;
    	exit(1); 
    }
    imageData = idata;
    width = iwidth;
    height = iheight;
}

void PPM_processor::ReadFile(const char *fileName)
{
    ifstream fin;
    OpenFile(fileName,fin);
    
    GetDimensions(fin);
    
    cerr << "Image Dimensions : " << width << " x " << height << '\n';
    
    unsigned int size = width * height;
    imageData.reserve(size*3);
    GetImageData(size,fin);
    fin.close();                
}    
    
ostream& operator<<(ostream &os, const PPM_processor &p)
{
   	 
    vector<double>::const_iterator i=p.imageData.begin();

    double maxval, minval;
    if (i != p.imageData.end())
    {
    	maxval=*i;
    	minval=*i;
    }
    	
    while (i!=p.imageData.end())
    {
    	double val = *i++;
    	if (val > maxval) maxval=val;    
    	if (val < minval) minval=val;    
    }
    cerr << "Max : " << maxval << " Min : " << minval << '\n';
    
    maxval -= minval;
    const long maxIntensity = 255;
    
    // Write header info
    os << "P3\n" 
       << "# CREATOR: The BTL's PPM_processor class DATE: " << __DATE__ 
       << '\n'
       << setw(3) << p.width << " " << p.height << '\n'
       << maxIntensity << '\n';
 
        
    i=p.imageData.begin();
    while (i!=p.imageData.end())
    {

    	int j=0;
    	while (j<15 && i!=p.imageData.end())
    	{
    	    double val = *i++;
   	    val = (val-minval)/maxval*maxIntensity;
    	    os << setw(3) << (long) val << " ";
    	    j++;
    	}
    	os << '\n';
    }
    return os;
}
    
_BTL_END_NAMESPACE

#endif
