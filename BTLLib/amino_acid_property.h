//
// amino_acid_property.h
//
// This file contains the BTL_AminoAcidProperty functor.
//
// These classes are part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2004,2005 University College London, Gower Street, London, U.K.
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

#if !defined (BTL_AMINOACIDPROPERTY_H)
#define BTL_AMINOACIDPROPERTY_H 1

#include <map>
#include <cctype>
#include <iostream>
_BTL_BEGIN_NAMESPACE

using namespace std;

enum BTL_Scale { KyteHydropathy, MillerAccessibility, DarbyVolume,
                 RadzickaHydrophilicity, RadzickaHydrophobicity };

/**#: [Description ="This class defines a function object that will, given an
    amino acid single letter code with return a property value for that residue
    type."]
    [Summary = "functor that returns amino acid property values"] 
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/amino_acid_property.h>amino_acid_property.h</A>"]
    [Dependencies="none"]

*/

class amino_acid_property
{
private:
    
    typedef map<char,float,less<char> > value_type;
    value_type data;

	    	/**#: [Hidden] */

// Amino acid volume in cubic Angstroms as given in "Protein Structure" 
// by Darby, N.J. and Creighton, T.E. (IRL Press, 1993)

    static void Use_DC_Volume(value_type& values){     
    values['A'] =   67; 
    values['C'] =   86; 
    values['D'] =   91; 
    values['E'] =  109; 
    values['F'] =  135; 
    values['G'] =   48; 
    values['H'] =  118; 
    values['I'] =  124; 
    values['K'] =  135; 
    values['L'] =  124; 
    values['M'] =  124; 
    values['N'] =   96; 
    values['P'] =   90; 
    values['Q'] =  114; 
    values['R'] =  148; 
    values['S'] =   73; 
    values['T'] =   93; 
    values['V'] =  105; 
    values['W'] =  163; 
    values['Y'] =  141; 
    values['X'] =    0; }

	    	/**#: [Hidden] */
// Accessible surface areas in Angstroms squared by Miller, S., Janin, J., 
// Lesk, A.M., and Chothia, C. (1987) J.Mol.Biol., 196, 641
    static void Use_MJLC_AccessibleSurfaceArea(value_type& values){     
    values['A'] =  113; 
    values['C'] =  140; 
    values['D'] =  151; 
    values['E'] =  183; 
    values['F'] =  218; 
    values['G'] =   85; 
    values['H'] =  194; 
    values['I'] =  182; 
    values['K'] =  211; 
    values['L'] =  180; 
    values['M'] =  204; 
    values['N'] =  158; 
    values['P'] =  143; 
    values['Q'] =  189; 
    values['R'] =  241; 
    values['S'] =  122; 
    values['T'] =  146; 
    values['V'] =  160; 
    values['W'] =  259; 
    values['Y'] =  229; 
    values['X'] =    0;}

	    	/**#: [Hidden] */
// Hydrophilicity of side chains in kJ/mol by Radzicka, A., and Wolfenden, R.
// (1988) Biochemistry, 27, 1664 
    static void Use_RW_Hydrophilicity(value_type& values){     
    values['A'] =  -1.89; 
    values['C'] = -15.24; 
    values['D'] = -56.02; 
    values['E'] = -53.05; 
    values['F'] = -13.23; 
    values['G'] =    0.0; 
    values['H'] = -53.17; 
    values['I'] =   -1.0; 
    values['K'] = -50.02; 
    values['L'] =  -0.46; 
    values['M'] = -16.25; 
    values['N'] = -50.69; 
    values['P'] =    0.0; 
    values['Q'] = -49.43; 
    values['R'] = -93.66; 
    values['S'] = -31.29; 
    values['T'] = -30.53; 
    values['V'] =  -1.68; 
    values['W'] = -34.73; 
    values['Y'] =  -35.7; 
    values['X'] =    0.0; }

	    	/**#: [Hidden] */
// Hydropathy index by Kyte, J., and Doolittle, R.F., (1982), J.Mol.Biol., 
// 157, 105-132
    static void Use_KD_Hydropathy(value_type& values){     
    values['A'] =  1.8;
    values['C'] =  2.5;
    values['D'] =  -3.5;
    values['E'] =  -3.5;
    values['F'] =  2.8;
    values['G'] =  -0.4;
    values['H'] =  -3.2;
    values['I'] =  4.5;
    values['K'] =  -3.9;
    values['L'] =  3.8;
    values['M'] =  1.9;
    values['N'] =  -3.5;
    values['P'] =  -1.6;
    values['Q'] =  -3.5;
    values['R'] =  -4.5;
    values['S'] =  -0.8;
    values['T'] =  -0.7;
    values['V'] =  4.2;
    values['W'] =  -0.9;
    values['Y'] =  -1.3;
    values['X'] =  0.0;}

	    	/**#: [Hidden] */
// Hydrophobicity of side chain analogues in kJ/mol by Radzicka, A., 
// and Wolfenden, R. (1988) Biochemistry, 27, 1664 
    static void Use_RW_Hydrophobicity(value_type& values){     
    values['A'] =  -3.65;
    values['C'] =  -1.42;
    values['D'] =  40.57;
    values['E'] =  32.55;
    values['F'] =  -8.57;
    values['G'] =    0.0;
    values['H'] =  23.52;
    values['I'] = -16.71;
    values['K'] =  27.25;
    values['L'] = -16.71;
    values['M'] =  -5.92;
    values['N'] =  21.83;
    values['P'] =    0.0;
    values['Q'] =  27.21;
    values['R'] =  66.61;
    values['S'] =  18.23;
    values['T'] =  14.74;
    values['V'] = -13.02;
    values['W'] =  -5.84;
    values['Y'] =   4.54;
    values['X'] =    0.0;}

        
public:

	    	/**#: [Description="Default constructor. The default property
		type is KyteHydropathy."] */
    amino_acid_property()
    {
        amino_acid_property(KyteHydropathy);
    }
    
	    	/**#: [Description="Construct function object which will return
		property values of a given type."] */
    amino_acid_property(const BTL_Scale whatScale)
    { 
    switch(whatScale)
        {
        case KyteHydropathy:
    
    	    Use_KD_Hydropathy(data);
    	    break;
    	
        case MillerAccessibility:
    
    	    Use_MJLC_AccessibleSurfaceArea(data);
    	    break;
    	
        case DarbyVolume:
    
    	    Use_DC_Volume(data);
    	    break;
    	
        case RadzickaHydrophilicity:
    
    	    Use_RW_Hydrophilicity(data);
    	    break;
    	
        case RadzickaHydrophobicity:
    
    	    Use_RW_Hydrophilicity(data);
    	    break;
        }
    }    

	    	/**#: [Description="Operator that will return a property value
		given a single letter amino acid code."] */
    float 
    operator()(const char aacode) const
    {
       char code = toupper(aacode);
       value_type::const_iterator i = data.find(code);
           if (i == data.end())
           {
    	        cerr << "!!!Property not available for :" << aacode << endl;
    	        return 0.0;
           }
       return (*i).second;
    }

};

_BTL_END_NAMESPACE

#endif
