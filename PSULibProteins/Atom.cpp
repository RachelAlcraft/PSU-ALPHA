#include "Atom.h"
#include <iostream>
#include <string>
#include <sstream>
#include "LogFile.h"
#include <StringManip.h>

Atom::Atom(string pdb_code, string atom_string)
{		
	pdbCode = pdb_code;
	// 7 - 11        Integer       serial       Atom  serial number.
	atomId = atol(StringManip::trim(atom_string.substr(6, 5)).c_str());
	//13 - 16        Atom          name         Atom name.
	elementName = StringManip::trim(atom_string.substr(12, 4));
	//18 - 20        Residue name  resName      Residue name.
	aminoCode = StringManip::trim(atom_string.substr(17, 3));
	//22             Character     chainID      Chain identifier.
	chainId = StringManip::trim(atom_string.substr(21, 1));
	//23 - 26        Integer       resSeq       Residue sequence number.
	aminoId = atol(StringManip::trim(atom_string.substr(22, 5)).c_str());
	//31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
	//39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
	//47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
	string x_c = StringManip::trim(atom_string.substr(30, 8));	
	string y_c = StringManip::trim(atom_string.substr(38, 8));
	string z_c = StringManip::trim(atom_string.substr(46, 8));
	coords = Coordinates(atof(x_c.c_str()), atof(y_c.c_str()), atof(z_c.c_str()));
	shifted_coords = Coordinates(atof(x_c.c_str()), atof(y_c.c_str()), atof(z_c.c_str()));
	//77 - 78        LString(2)    element      Element symbol, right-justified.
	elementType = StringManip::trim(atom_string.substr(76, 2));

	stringstream di;
	di << pdbCode << chainId << aminoId;
	dataId = di.str();
}
Atom::~Atom()
{
	_bonds.clear();
}

void Atom::applyShift(double x_shift, double y_shift, double z_shift, bool applyToOriginal)
{
	if (applyToOriginal)
	{
		shifted_coords.x = coords.x + x_shift;
		shifted_coords.y = coords.y + y_shift;
		shifted_coords.z = coords.z + z_shift;
	}
	else
	{
		shifted_coords.x = shifted_coords.x + x_shift;
		shifted_coords.y = shifted_coords.y + y_shift;
		shifted_coords.z = shifted_coords.z + z_shift;
	}

}

/*string Atom::trim(string string_to_trim)
{
	string string_trimmed = string_to_trim;
	size_t startpos = string_trimmed.find_first_not_of(" ");
	size_t endpos = string_trimmed.find_last_not_of(" ");
	if (startpos == string::npos)
		string_trimmed = "";
	else if (endpos == string::npos)
		string_trimmed = string_trimmed.substr(startpos);
	else
		string_trimmed = string_trimmed.substr(startpos, endpos - startpos + 1);
	return string_trimmed;
}*/

void Atom::printAtom()
{
	LogFile::getInstance()->writeMessage(getDescription());
}

string Atom::getDescription()
{
	stringstream ss;
	ss << atomId << "\t";
	ss << elementType << "=" << elementName << "\t";
	ss << aminoCode << "=" << aminoId << "\t";
	ss << "Chain=" << chainId << "\t";
	ss << "(x,y,z)=(" << coords.x << "," << coords.y << "," << coords.z << ")";
	return ss.str();
}

GeoVector Atom::vectorDifference(Atom* comp)
{
	return GeoVector(comp->coords.x - coords.x, comp->coords.y - coords.y, comp->coords.z - coords.z);
}

double Atom::atomicDistance(Atom* comp)
{
	GeoVector v = GeoVector(comp->coords.x - coords.x, comp->coords.y - coords.y, comp->coords.z - coords.z);
	return v.getMagnitude();
}