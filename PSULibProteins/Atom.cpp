#include "Atom.h"
#include <iostream>
#include <string>
#include <sstream>
#include "LogFile.h"
#include <StringManip.h>
#include <Torsion.h>

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
	coords = GeoCoords(atof(x_c.c_str()), atof(y_c.c_str()), atof(z_c.c_str()));
	shifted_coords = GeoCoords(atof(x_c.c_str()), atof(y_c.c_str()), atof(z_c.c_str()));
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
	return GeoVector(coords, comp->coords);
}

double Atom::atomicDistance(Atom* comp, bool shifted)
{
	if (shifted)
	{
		GeoVector v = GeoVector(shifted_coords, comp->shifted_coords);
		return v.getMagnitude();
	}
	else
	{
		GeoVector v = GeoVector(coords, comp->coords);
		return v.getMagnitude();
	}
}

void Atom::applyTransformation(GeoTransformations* trans)
{
	GeoCoords newCoords = trans->applyTransformation(coords);
	shifted_coords = newCoords;
}

// BONDS ################################################################
AtomBond::AtomBond(string ss, Atom* a1, Atom* a2)
{
	_A1 = a1;
	_A2 = a2;	
	_SS = ss;
	_atomString = a1->elementName + "-" + a2->elementName;

}
double AtomBond::getValue()
{
	GeoVector a = _A1->vectorDifference(_A2);
	return a.getMagnitude();
}
// ANGLES ################################################################
AtomAngle::AtomAngle(string ss, Atom* a1, Atom* a2, Atom* a3):AtomBond(ss,a1,a2)
{
	_A1 = a1;
	_A2 = a2;
	_A3 = a3;
	_atomString += "-" + a3->elementName;

}

double AtomAngle::getValue()
{
	GeoVector a = _A2->vectorDifference(_A1);
	GeoVector b = _A2->vectorDifference(_A3);
	double angle = a.angle(b);
	return angle;
}

// TORSIONS ################################################################
AtomTorsion::AtomTorsion(string ss, Atom* a1, Atom* a2, Atom* a3, Atom* a4) :AtomBond( ss, a1, a2)
{
	_A1 = a1;
	_A2 = a2;
	_A3 = a3;
	_A4 = a4;
	_atomString += "-" + a3->elementName + "-" + a4->elementName;
}

double AtomTorsion::getValue()
{
	vector < Atom *> atoms; 
	atoms.push_back(_A1);
	atoms.push_back(_A2);
	atoms.push_back(_A3);
	atoms.push_back(_A4);
	return Torsion::getDihedralAngle(atoms);
}
