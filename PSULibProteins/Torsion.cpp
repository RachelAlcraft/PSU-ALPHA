#include "Torsion.h"
#include "GeoVector.h"
#include "GeoPlane.h"
/*
	Bonds are formed in the order
	C'-N-CA-C-N''-CA''
*/

Torsion::Torsion(string aa_name, int aa_id, Atom* N, Atom* CA, Atom* C)
{
	aa = aa_name;
	id = aa_id;
	_N = N;
	_CA = CA;
	_C = C;
}
BackboneTorsion::BackboneTorsion(string aa_name, int aa_id, Atom* Cp, Atom* N, Atom* CA, Atom* C, Atom* Npp, Atom* CApp) :Torsion(aa, id, N, CA, C)
{
	_Cp = Cp;
	_Npp = Npp;
	_CApp = CApp;
}
SidechainTorsion::SidechainTorsion(string aa_name, int aa_id, Atom* N, Atom* CA, Atom* C, Atom* AG1, Atom* AD1) :Torsion(aa, id, N, CA, C)
{
	_AG1 = AG1;
	_AD1 = AD1;
}

double Torsion::getDihedralAngle(Atom* atmA, Atom* atmB, Atom* atmC, Atom* atmD)
{
	//http://xiang-jun.blogspot.com/2009/10/how-to-calculate-torsion-angle.html
	GeoVector ab = atmA->vectorDifference(atmB);
	GeoVector cb = atmC->vectorDifference(atmB);
	GeoVector cd = atmC->vectorDifference(atmD);

	GeoVector p = ab.getCrossProduct(cb);
	GeoVector q = cb.getCrossProduct(cd);

	double dot = p.getDotProduct(q);
	double magP = p.getMagnitude();
	double magQ = q.getMagnitude();

	double theta = acos(dot / (magP * magQ));

	//Now check the sign
	GeoVector r = p.getCrossProduct(q);
	double dotsign = r.getDotProduct(cb);
	if (dotsign > 0)
		theta *= -1;
	double theta_deg = (theta / 3.141592653589793238463) * 180;//convert to degrees 		
	return round(theta_deg);
}

double BackboneTorsion::getPhi()
{//Cp:N:CA:C
	return getDihedralAngle(_Cp, _N, _CA, _C);
}
double BackboneTorsion::getPsi()
{//N:CA:C:Npp
	return getDihedralAngle(_N, _CA, _C, _Npp);
}
double BackboneTorsion::getOmega()
{//CA:C:Npp:CApp
	return getDihedralAngle(_CA, _C, _Npp, _CApp);
}
double SidechainTorsion::getChi1()
{//N:CA:C:AG1
	return 0;// getDihedralAngle(_N, _CA, _C, _AG1);
}
double SidechainTorsion::getChi2()
{//CA:C:AG1:AD1		
	return 0;// getDihedralAngle(_CA, _C, _AG1, _AD1);
}


