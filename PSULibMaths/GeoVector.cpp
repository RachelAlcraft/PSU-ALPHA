#include "GeoVector.h"
#include <math.h>

GeoVector::GeoVector()
{
	x = 0;
	y = 0;
	z = 0;
}
GeoVector::GeoVector(double X, double Y, double Z)
{
	x = X;
	y = Y;
	z = Z;
}

GeoVector::GeoVector(GeoCoords a, GeoCoords b)
{
	x = b.x - a.x;
	y = b.y - a.y;
	z = b.z - a.z;
}

double GeoVector::angle(GeoVector b)
{
	//angle betwen 2 vectors in 3 dimensions on the plain formed (between outer pairs of dihedral angles)
	GeoVector a = *this;
	//dot product
	double dot = (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
	//magnitude
	double ma2 = a.magnitude();
	double mb2 = b.magnitude();
	double div = pow(ma2, 0.5) * pow(mb2, 0.5);
	//cos
	double cos_theta = dot / div;
	//inverse cos
	double theta = acos(cos_theta);
	theta = (theta / PI) * 180;//convert to degrees
	//if (theta > 180)
	//	theta = theta - 180;
	theta = round(theta);
	//ATAN2 version?	I prefer this but the decision on sign is effectively the same thing.
	return theta;
}

double GeoVector::magnitude()
{
	double mag = (x * x) + (y * y) + (z * z);
	return sqrt(mag);
}

GeoVector GeoVector::getCrossProduct(GeoVector B)
{
	GeoVector A = *this;
	double px = (A.y * B.z) - (A.z * B.y);
	double py = (A.z * B.x) - (A.x * B.z);
	double pz = (A.x * B.y) - (A.y * B.x);
	return GeoVector(px, py, pz);
}

double GeoVector::getDotProduct(GeoVector B)
{
	GeoVector A = *this;
	double dot = (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
	return dot;
}

double GeoVector::getOrthogonalDistance(GeoCoords p)
{
	return 0.0;
}

/*double GeoVector::getOrthogonalDistance(GeoCoord a, GeoCoord b, GeoCoord c)
{//http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
	GeoVector ca = GeoVector(c,a);
	GeoVector cb = GeoVector(c,b);
	GeoVector ba = GeoVector(b,a);
	GeoVector caDOTcb = ca.getCrossProduct(cb);
	double mag = ba.getMagnitude();
	double distance = caDOTcb.getMagnitude();
	distance = distance / mag;
	return distance;
}*/

double GeoVector::getMagnitude()
{
	GeoVector A = *this;
	double mag = (A.x * A.x) + (A.y * A.y) + (A.z * A.z);
	return sqrt(mag);
}
