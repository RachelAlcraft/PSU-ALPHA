#pragma once
#include <GeoCoords.h>
class GeoVector
{
	//get angle given another vector
	//apply vector to a point
	//some shifts around turns??
public://lazy public interface
	double x;
	double y;
	double z;
	const double PI = 3.141592653589793238463;
public:
	GeoVector();
	GeoVector(double, double, double);
	GeoVector(GeoCoords a, GeoCoords b);
	double angle(GeoVector);//the angle this vector makes with another vector (can start anywhere)
	double magnitude();
	GeoVector getCrossProduct(GeoVector);
	double getDotProduct(GeoVector);
	double getOrthogonalDistance(GeoCoords p);//shortest distance from p to the vector
	double getMagnitude();

};

