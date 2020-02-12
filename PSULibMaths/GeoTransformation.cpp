#include "GeoTransformation.h"
#include <math.h>


GeoTransformations::GeoTransformations()
{
	

}
GeoTransformations::~GeoTransformations()
{	
	for (unsigned int i = 0; i < transformations.size(); ++i)
	{
		if (transformations[i])
			delete transformations[i];
	}
	transformations.clear();
}
GeoTransformations::GeoTransformations(GeoCoords A, GeoCoords B, GeoCoords C)
{
	//This is the constructor to create a transformation that maps the 3 given points 
	//onto the origin and flat against the plane xz.
	MapA = A;
	MapB = B;
	MapC = C;
	TranslateRelativeToOrigin* t1 = new TranslateRelativeToOrigin(A);
	A = t1->applyTransformation(A);
	B = t1->applyTransformation(B);
	C = t1->applyTransformation(C);
	transformations.push_back(t1);

	RotateTo_Y_Is_Zero_AboutOrigin* t2 = new RotateTo_Y_Is_Zero_AboutOrigin(B);
	B = t2->applyTransformation(B);
	C = t2->applyTransformation(C);
	transformations.push_back(t2);

	RotateTo_Z_Is_Zero_AboutOrigin* t3 = new RotateTo_Z_Is_Zero_AboutOrigin(B);
	B = t3->applyTransformation(B);
	C = t3->applyTransformation(C);
	transformations.push_back(t3);

	RotateTo_Y_Is_Zero_OverX_Axis* t4 = new RotateTo_Y_Is_Zero_OverX_Axis(C);
	C = t4->applyTransformation(C);
	transformations.push_back(t4);	
}

GeoCoords GeoTransformations::applyTransformation(GeoCoords point)
{
	GeoCoords movedPoint = point;
	for (unsigned int i = 0; i < transformations.size(); ++i)
	{
		point = transformations[i]->applyTransformation(point);
	}
	return point;
}

//TRANSLATION///////////////////////////////////////////////////////////////////////////////////////////////////////////////
TranslateRelativeToOrigin::TranslateRelativeToOrigin(GeoCoords A):GeoTransform()
{
	//This finds the vector that translates a fixed boddy relative to tge origin for the reference point p
	V = GeoVector(A, GeoCoords(0, 0, 0));
}
GeoCoords TranslateRelativeToOrigin::applyTransformation(GeoCoords point)
{	
	return V.movePoint(point);
}
//ROTATION ABOUT THE ORIGIN////////////////////////////////////////////////////////////////////////////////////////////////////////////
RotateTo_Y_Is_Zero_AboutOrigin::RotateTo_Y_Is_Zero_AboutOrigin(GeoCoords A) :GeoTransform()
{
	//Z will remain unchanged so make a temporary no z vector
	GeoCoords pNoZ(A.x, A.y, 0);
	//We are mapping from A to 0, then 0 to B, so we have an iscoseles triangle |AO|==|OB| and AB
	GeoVector AO(pNoZ, GeoCoords(0, 0, 0));
	double magAO = AO.getMagnitude();
	GeoVector OB(magAO,0,0);//moving into +ve quadrant
	GeoVector AB = OB - AO;
	double magAB = AB.getMagnitude();
	//use cosine rule
	//Find theta with the cosine rule	
	double magAO2 = pow(magAO, 2);
	double magAB2 = pow(magAB, 2);
	double costheta = ((2 * magAO2) - magAB2) / (2 * magAO2);
	theta = acos(costheta);// in radians		
}
GeoCoords RotateTo_Y_Is_Zero_AboutOrigin::applyTransformation(GeoCoords point)
{
	//Triangle has 2 sides length a and a side length b
	GeoCoords pNoZ(point.x, point.y, 0);
	GeoVector AO(pNoZ, GeoCoords(0, 0, 0));
	double magAO = AO.getMagnitude();
	if (abs(magAO) > 0.0000001)
	{
		int Qfrom = 1;
		if (point.y < 0 && point.x > 0)
			Qfrom = 4;
		else if (point.y < 0 && point.x < 0)
			Qfrom = 3;
		else if (point.y > 0 && point.x < 0)
			Qfrom = 2;

		//If this is now the hypotanuse of a right-angled triangle with the x-axis
		double sinT = abs(point.y) / magAO;
		double T = asin(sinT);
		//Now we can subtract theta and we have the angle with the x-axis fort the lower side of the triangle (a diagram would help!)
		double t = T - theta;
		if (Qfrom == 2 || Qfrom == 4)
			t = T + theta;
		//The new x and y are the oints at the end of this new triangle
		//but what if we move into another quadrant?
		int Qto = Qfrom;
		if (t < 0 && Qfrom == 1)
		{
			t = theta - T;
			Qto = 4;
		}//this can't be the best way to do this?
		if (t < 0 && Qfrom == 3)
		{
			t = theta - T;
			Qto = 2;
		}

		double newY = sin(t) * magAO;		
		double newX = cos(t) * magAO;
		if (Qto == 3 || Qto == 4)
			newY *= 1;
		if (Qto == 2 || Qto == 3)
			newX *= 1;
				
		return GeoCoords(newX, newY, point.z);
	}
	else//TODO haven't sorted out the quadrants when x is negative
	{
		return point;
	}
}
RotateTo_Z_Is_Zero_AboutOrigin::RotateTo_Z_Is_Zero_AboutOrigin(GeoCoords A) :GeoTransform()
{
	//Z will remain unchanged so make a temporary no z vector
	GeoCoords pNoZ(A.x, 0, A.z);
	//We are mapping from A to 0, then 0 to B, so we have an iscoseles triangle |AO|==|OB| and AB
	GeoVector AO(pNoZ, GeoCoords(0, 0, 0));
	double magAO = AO.getMagnitude();
	GeoVector OB(magAO, 0, 0);//moving into +ve quadrant
	GeoVector AB = OB - AO;
	double magAB = AB.getMagnitude();
	//use cosine rule
	//Find theta with the cosine rule	
	double magAO2 = pow(magAO, 2);
	double magAB2 = pow(magAB, 2);
	double costheta = ((2 * magAO2) - magAB2) / (2 * magAO2);
	theta = acos(costheta);// in radians	
}
GeoCoords RotateTo_Z_Is_Zero_AboutOrigin::applyTransformation(GeoCoords point)
{
	//Triangle has 2 sides length a and a side length b
	GeoCoords pNoY(point.x, 0, point.z);
	GeoVector AO(pNoY, GeoCoords(0, 0, 0));
	double magAO = AO.getMagnitude();
	if (abs(magAO) > 0.0000001)
	{
		int Qfrom = 1;
		if (point.z < 0 && point.x > 0)
			Qfrom = 4;
		else if (point.z < 0 && point.x < 0)
			Qfrom = 3;
		else if (point.z > 0 && point.x < 0)
			Qfrom = 2;

		//If this is now the hypotanuse of a right-angled triangle with the x-axis
		double sinT = abs(point.z) / magAO;
		double T = asin(sinT);
		//Now we can subtract theta and we have the angle with the x-axis fort the lower side of the triangle (a diagram would help!)
		double t = T - theta;
		if (Qfrom == 2 || Qfrom == 4)
			t = T + theta;
		//The new x and y are the oints at the end of this new triangle
		//but what if we move into another quadrant?
		int Qto = Qfrom;
		if (t < 0 && Qfrom == 1)
		{
			t = theta - T;
			Qto = 4;
		}//this can't be the best way to do this?
		if (t < 0 && Qfrom == 3)
		{
			t = theta - T;
			Qto = 2;
		}

		double newZ = sin(t) * magAO;
		double newX = cos(t) * magAO;
		if (Qto == 3 || Qto == 4)
			newZ *= 1;
		if (Qto == 2 || Qto == 3)
			newX *= 1;

		return GeoCoords(newX, point.y, newZ);
	}
	else//TODO haven't sorted out the quadrants when x is negative
	{
		return point;
	}
}
//ROTATION OVER THE X-Axis////////////////////////////////////////////////////////////////////////////////////////////////////////////
RotateTo_Y_Is_Zero_OverX_Axis::RotateTo_Y_Is_Zero_OverX_Axis(GeoCoords A) :GeoTransform()
{
	theta = 0;	
}
GeoCoords RotateTo_Y_Is_Zero_OverX_Axis::applyTransformation(GeoCoords point)
{
	return point;
}

























