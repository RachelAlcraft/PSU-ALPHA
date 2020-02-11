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
	double costheta = 2 * magAO2 - magAB2 / (2 * magAO2);
	theta = acos(costheta);// in radians	
}
GeoCoords RotateTo_Y_Is_Zero_AboutOrigin::applyTransformation(GeoCoords point)
{
	return point;
}
RotateTo_Z_Is_Zero_AboutOrigin::RotateTo_Z_Is_Zero_AboutOrigin(GeoCoords A) :GeoTransform()
{
	theta = 0;	
}
GeoCoords RotateTo_Z_Is_Zero_AboutOrigin::applyTransformation(GeoCoords point)
{
	return point;
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























/*GeoCoords GeoTransformations::rotateAboutPoint(GeoCoords point, GeoCoords movePoint, GeoCoords fixPointA, GeoCoords fixPointB)
{
	//It will rotate theta about the perpendicular axis to the plane by a length propertional to the distances of the template to the new from the perpendicular
	GeoVector toCentre = GeoVector(movePoint, fixPointA);
	double mag = toCentre.getMagnitude();
	GeoVector fromCentre = GeoVector(fixPointA, fixPointB) * mag;
	GeoVector moved = toCentre + point;
	moved = moved + fromCentre;
	//We now have 3 points of the template isoceles triangle
	// 2 of the lengths are mag the top is
	double tri_top = moved.getMagnitude();
	//Find theta with the cosine rule
	double mag2 = pow(mag, 2);
	double tri2 = pow(tri_top, 2);
	double costheta = 2*mag2 - tri2 / (2 * mag2);
	double theta = acos(costheta);// in radians
	//find the relative magnitudes to the perpendicular
	GeoVector perp = GeoPlane(toCentre, fromCentre).getPerpendicular();
	double templateMag = perp.getOrthogonalDistance(movePoint);
	double newMag = perp.getOrthogonalDistance(point);
	//move theta about perpendicular which will give newMag
	return point;

}*/

