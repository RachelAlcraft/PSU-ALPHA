#include "GeoTransformation.h"
#include <math.h>


GeoTransformation::GeoTransformation()
{
	

}
GeoTransformation::GeoTransformation(GeoCoords fixPointA, GeoCoords fixPointB, GeoCoords fixPointC, GeoCoords movPointA, GeoCoords movPointB, GeoCoords movPointC)
{
	FixedPointA = fixPointA;
	FixedPointB = fixPointB;
	FixedPointC = fixPointC;
	MovePointA = movPointA;
	MovePointB = movPointB;
	MovePointC = movPointC;
}

GeoCoords GeoTransformation::applyTransformation(GeoCoords point)
{
	/*To describe the transformation all I need is the tripod.
	To apply the transformation.... well...*/

	//1) Total freedom, a translation to the anchor point
	GeoCoords gc = translate(point, MovePointA,FixedPointA);

	//2) 1 fixed point, a rotation about a fixed point onto a line
	gc = rotateAboutPoint(gc, MovePointB, FixedPointA,FixedPointB);
	
	
	//3) 2 fixed points, a rotation about a fixed line onto a plane
	gc = rotateAboutLine(gc, MovePointC, FixedPointA, FixedPointB, FixedPointC);
	
	return gc;
}

GeoCoords GeoTransformation::translate(GeoCoords point, GeoCoords movePoint, GeoCoords fixPoint)
{
	GeoVector v = GeoVector(movePoint, fixPoint);
	v = v + point;
	return GeoCoords(v.x, v.y, v.z);
}
GeoCoords GeoTransformation::rotateAboutPoint(GeoCoords point, GeoCoords movePoint, GeoCoords fixPointA, GeoCoords fixPointB)
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

}
GeoCoords GeoTransformation::rotateAboutLine(GeoCoords point, GeoCoords movePoint, GeoCoords fixPointA, GeoCoords fixPointB, GeoCoords fixPointC)
{
	return point;
}