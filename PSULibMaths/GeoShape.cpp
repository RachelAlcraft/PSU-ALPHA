#include "GeoShape.h"

GeoTripod::GeoTripod(GeoCoords anchor, GeoCoords axisFar, GeoCoords perpFar)
{
	anchor = anchor;
	axisFar = axisFar;
	perpFar = perpFar;
}

GeoTransformation GeoTripod::getTransformation(GeoTripod tri, int orientation)
{
	GeoCoords rotationPoint;
	GeoVector baseTrans;
	GeoPlane basePlane;
	/*To describe the transformation all I need is the tripod.
	To apply the transformation.... well
	
	Looking at similar triangles in a parallel plane	
	So I am looking at the distance from the perpendicular to the plane
	And that proportional distance is equal to the propertional distance I move
	In the direction of the vector on my parallel plane	 */

	//1) Total freedom, a translation to the anchor point
	GeoVector trans1 = GeoVector(anchor, tri.anchor);
	GeoVector trans2 = GeoVector(anchor, tri.axisFar);
	GeoVector trans3 = GeoVector(axisFar, tri.anchor);
	GeoVector trans4 = GeoVector(axisFar, tri.axisFar);
	

	//2) 1 fixed point, a rotation about a fixed point
	GeoVector rotPoint1a = GeoVector(tri.axisFar, tri.anchor);
	double mag1 = rotPoint1a.getMagnitude();
	GeoVector rotPoint1b = GeoVector(anchor, axisFar) * mag1;
	GeoVector rotPoint1c = rotPoint1a + rotPoint1b;

	// etc * 4
	
	//3) 2 fixed points, a rotation about a fixed line (and so on over the dimensions :-))
	
	
	//There are 4 possible orientations, each end and then a rotation of 180 degrees	
	if (orientation == 1)	
		return GeoTransformation(anchor,axisFar,perpFar, tri.anchor, tri.axisFar, tri.perpFar);	
	else if (orientation == 2)	
		return GeoTransformation(anchor, axisFar, perpFar, tri.anchor, tri.axisFar, tri.perpFar);
	else if (orientation == 3)
		return GeoTransformation(anchor, axisFar, perpFar, tri.anchor, tri.axisFar, tri.perpFar);
	else
		return GeoTransformation(anchor, axisFar, perpFar, tri.anchor, tri.axisFar, tri.perpFar);

	//TODO currently these are all the same!!!
}

GeoTripod GeoTripod::operator=(GeoTripod const& obj)
{
	anchor = obj.anchor;
	axisFar = obj.axisFar;
	perpFar = obj.perpFar;
	return (GeoTripod(anchor, axisFar, perpFar));
}
