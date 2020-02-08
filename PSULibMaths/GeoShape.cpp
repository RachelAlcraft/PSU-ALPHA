#include "GeoShape.h"

GeoTripod::GeoTripod(GeoCoords anchor, GeoCoords axisFar, GeoCoords perpFar)
{
	_anchor = anchor;
	_axisFar = axisFar;
	_perpFar = perpFar;
}

GeoTransformation GeoTripod::getTransformation(GeoTripod tri)
{
	//To describe the transformation all I need is the tripod.
	//To apply the transformation.... well
	
	//Looking at similar triangles in a parallel plane
	//So I am looking at the distance from the perpendicular to the plane
	//And that proportional distance is equal to the propertional distance I move
	//In the direction of the vector on my parallel plane

	// There will be a translation first

	// and then a rotation about a point????
	return GeoTransformation();
}
