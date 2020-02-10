#include "GeoShape.h"

GeoTripod::GeoTripod(GeoCoords anchor, GeoCoords axisFar, GeoCoords perpFar)
{
	anchor = anchor;
	axisFar = axisFar;
	perpFar = perpFar;
}

GeoTransformations* GeoTripod::getTransformation(GeoTripod tri, int orientation)
{
	/*
	TODO
	I ought to be getting variious different orientations of the transformations
	For now they are all the same
	*/
	if (orientation >0)
	{
		return new GeoTransformations(tri.A, tri.B, tri.C);
	}
	else
	{
		return new GeoTransformations(tri.A, tri.B, tri.C);
	}

}

GeoTripod GeoTripod::operator=(GeoTripod const& obj)
{
	A = obj.A;
	B = obj.B;
	C = obj.C;
	return (GeoTripod(A, B, C));
}
