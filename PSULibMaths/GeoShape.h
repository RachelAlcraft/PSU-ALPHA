#pragma once
#include <GeoTransformation.h>
#include <GeoCoords.h>



class GeoCuboid
{
	public:
		
};

class GeoTripod
{
private:
	GeoCoords _anchor;
	GeoCoords _axisFar;
	GeoCoords _perpFar;
public:
	GeoTripod(GeoCoords anchor, GeoCoords axisFar, GeoCoords perpFar);
	GeoTransformation getTransformation(GeoTripod tri);

};
