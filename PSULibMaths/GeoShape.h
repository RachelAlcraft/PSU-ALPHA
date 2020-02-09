#pragma once
#include <GeoTransformation.h>
#include <GeoCoords.h>



class GeoCuboid
{
	public:
		
};

class GeoTripod
{
public://lazy public interface TODO
	GeoCoords anchor;
	GeoCoords axisFar;
	GeoCoords perpFar;
public:
	GeoTripod() {}
	GeoTripod(GeoCoords anchor, GeoCoords axisFar, GeoCoords perpFar);
	GeoTransformation getTransformation(GeoTripod tri, int orientation);
	GeoTripod operator = (GeoTripod const& obj);

};
