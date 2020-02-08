#include "GeoTransformation.h"


GeoTransformation::GeoTransformation()
{

}

GeoCoords GeoTransformation::applyTransformation(GeoCoords point)
{
	return GeoCoords(point.x, point.y, point.z);
}