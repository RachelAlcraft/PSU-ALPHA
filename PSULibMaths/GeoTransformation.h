#pragma once
#include <GeoPlane.h>
#include <GeoVector.h>
#include <GeoCoords.h>



class GeoTransformation
{
	public:
		GeoCoords RotationPoint;
		GeoVector BaseTranslation;
		GeoPlane BasePlane;
	public:
		GeoTransformation();
		GeoCoords applyTransformation(GeoCoords point);

};

