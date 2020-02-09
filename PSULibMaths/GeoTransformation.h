#pragma once
#include <GeoPlane.h>
#include <GeoVector.h>
#include <GeoCoords.h>



class GeoTransformation
{
	public:
		GeoCoords FixedPointA;
		GeoCoords FixedPointB;
		GeoCoords FixedPointC;
		GeoCoords MovePointA;
		GeoCoords MovePointB;
		GeoCoords MovePointC;
private:


	public:
		GeoTransformation();
		GeoTransformation(GeoCoords,GeoCoords, GeoCoords, GeoCoords, GeoCoords, GeoCoords);
		GeoCoords applyTransformation(GeoCoords point);
private:
	GeoCoords translate(GeoCoords point, GeoCoords movePoint, GeoCoords fixPoint);
	GeoCoords rotateAboutPoint(GeoCoords point, GeoCoords movePoint, GeoCoords fixPointA, GeoCoords fixPointB);
	GeoCoords rotateAboutLine(GeoCoords point, GeoCoords movePoint, GeoCoords fixPointA, GeoCoords fixPointB, GeoCoords fixPointC);		
};

