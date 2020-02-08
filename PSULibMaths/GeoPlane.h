#pragma once

#include "GeoCoords.h"


class GeoPlane
{
public:
	double getOrthogonalDistance(GeoCoords p);//shortest distance from p to the plane (via the orthogonal vector)
};

