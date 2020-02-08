#include "GeoCloud.h"
#include <GeoCoords.h>
#include <GeoVector.h>
#include <set>
#include <map>


GeoCloud::GeoCloud()
{
}

void GeoCloud::addCoords(GeoCoords coord)
{
}

GeoTripod GeoCloud::getTripod(unsigned int best1, unsigned int best2)
{
	//First we need to do a distance map and find the furthest points :-(
	vector<double> farmag;
	vector<pair<GeoCoords, GeoCoords >> farcoords;
	
	for (unsigned int i = 0; i < _coords.size(); ++i)
	{
		for (unsigned int j = i+1; j < _coords.size(); ++j)
		{
			GeoVector gv(_coords[i], _coords[j]);
			double mag = gv.getMagnitude();
			bool inserted = true;
			for (unsigned int j = 0; j < farmag.size(); ++j)
			{
				double dd = farmag[j];
				if (mag > dd && !inserted)
				{//biggest at the front
					farmag.insert(farmag.begin() + j, mag);
					farcoords.insert(farcoords.begin() + j, pair<GeoCoords, GeoCoords>(_coords[i], _coords[j]));
				}
				if (!inserted && farmag.size() < best1)
				{
					farmag.insert(farmag.end() + j, mag);
					farcoords.insert(farcoords.end() + j, pair<GeoCoords, GeoCoords>(_coords[i], _coords[j]));
				}
				if (farmag.size() > best1)
				{
					farmag.pop_back();
					farcoords.pop_back();
				}
			}						
		}
	}

	if (farmag.size() > 0)
	{		
		_furthestPoints1.first = farcoords[0].first;
		_furthestPoints1.second = farcoords[0].second;
	}

	//Then we need to find the furtherst piunts on a distance orthogonal at some rotation
	GeoVector ortho1(_furthestPoints1.first, _furthestPoints1.second);
	//we go through all the points and find the one that is furtherst from this vector to make the next orthoginal vector
	vector<double> farmag2;
	vector<GeoCoords> farcoord2;
	double furthest2 = 0;
	for (unsigned int i = 0; i < _coords.size(); ++i)
	{
		double distance = 0;// = GeoVector::getOrthogonalDistance(_furthestPoints1.first, _furthestPoints1.second, _coords[i]);
		bool inserted = true;
		for (unsigned int j = 0; j < farmag2.size(); ++j)
		{
			double dd = farmag2[j];
			if (distance > dd && !inserted)
			{//biggest at the front
				farmag2.insert(farmag2.begin() + j, distance);
				farcoord2.insert(farcoord2.begin() + j, _coords[i]);
			}
			if (!inserted && farmag2.size() < best2)
			{
				farmag2.insert(farmag2.end() + j, distance);
				farcoord2.insert(farcoord2.end() + j, _coords[i]);
			}
			if (farmag2.size() > best2)
			{
				farmag2.pop_back();
				farcoord2.pop_back();
			}
		}
	}
	_furthestPoint2 = farcoord2[0];
	//This should fully define what we need for a transformation for this cloud. We have an axis along the furthyest 2 pointsd
	// Orthogonal to that we have the next furthest point
	// Our transformations of clouds will match anchor points along the axis and then rotate to the plane
	return GeoTripod(_furthestPoints1.first, _furthestPoints1.second, _furthestPoint2);
}
