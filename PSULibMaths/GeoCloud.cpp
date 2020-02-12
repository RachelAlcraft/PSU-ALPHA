#include "GeoCloud.h"
#include <GeoCoords.h>
#include <GeoVector.h>
#include <set>
#include <map>

/*
TODO
CENTRE OF GEOMETRY
I'm not convinved that a centre of geometry (which I've since discovered is the usual way) is any better.
For a centre of geometry I would then want to find the furtherst point for an axis, and the furtherst othogonal
for another, it seems it would come roughly to the same thing and since the centre of geometry has no topological meaning
doesn't seem intrinsically better. I will experiment with this method and see what I get before deciding.

INSIDE OUTSIDE
For the problem of deciding what is inside and what is outside - there is no way of doing it without the bonds so the GeoCloud would have to become a connected cloud.
Then it would involve finding no particles in certain directons away from bonds. Surface particles could be labelled with a value, and the distance from the surface 
could be recorded for each particle.
*/

GeoCloud::GeoCloud()
{
}

void GeoCloud::addCoords(GeoCoords coord)
{
	_coords.push_back(coord);
}

void GeoCloud::makeTripod(GeoTripod& geo, unsigned int best1, unsigned int best2)
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
			bool inserted = false;
			for (unsigned int k = 0; k < farmag.size(); ++k)
			{
				double dd = farmag[k];
				if (mag > dd && !inserted)
				{//biggest at the front					
					farmag.insert(farmag.begin() + k, mag);
					farcoords.insert(farcoords.begin() + k, pair<GeoCoords, GeoCoords>(_coords[i], _coords[j]));
					inserted = true;
				}
			}
			if (!inserted && farmag.size() < best1)
			{
				farmag.push_back(mag);
				farcoords.push_back(pair<GeoCoords, GeoCoords>(_coords[i], _coords[j]));
			}
			if (farmag.size() > best1)
			{
				farmag.pop_back();
				farcoords.pop_back();
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
		double distance = ortho1.getOrthogonalDistance(_coords[i]);
		bool inserted = false;
		for (unsigned int j = 0; j < farmag2.size(); ++j)
		{
			double dd = farmag2[j];
			if (distance > dd && !inserted)
			{//biggest at the front
				farmag2.insert(farmag2.begin() + j, distance);
				farcoord2.insert(farcoord2.begin() + j, _coords[i]);
				inserted = true;
			}
		}
		if (!inserted && farmag2.size() < best2)
		{
			farmag2.push_back(distance);
			farcoord2.push_back(_coords[i]);
		}
		if (farmag2.size() > best2)
		{
			farmag2.pop_back();
			farcoord2.pop_back();
		}		
	}
	_furthestPoint2 = farcoord2[0];
	//This should fully define what we need for a transformation for this cloud. We have an axis along the furthyest 2 pointsd
	// Orthogonal to that we have the next furthest point
	// Our transformations of clouds will match anchor points along the axis and then rotate to the plane
	geo.A = _furthestPoints1.first;
	geo.B = _furthestPoints1.second;
	geo.C = _furthestPoint2;
}
