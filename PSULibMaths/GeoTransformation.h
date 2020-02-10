#pragma once
#include <GeoPlane.h>
#include <GeoVector.h>
#include <GeoCoords.h>
#include <vector>




using namespace std;

/* I have failed to work with 3d transformations, so for now I am going to go via 0,0,0
and map onto the x and y axes
TODO put it into 3d space properly
*/

class GeoTransform
{
public:
	//virtual GeoCoords applyTransformation(GeoCoords point) { return GeoCoords(0, 0, 0); }; //won't let me =0 it ??? TODO
	GeoTransform() {}
	virtual GeoCoords applyTransformation(GeoCoords point) = 0;
private:
	

};


class TranslateRelativeToOrigin :public GeoTransform
{
private:
	GeoVector V;
public:
	TranslateRelativeToOrigin(GeoCoords p);
	GeoCoords applyTransformation(GeoCoords point) override;
};

class RotateTo_Y_Is_Zero_AboutOrigin :public GeoTransform
{
private:
	double theta; // in radians
public:
	RotateTo_Y_Is_Zero_AboutOrigin(GeoCoords p);
	GeoCoords applyTransformation(GeoCoords point) override;
};
class RotateTo_Z_Is_Zero_AboutOrigin :public GeoTransform
{
private:
	double theta; // in radians
public:
	RotateTo_Z_Is_Zero_AboutOrigin(GeoCoords p);
	GeoCoords applyTransformation(GeoCoords point) override;
};

class RotateTo_Y_Is_Zero_OverX_Axis :public GeoTransform
{
private:
	double theta; // in radians
public:
	RotateTo_Y_Is_Zero_OverX_Axis(GeoCoords p);
	GeoCoords applyTransformation(GeoCoords point) override;
};


class GeoTransformations
{
private:
	vector<GeoTransform*> transformations;
	GeoCoords MapA;
	GeoCoords MapB;
	GeoCoords MapC;
public:
	GeoTransformations();
	~GeoTransformations();
	GeoTransformations(GeoCoords a, GeoCoords b, GeoCoords c);
	GeoCoords applyTransformation(GeoCoords point);
private:	
};


