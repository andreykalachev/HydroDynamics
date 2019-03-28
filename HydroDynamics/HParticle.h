#pragma once
#ifndef __HPARTICLE_H__
#define __HPARTICLE_H__

#include "fade3d/include_fade3d/Fade_3D.h"

using namespace FADE3D;
using namespace std;

struct HParticle
{
public:
	HParticle(double x, double y, double z) :  coordinates(x, y, z), density(1400), temperature(400)
	{
	}

	Point3 coordinates;
	Point3 velocity;
	double mass;
	double density;
	double volume;
	double temperature;
	vector<double> tetsDensities;
	vector<double> tetsVolumes;
	vector<double> tetsTemperature;
	vector<Point3> tetsVelocities;
	vector<Point3> vectorsB;
	vector<HParticle*> neighbours_points;


	void clear()
	{
		vectorsB.clear();
		neighbours_points.clear();
		tetsVelocities.clear();
		tetsDensities.clear();
		tetsVolumes.clear();
	}
};



#endif