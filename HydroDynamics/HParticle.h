#pragma once
#ifndef __HPARTICLE_H__
#define __HPARTICLE_H__

#include "fade3d/include_fade3d/Fade_3D.h"

using namespace FADE3D;
using namespace std;

struct HParticle
{
public:
	HParticle(double x, double y, double z) : mass(0.1), velocity(0, 0, 0), coordinates(x, y, z), density(1000), volume(0.1)
	{
	}

	Point3 coordinates;
	Point3 velocity;
	double mass;
	double density;
	double volume;
	vector<double> tetsDensities;
	vector<double> tetsVolumes;
	vector<Point3> tetsVelocities;
	vector<Point3> vectorsB;
	vector<HParticle*> neighbours_points;

};



#endif