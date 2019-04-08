#pragma once
#include "fade3d/include_fade3d/Fade_3D.h"

using namespace FADE3D;
using namespace std;

class Tetrahedron
{
public:
	Tetrahedron(double Volume, double Density, double Temperature, Point3 Velocity, Point3 VectorB)
	{
		density = Density;
		volume = Volume;
		temperature = Temperature;
		velocity = Velocity;
		vectorB = VectorB;
	}

	Tetrahedron(): density(0), volume(0), temperature(0), velocity(0,0,0), vectorB(0,0,0)
	{
	}

	double density;
	double volume;
	double temperature;
	Point3 velocity;
	Point3 vectorB;

	void increase_velocity(double x, double y, double z)
	{
		velocity = Point3(velocity.x() + x, velocity.y() + y, velocity.z() + z);
	}

	bool is_equal(Tetrahedron tet)
	{
		return tet.density == density && tet.velocity.x() == velocity.x() && tet.velocity.y() == velocity.y() &&
			tet.velocity.z() == velocity.z() && tet.volume == volume;
	}

};