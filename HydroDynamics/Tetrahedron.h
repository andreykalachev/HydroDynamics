#pragma once
#include "fade3d/include_fade3d/Fade_3D.h"

using namespace FADE3D;
using namespace std;

class Tetrahedron
{
public:
	Tetrahedron(double Volume, Point3 Circumcenter, Point3 VectorB): density(0), temperature(0), velocity(0, 0, 0), vectorB(0,0,0)
	{
		volume = Volume;
		circumcenter = Circumcenter;
		vectorB = VectorB;
	}

	Tetrahedron(double Volume, Point3 Circumcenter) : density(0), temperature(0), velocity(0, 0, 0), vectorB(0,0,0), circumcenter(0,0,0)
	{
		volume = Volume;
		circumcenter = Circumcenter;
	}

	Tetrahedron(): density(0), volume(0), temperature(0), velocity(0,0,0), vectorB(0, 0, 0), circumcenter(0, 0, 0)
	{
	}

	Point3 circumcenter;
	double density;
	double volume;
	double temperature;
	Point3 velocity;
	Point3 vectorB;

	bool is_equal(Tetrahedron tet)
	{
		return tet.circumcenter == this->circumcenter;
	}

	Tetrahedron get_copy(Point3 shift)
	{
		auto new_tet = *this;
		new_tet.circumcenter = Point3 (new_tet.circumcenter.x() + shift.x(), new_tet.circumcenter.y() + shift.y(), new_tet.circumcenter.z() + shift.z());
		return new_tet;
	}

};