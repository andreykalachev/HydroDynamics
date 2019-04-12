#pragma once
#include "fade3d/include_fade3d/Fade_3D.h"
#include "Tetrahedron.h"

using namespace FADE3D;
using namespace std;

class HParticle
{
public:

	HParticle()
	{
		
	}

	HParticle(double x, double y, double z) : coordinates(x, y, z), mass(0), density(996.3 / 1.66e3), volume(0), temperature(300)
	{
	}

	Point3 coordinates;
	Point3 velocity;
	Point3 momentum;
	double mass;
	double density;
	double volume;
	double temperature;
	double velocity_absolute;
	double momentum_absolute;
	vector<Tetrahedron*> tets;
	vector<HParticle*> neighbours_points;

	void clear()
	{
		volume = 0;
		tets.clear();
		neighbours_points.clear();
	}

	//return a copy with new coordinates
	HParticle copy(double x, double y, double z)
	{
		auto new_particle = *this;
		new_particle.coordinates = Point3 (x,y,z);
		return new_particle;
	}

	//return a copy with new coordinates
	HParticle copy(Point3 new_coordinates)
	{
		auto new_particle = *this;
		new_particle.coordinates.init(new_coordinates);
		return new_particle;
	}

	void display()
	{
		cout << "coords : ( " << this->coordinates.x() << ", " << this->coordinates.y() << ", " << this->coordinates.z() << " )   \t"
			<< " vel : " << this->velocity_absolute << " \t" << " mom : " << this->momentum_absolute << "   \t"
			<< " mass :" << this->mass << " \t" << " temp:" << this->temperature << " \t" << "density:" << this->density << endl
			<< "---------" << endl;
	}
};
