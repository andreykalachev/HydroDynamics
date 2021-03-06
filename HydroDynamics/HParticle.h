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
		tets = vector<Tetrahedron>();
		neighbours_points = vector<HParticle*>();
	}

	HParticle(double x, double y, double z) : coordinates(x, y, z), mass(0), density(996.3 / 1.66e3), volume(0), temperature(0)
	{
		tets = vector<Tetrahedron>();
		neighbours_points = vector<HParticle*>();
	}

	Point3 coordinates;
	Point3 velocity;
	double mass;
	double density;
	double volume;
	double temperature;
	vector<Tetrahedron> tets;
	vector<HParticle*> neighbours_points;

	void clear()
	{
		tets = vector<Tetrahedron>();
		neighbours_points = vector<HParticle*>();
	}

	//return a copy with new coordinates
	HParticle get_copy(double x, double y, double z)
	{
		auto new_particle = *this;
		new_particle.coordinates.init(Point3(x, y, z));
		return new_particle;
	}

	//return a copy with new coordinates
	HParticle get_copy(Point3 new_coordinates)
	{
		auto new_particle = *this;
		new_particle.coordinates.init(new_coordinates);
		return new_particle;
	}

	void copy(HParticle *particle)
	{
		this->density = particle->density;
		this->velocity = particle->velocity;
		this->coordinates = particle->coordinates;
	}

	void display(ofstream &file)
	{
		file << fixed;
		file << "coords: (" << this->coordinates.x() << ", " << this->coordinates.y() << ", " << this->coordinates.z() << ")  \t";
		file << scientific;
		file << " mass:" << this->mass << " \t" << " temp:" << this->temperature << " \t" << "density:" << this->density << " \t"
			<< "vel: (" << this->velocity.x() << ", " << this->velocity.y() << ", " << this->velocity.z() << ")"
			<< endl << "---------" << endl;
	}

	void display_for_plot(ofstream &file)
	{
		//file << scientific << this->coordinates.x() << " \t" << this->coordinates.y() << " \t" << this->coordinates.z() << " \t";
		file << scientific << this->velocity.x() << " \t" << this->velocity.y() << " \t" << this->velocity.z() << " \t";
		//file << scientific << this->velocity.x() << " \t" << this->velocity.y() << " \t" << this->velocity.z() << " \t" << this->temperature << " \t";

	}
};


static bool isInside(vector<HParticle*> particles, HParticle *p)
{
	for (auto particle : particles)
	{
		if (particle == p) return true;
	}
	return false;
}
