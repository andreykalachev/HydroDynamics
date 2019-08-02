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

	HParticle(double x, double y, double z) : coordinates(x, y, z), mass(0), density(0.844), volume(0), temperature(86.5)
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

	void display_coords(ofstream &file)
	{
		file << fixed << this->coordinates.x() << " \t" << this->coordinates.y() << " \t" << this->coordinates.z() << " \t";
	}

	void display_velocity(ofstream &file)
	{
		file << fixed << this->velocity.x() * 1e2 << " \t" << this->velocity.y() * 1e2 << " \t" << this->velocity.z() * 1e2 << " \t";
	}

	void display_density(ofstream &file)
	{
		file << fixed << this->density * 1.66e3 << " \t";
	}

	void display_velocity_fluctuations(ofstream &file, Point3 *system_velocity)
	{
		file << fixed << (this->velocity.x() - system_velocity->x()) * 1e2 << " \t" << 
			(this->velocity.y() - system_velocity->y()) * 1e2 << " \t" << (this->velocity.z() - system_velocity->z()) * 1e2 << " \t";
	}

	void display_density_fluctuations(ofstream &file, double system_density)
	{
		file << fixed << (this->density - system_density) * 1.66e3 << " \t";
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