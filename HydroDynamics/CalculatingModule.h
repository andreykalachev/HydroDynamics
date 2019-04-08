#pragma once
#include <math.h>
#include "Math_operations.h"

//number of particles, new_particles
const int number_of_particles = 50;
const double coef = 1e24;
const double volume = 37.5 * 1e-27 * coef;
const double box_size = cbrt(volume);
const double h = 1e-4;
const double shear_viscosity = 9.0898e-5;
const double bulk_viscosity = 3.0272e-5;
const double bolzmana = 1.38 * 1e-23 * coef;
vector<HParticle*> particles(number_of_particles);
vector<HParticle> particles_image(number_of_particles);
vector<HParticle*> new_particles(number_of_particles);

static bool compare_coords(vector<Tet3*>::value_type& tet, int tet_corner, HParticle* particle)
{
	return tet->getCorner(tet_corner)->x() == particle->coordinates.x() && tet->getCorner(tet_corner)->y() == particle->coordinates.y() &&
		tet->getCorner(tet_corner)->z() == particle->coordinates.z();
}

double calcTetVolume(vector<Tet3*>::value_type& tet)
{
	double tet_volume = (tet->getCorner(0)->x() - tet->getCorner(1)->x()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z())
		+ (tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->z() - tet->getCorner(1)->z()) +
		(tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(0)->z() - tet->getCorner(1)->z()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z());

	return fabs(tet_volume);
}

Point3 calculateVectorB(vector<Tet3*>::value_type& tet, int corner, double tet_volume)
{
	int sign = corner % 2 == 0 ? 1 : -1;

	int _1 = positive_mod((corner + sign * 1), 4);
	int _2 = positive_mod((corner + sign * 2), 4);
	int _3 = positive_mod((corner + sign * 3), 4);

	double x = -(tet->getCorner(_2)->y() - tet->getCorner(_1)->y())*(tet->getCorner(_3)->z() - tet->getCorner(_1)->z()) +
		(tet->getCorner(_3)->y() - tet->getCorner(_1)->y())*(tet->getCorner(_2)->z() - tet->getCorner(_1)->z());

	double y = -(tet->getCorner(_2)->z() - tet->getCorner(_1)->z())*(tet->getCorner(_3)->x() - tet->getCorner(_1)->x()) +
		(tet->getCorner(_3)->z() - tet->getCorner(_1)->z())*(tet->getCorner(_2)->x() - tet->getCorner(_1)->x());

	double z = -(tet->getCorner(_2)->x() - tet->getCorner(_1)->x())*(tet->getCorner(_3)->y() - tet->getCorner(_1)->y()) +
		(tet->getCorner(_3)->x() - tet->getCorner(_1)->x())*(tet->getCorner(_2)->y() - tet->getCorner(_1)->y());

	return Point3(x / tet_volume, y / tet_volume, z / tet_volume);
}

//calculate tetrahedrons' neighbors, volume, density, velocity and vector B
void tetsCount(int index)
{
	//loop for each tetrahedron
	for (auto& tetrahedron : tets)
	{
		for (int corner = 0; corner < 4; corner++)
		{
			//if index particle lie in the nth corner of the tetrahedron
			if (compare_coords(tetrahedron, corner, particles[index]))
			{
				auto new_tet = new Tetrahedron();
				//check all other corners of the tetrahedrons and add if there are any particles in these corners add them to the neighbors
				for (int j = 0; j < particles_image.size(); j++)
				{
					auto particle = &particles_image[j];
					if (tetrahedron->hasVertex(particle->coordinates)) 
					{
						particles[index]->neighbours_points.push_back(particle);
						new_tet->density += particle->density / 4;
						new_tet->temperature += particle->temperature / 4;
						new_tet->increase_velocity(particle->velocity.x() / 4, particle->velocity.y() / 4, particle->velocity.z() / 4);
					}
				}
				new_tet->volume = calcTetVolume(tetrahedron);
				new_tet->vectorB = calculateVectorB(tetrahedron, corner, new_tet->volume);
				particles[index]->tets.push_back(new_tet); break;
			}
		}
	}
}

void calcVolume(int index)
{
	double volume = 0;

	for (int i = 0; i < particles[index]->tets.size(); i++)
	{
		volume += particles[index]->tets[i]->volume;
	}
	particles[index]->volume = volume;
}

void calcTempreture(int index)
{
	auto particle = particles[index];
	double velocity_squared = pow(particle->velocity.x(), 2) + pow(particle->velocity.y(), 2) + pow(particle->velocity.z(), 2);
	new_particles[index]->temperature = particle->density * particle->volume * velocity_squared / (bolzmana * 3);
}

void calcNewDensity(int index)
{
	double sum = 0;
	for (int i = 0; i < particles[index]->tets.size(); i++)
	{
		auto tet = particles[index]->tets[i];
		sum += tet->volume * tet->density * (tet->vectorB.x() * tet->velocity.x() + tet->vectorB.y() * tet->velocity.y() + tet->vectorB.z() * tet->velocity.z());
	}
	new_particles[index]->density = particles[index]->density + h * sum / particles[index]->volume;
}

double calcPressure(double density)
{
	static double a, b, c, d, e, f;
	a = -1.86789e-004;
	b = 1.9878e-001;
	c = -9.6160e001;
	d = 8.2666e004;
	e = -1.5356e006;
	f = 1.1584e-007;

	return f * pow(density, 5) + a * pow(density, 4) + b * pow(density, 3) + c * pow(density, 2) + d * density + e;
}

Point3 calcForce(int index)
{
	auto particle = particles[index];
	double term1_x = 0, term1_y = 0, term1_z = 0, term2_x = 0, term2_y = 0, term2_z = 0;
	particles[index]->mass = particle->volume * particle->density;

	for (int i = 0; i < particle->tets.size(); i++)
	{
		auto tet = particle->tets[i];
		auto pressure = calcPressure(tet->density);

		term1_x += tet->volume * (tet->vectorB.x()*(pressure + pow(tet->velocity.x(), 2)*tet->density)
			+ tet->vectorB.y()*tet->velocity.y()*tet->velocity.x()*tet->density + tet->vectorB.z()*tet->velocity.z()*tet->velocity.x()*tet->density);

		term1_y += tet->volume * (tet->vectorB.x()*tet->velocity.x()*tet->velocity.y()*tet->density
			+ tet->vectorB.y()*(pressure + pow(tet->velocity.y(), 2)*tet->density) + tet->vectorB.z()*tet->velocity.y()*tet->velocity.z()*tet->density);

		term1_z += tet->volume * (tet->vectorB.x()*tet->velocity.z()*tet->velocity.x()*tet->density
			+ tet->vectorB.y()*tet->velocity.z()*tet->velocity.y()*tet->density + tet->vectorB.z()*(pressure + pow(tet->velocity.z(), 2)*tet->density));
	}

	for (int i = 0; i < particle->neighbours_points.size(); i++)
	{
		auto neighbor = particle->neighbours_points[i];

		for (int j = 0; j < neighbor->tets.size(); j++)
		{
			for (int k = 0; k < particle->tets.size(); k++)
			{
				auto n_tet = neighbor->tets[j];
				auto tet = particle->tets[k];

				if (tet->is_equal(*n_tet))
				{
					double scalar_product1 = (tet->velocity.x() - neighbor->velocity.x()) * n_tet->vectorB.x() +
						(tet->velocity.y() - neighbor->velocity.y()) * n_tet->vectorB.y() +
						(tet->velocity.z() - neighbor->velocity.z()) * n_tet->vectorB.z();

					double scalar_product2 = tet->vectorB.x() * n_tet->vectorB.x() +
						tet->vectorB.y() * n_tet->vectorB.y() + tet->vectorB.z() * n_tet->vectorB.z();

					double scalar_product3 = (tet->velocity.x() - neighbor->velocity.x()) * tet->vectorB.x() +
						(tet->velocity.y() - neighbor->velocity.y()) * tet->vectorB.y() +
						(tet->velocity.z() - neighbor->velocity.z()) * tet->vectorB.z();

					term2_x += tet->temperature * tet->volume * (bulk_viscosity * scalar_product1 * tet->vectorB.x() +
						(tet->velocity.x() - neighbor->velocity.x()) * scalar_product2 * shear_viscosity +
						scalar_product3 * n_tet->vectorB.x() * shear_viscosity - 2.0 / 3.0 * scalar_product1 * tet->vectorB.x() * shear_viscosity) / neighbor->temperature;

					term2_y += tet->temperature * tet->volume * (bulk_viscosity * scalar_product1 * tet->vectorB.y() +
						(tet->velocity.y() - neighbor->velocity.y()) * scalar_product2 * shear_viscosity +
						scalar_product3 * n_tet->vectorB.y() * shear_viscosity - 2.0 / 3.0 * scalar_product1 * tet->vectorB.y() * shear_viscosity) / neighbor->temperature;

					term2_z += tet->temperature * tet->volume * (bulk_viscosity * scalar_product1 * tet->vectorB.z() +
						(tet->velocity.z() - neighbor->velocity.z()) * scalar_product2 * shear_viscosity +
						scalar_product3 * n_tet->vectorB.z() * shear_viscosity - 2.0 / 3.0 * scalar_product1 * tet->vectorB.z() * shear_viscosity) / neighbor->temperature;
				}
			}
		}
	}

	return Point3(term1_x + term2_x, term1_y + term2_y, term1_z + term2_z);
}

void calcNewVelocity(int index)
{
	Point3 force = calcForce(index);

	new_particles[index]->velocity = Point3(particles[index]->velocity.x() + h * force.x() / particles[index]->mass,
		particles[index]->velocity.y() + h * force.y() / particles[index]->mass,
		particles[index]->velocity.z() + h * force.z() / particles[index]->mass);
}

void calcNewPosition(int index)
{
	new_particles[index]->coordinates = Point3(particles[index]->coordinates.x() + h * particles[index]->velocity.x(),
		particles[index]->coordinates.y() + h * particles[index]->velocity.y(),
		particles[index]->coordinates.z() + h * particles[index]->velocity.z());
}