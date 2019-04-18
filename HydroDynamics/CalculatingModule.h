#pragma once
#include <math.h>
#include "Math_operations.h"

//current iteration
int iteration = 0;
const int number_of_particles = 10;
const double volume = 37.5 * 1e-27 * pow(1e9,3);
const double box_size = cbrt(volume);
double time_step = 1e-12 * 1e12;
const double shear_viscosity = 9.0898e-5 / 1.66e-3;
const double bulk_viscosity = 3.0272e-5 / 1.66e-3;
const double bolzmana = 1.38064852 * 1e-23 / 1.66e-18;
//total system velocity
Point3 system_velocity = Point3(0, 0, 0);
vector<HParticle*> particles(number_of_particles);
//images of particles to clrate periodic conditions
vector<HParticle> particles_image(0);
//values for next iteration are stored in this vector
vector<HParticle*> new_particles(number_of_particles);

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

//copy all tetrahedrons from particles to all their images
void copyTets()
{
	for (int i = 0, n = number_of_particles; i < n; i++)
	{
		for (int j = 0; j < 27; j++)
		{
			particles_image[j + i * 27].tets = particles[i]->tets;
		}
	}
}

//calculate particles' tetrahedrons, neighbors, volume, mass 
void tetsCount()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		//loop for each tetrahedron
		for (auto& tetrahedron : tets)
		{
			for (int corner = 0; corner < 4; corner++)
			{
				//if index particle lie in the nth corner of the tetrahedron
				if (compare_coords(tetrahedron, corner, particles[i]))
				{
					auto new_tet = new Tetrahedron();
					//check all other corners of the tetrahedrons and add if there are any particles in these corners add them to the neighbors
					for (int j = 0; j < particles_image.size(); j++)
					{
						auto particle = &particles_image[j];
						//number of vertexes that we have taken into account (we need only 4 for each tetrahedron)
						auto vertex_count = 0;
						if (tetrahedron->hasVertex(particle->coordinates))
						{
							vertex_count++;
							//check that we don't add the particle itself to its neighbours and that we don't add one neighbours more than one time
							if (!(particle->coordinates == particles[i]->coordinates) ) 
									particles[i]->neighbours_points.push_back(particle);
							new_tet->density += particle->density / 4;
							new_tet->temperature += particle->temperature / 4;
							new_tet->increase_velocity(particle->velocity.x() / 4, particle->velocity.y() / 4, particle->velocity.z() / 4);
							//we exit the loop when we have 4 vertexes
							if (vertex_count == 4) break;
						}
					}
					new_tet->volume = calcTetVolume(tetrahedron);
					new_tet->vectorB = calculateVectorB(tetrahedron, corner, new_tet->volume);
					particles[i]->tets.push_back(new_tet); 
					particles[i]->volume += new_tet->volume;
					break;
				}
			}
		}
		particles[i]->mass = particles[i]->volume * particles[i]->density;
		if (iteration == 0) particles[i]->momentum = Point3(particles[i]->velocity.x() * particles[i]->mass, particles[i]->velocity.y() * particles[i]->mass, particles[i]->velocity.z() * particles[i]->mass);
	}
	copyTets();
}

void calcTempreture(int index)
{
	new_particles[index]->temperature = particles[index]->density * particles[index]->volume * pow(particles[index]->velocity_absolute, 2) / (bolzmana * 3);
}

void calcNewDensity(int index)
{
	double sum = 0;
	for (int i = 0; i < particles[index]->tets.size(); i++)
	{
		auto tet = particles[index]->tets[i];
		sum += tet->volume * tet->density * (tet->vectorB.x() * tet->velocity.x() + tet->vectorB.y() * tet->velocity.y() + tet->vectorB.z() * tet->velocity.z());
	}
	new_particles[index]->density = particles[index]->density + time_step * sum / particles[index]->volume;
}

double calcPressure(double density)
{
	static double a, b, c, d, e, f;
	a = -1.86789e-4;
	b = 1.9878e-1;
	c = -9.6160e1;
	d = 8.2666e4;
	e = -1.5356e6;
	f = 1.1584e-7;

	return (f * pow(density, 5) + a * pow(density, 4) + b * pow(density, 3) + c * pow(density, 2) + d * density + e) / 1.66e9;
}

Point3 calcMomentum(int index)
{
	auto particle = particles[index];
	double term1_x = 0, term1_y = 0, term1_z = 0, term2_x = 0, term2_y = 0, term2_z = 0;

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
	double result_x = (term1_x + term2_x);
	double result_y = (term1_y + term2_y);
	double result_z = (term1_z + term2_z);

	new_particles[index]->momentum =
		Point3(particles[index]->momentum.x() + result_x, particles[index]->momentum.y() + result_y, particles[index]->momentum.z() + result_z);
	new_particles[index]->momentum_absolute = calculate_absolute_value(new_particles[index]->momentum);
	return Point3(result_x, result_y, result_z);
}

void calcNewVelocity(int index)
{
	auto momentum = calcMomentum(index);

	new_particles[index]->velocity = Point3(particles[index]->velocity.x() + time_step * momentum.x() / particles[index]->mass,
		particles[index]->velocity.y() + time_step * momentum.y() / particles[index]->mass,
		particles[index]->velocity.z() + time_step * momentum.z() / particles[index]->mass);
	new_particles[index]->velocity_absolute = calculate_absolute_value(new_particles[index]->velocity);
	system_velocity = calculate_sum(system_velocity, new_particles[index]->velocity);
}

void calcNewPosition(int index)
{
	new_particles[index]->coordinates = Point3(particles[index]->coordinates.x() + time_step * particles[index]->velocity.x(),
		particles[index]->coordinates.y() + time_step * particles[index]->velocity.y(),
		particles[index]->coordinates.z() + time_step * particles[index]->velocity.z());
}