#pragma once
#include <math.h>
#include "Math_operations.h"

//current iteration
int iteration = 0;
vector<Tet3*> tets;
const int number_of_particles = 10;
const double volume = 37.5;
const double box_size = cbrt(volume);
double time_step = 0.02;
const double shear_viscosity = 9.0898 / 166;
const double bulk_viscosity = 3.0272 / 166;
const double bolzmana = 1.38064852 / 1.66e5;
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

	return fabs(tet_volume / 6);
}

Point3 calculateVectorB(vector<Tet3*>::value_type& tet, int corner, double tet_volume)
{
	int sign = corner % 2 == 0 ? 1 : -1;

	auto _1 = tet->getCorner(positive_mod((corner + sign * 1), 4));
	auto _2 = tet->getCorner(positive_mod((corner + sign * 2), 4));
	auto _3 = tet->getCorner(positive_mod((corner + sign * 3), 4));

	auto x = -(_2->y() - _1->y())*(_3->z() - _1->z()) + (_3->y() - _1->y())*(_2->z() - _1->z());
	auto y = -(_2->z() - _1->z())*(_3->x() - _1->x()) + (_3->z() - _1->z())*(_2->x() - _1->x());
	auto z = -(_2->x() - _1->x())*(_3->y() - _1->y()) + (_3->x() - _1->x())*(_2->y() - _1->y());

	auto vector_b = Point3(x, y, z) / (6 * tet_volume);
	return vector_b / calculate_absolute_value(vector_b) * distance(*_1, (*_2 + *_3) / 2) * distance(*_2, *_3) / 2;
}

//calculate particles' tetrahedrons, neighbors, volume, mass 
void analyzeTets()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		double particle_volume = 0;
		//loop for each tetrahedron
		for (auto& tetrahedron : tets)
		{
			for (int corner = 0; corner < 4; corner++)
			{
				//if index particle lie in the nth corner of the tetrahedron
				if (compare_coords(tetrahedron, corner, particles[i]))
				{
					//number of vertexes that we have taken into account (we need only 4 for each tetrahedron)
					auto vertex_count = 0;
					auto new_tet = new Tetrahedron();
					new_tet->volume = calcTetVolume(tetrahedron);
					new_tet->circumcenter = tetrahedron->getCircumcenter();
					new_tet->vectorB = calculateVectorB(tetrahedron, corner, new_tet->volume);
					//check all other corners of the tetrahedrons and add if there are any particles in these corners add them to the neighbors
					for (int j = 0; j < particles_image.size(); j++)
					{
						auto particle = &particles_image[j];
						auto image_corner = 0;
						if (has_corner(tetrahedron, particle, &image_corner))
						{
							if (vertex_count == 4) break;
							vertex_count++;
							//check that we don't add the particle itself to its neighbours and that we don't add one neighbours more than one time
							if (particle->coordinates != particles[i]->coordinates && !isInside(particles[i]->neighbours_points, particle))
								particles[i]->neighbours_points.push_back(particle);
							new_tet->density += particle->density / 4;
							new_tet->temperature += particle->temperature / 4;
							new_tet->velocity += particle->velocity / 4;
							//we exit the loop when we have 4 vertexes

							if (image_corner == corner) continue;
							auto tet = Tetrahedron();
							tet.volume = new_tet->volume;
							tet.circumcenter = new_tet->circumcenter;
							tet.vectorB = calculateVectorB(tetrahedron, image_corner, new_tet->volume);
							particle->tets.push_back(tet);
						}
					}
					particles[i]->tets.push_back(*new_tet); 
					particle_volume += new_tet->volume;
					break;
				}
			}
		}
		particles[i]->volume = particle_volume;
		particles[i]->mass = particles[i]->volume * particles[i]->density;
		if (iteration == 0) particles[i]->momentum = particles[i]->velocity * particles[i]->mass;
		//copy missing data from particle to all it's images
		for (int j = 0; j < 27; j++) {
			particles_image[j + i * 27].mass = particles[i]->mass;
			particles_image[j + i * 27].volume = particles[i]->volume;
		}
	}
}

void calcTempreture(int index)
{
	new_particles[index]->temperature = new_particles[index]->density * particles[index]->volume * pow(new_particles[index]->velocity_absolute, 2) / (bolzmana * 3);
}

void calcNewDensity(int index)
{
	double sum = 0;
	for (auto tet : particles[index]->tets)
	{
		sum += tet.volume * (tet.vectorB * tet.velocity) * tet.density;
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

Point3 calcForce(int index)
{
	auto particle = particles[index];
	auto term1 = Point3(0,0,0), term2 = Point3(0, 0, 0);

	for (auto tet : particle->tets)
	{
		term1 += tet.volume * (calcPressure(tet.density) * tet.vectorB + tet.density * (tet.vectorB * tet.velocity) * tet.velocity);

		for (auto neighbor : particle->neighbours_points)
		{
			for (auto n_tet : neighbor->tets)
			{
				if (tet.is_equal(n_tet))
				{
					const auto scalar_product1 = (tet.velocity - neighbor->velocity) * n_tet.vectorB;
					const auto scalar_product2 = tet.vectorB * n_tet.vectorB;
					const auto scalar_product3 = (tet.velocity - neighbor->velocity) * tet.vectorB;

					term2 += tet.temperature * tet.volume / neighbor->temperature *
						(bulk_viscosity * scalar_product1 * tet.vectorB +
							shear_viscosity * ((tet.velocity - neighbor->velocity) * scalar_product2 +
								scalar_product3 * n_tet.vectorB - 2.0 / 3.0 * scalar_product1 * tet.vectorB));
					break;
				}
			}
		}
	}

	auto result = term1 + term2;

	new_particles[index]->momentum = particle->momentum + result * time_step;
	new_particles[index]->momentum_absolute = calculate_absolute_value(new_particles[index]->momentum);
	return result;
}

void calcNewVelocity(int index)
{
	new_particles[index]->velocity = particles[index]->velocity + calcForce(index) * time_step / particles[index]->mass;
	new_particles[index]->velocity_absolute = calculate_absolute_value(new_particles[index]->velocity);
	system_velocity += new_particles[index]->velocity;
}

void calcNewPosition(int index)
{
	new_particles[index]->coordinates = particles[index]->coordinates + (particles[index]->velocity + new_particles[index]->velocity) * time_step / 2;
}