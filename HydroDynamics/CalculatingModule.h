#pragma once
#include <math.h>
#include "Math_operations.h"
#include "Tetrahedron.h"

//current iteration
int iteration = 0;
vector<Tet3*> delaunay_tets;
const int number_of_particles = 10;
const double volume = 37.5;
const double box_size = cbrt(volume);
double time_step = 1;
double elapsed_time = 0;
const double shear_viscosity = 9.0898 / 166;
const double bulk_viscosity = 3.0272 / 166;
const double b_const = 1.38064852 / 160000 ;
//total system velocity
auto system_velocity = new Point3(0, 0, 0);
vector<HParticle*> particles(number_of_particles);
//images of particles to clrate periodic conditions
vector<HParticle> particles_image(0);
//values for next iteration are stored in this vector
vector<HParticle*> new_particles(number_of_particles);
//random for thermal fluctiations
default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);

double calcTetVolume(vector<Tet3*>::value_type& tet)
{
	double tet_volume = (tet->getCorner(0)->x() - tet->getCorner(1)->x()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z())
		+ (tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->z() - tet->getCorner(1)->z()) +
		(tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(0)->z() - tet->getCorner(1)->z()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z());

	return fabs(tet_volume / 6.0);
}

Point3 calculateVectorB(vector<Tet3*>::value_type& tet, int corner, double tet_volume)
{
	auto _0 = *tet->getCorner(corner);
	auto _1 = *tet->getCorner(positive_mod((corner + 1), 4));
	auto _2 = *tet->getCorner(positive_mod((corner + 2), 4));
	auto _3 = *tet->getCorner(positive_mod((corner + 3), 4));

	auto vector_b = _0 - (_1 + _2 + _3) / 3.0;
	auto s = calculate_area(_0, _1, _2) + calculate_area(_0, _2, _3) + calculate_area(_0, _1, _3) + calculate_area(_1, _2, _3);
	return vector_b / (2 * s);
}

//calculate particles' tetrahedrons, neighbors, volume, mass 
void analyzeTets()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		double particle_volume = 0;
		//loop for each tetrahedron
		for (auto tetrahedron : delaunay_tets)
		{
			for (int corner = 0; corner < 4; corner++)
			{
				//if index particle lie in the nth corner of the tetrahedron
				if (compare_coords(tetrahedron, corner, particles[i]))
				{
					//number of vertexes that we have taken into account (we need only 4 for each tetrahedron)
					auto vertex_count = 0;
					auto new_tet = Tetrahedron(calcTetVolume(tetrahedron), tetrahedron->getCircumcenter());
					new_tet.vectorB = calculateVectorB(tetrahedron, corner, new_tet.volume);
					new_tet.gaussMatrix.generateRandom(generator, distribution);
					//check all other corners of the tetrahedrons and add if there are any particles in these corners add them to the neighbors
					for (int j = 0; j < particles_image.size(); j++)
					{
						auto particle = &particles_image[j];
						auto image_corner = 0;
						if (has_corner(tetrahedron, particle, &image_corner))
						{
							//check that we don't add the particle itself to its neighbours and that we don't add one neighbours more than one time
							if (particles_image[j].coordinates != particles[i]->coordinates)
							{
								if (!isInside(particles[i]->neighbours_points, particle)) particles[i]->neighbours_points.push_back(particle);

								auto tt = Tetrahedron(new_tet.volume, new_tet.circumcenter, calculateVectorB(tetrahedron, image_corner, new_tet.volume));
								if (!isTetInside(particle->tets, tt)) particle->tets.push_back(tt);
							}
							new_tet.density += particle->density / 4.0;
							new_tet.temperature += particle->temperature / 4.0;
							new_tet.velocity += particle->velocity / 4.0;
							if (++vertex_count == 4) break;
						}
					}
					particle_volume += new_tet.volume;
					particles[i]->tets.push_back(new_tet);
					break;
				}
			}
		}
		particles[i]->volume = particle_volume;
		particles[i]->mass = particles[i]->volume * particles[i]->density;
	}
}

double calcPressure(double density, Point3 velocity)
{
	static double a, b, c, d, e, f;
	a = -1.86789e-4;
	b = 1.9878e-1;
	c = -9.6160e1;
	d = 8.2666e4;
	e = -1.5356e6;
	f = 1.1584e-7;

	density *= 1.66e3;
	velocity *= 1000;

	/*return (f * pow(density, 5) + a * pow(density, 4) + b * pow(density, 3) +
		c * pow(density, 2) + d * density + e) / 1.66e9;*/

	return density * pow(calculate_absolute_value(velocity), 2) / 3.0 / 1.66e9;
}

void calcTempreture()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->temperature = calcPressure(particles[i]->density, particles[i]->velocity) * 1e6 / particles[i]->density / 208.13;
	}
}

void calcNewDensity()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		double sum = 0;
		for (const auto tet : particles[i]->tets)
		{
			sum += tet.volume * (tet.vectorB * tet.velocity) * tet.density;
		}
		new_particles[i]->density = particles[i]->density + time_step * sum / particles[i]->volume;
	}
}

void calcNewPosition()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->coordinates = (particles[i]->coordinates + particles[i]->velocity * time_step / 2.0) % box_size;
	}
}

Point3 calcForce1(int index)
{
	auto term1 = Point3(0, 0, 0), term2 = Point3(0, 0, 0), thermal_fluct = Point3(0, 0, 0);

	for (auto tet : particles[index]->tets)
	{
		term1 += tet.volume *  (calcPressure(tet.density, particles[index]->velocity) * tet.vectorB +
			particles[index]->density * tet.vectorB * tet.velocity * (4 * tet.velocity - particles[index]->velocity) / 4.0);

		thermal_fluct += tet.vectorB * (sqrt(4 * b_const*tet.temperature*shear_viscosity*tet.volume) *
			((tet.gaussMatrix + tet.gaussMatrix.getTransported()) / 2.0 - multiply_by_identity_matrix(tet.gaussMatrix.getDiagonalSum() / 3)) +
			sqrt(3 * b_const*tet.temperature*bulk_viscosity*tet.volume) * multiply_by_identity_matrix(tet.gaussMatrix.getDiagonalSum() / 3)) / (tet.volume * sqrt(time_step / 2));

		for (auto neighbor : particles[index]->neighbours_points)
		{
			for (const auto n_tet : neighbor->tets)
			{
				if (tet.has_equal_circumcenter(n_tet))
				{
					const auto scalar_product1 = (tet.velocity - neighbor->velocity) * n_tet.vectorB;
					const auto scalar_product3 = (tet.velocity - neighbor->velocity) * tet.vectorB;

					term2 += tet.temperature * tet.volume / neighbor->temperature *
						(bulk_viscosity * scalar_product1 * tet.vectorB +
							shear_viscosity * (scalar_product3 * n_tet.vectorB - 2.0 / 3.0 * scalar_product1 * tet.vectorB) +
							((4 * tet.velocity - particles[index]->velocity) / 4.0 - n_tet.velocity) * n_tet.vectorB * tet.vectorB);
					break;
				}
			}
		}
	}

	return (term1 + term2 + thermal_fluct) / particles[index]->mass;
}

double calcForce2(int index)
{
	double term1 = 0, term2 = 0;

	for (auto tet : particles[index]->tets)
	{
		term1 += tet.volume * particles[index]->density * tet.vectorB * tet.vectorB;

		for (auto neighbor : particles[index]->neighbours_points)
		{
			for (const auto n_tet : neighbor->tets)
			{
				if (tet.has_equal_circumcenter(n_tet))
				{
					term2 += shear_viscosity * n_tet.vectorB * tet.vectorB * tet.temperature * tet.volume / neighbor->temperature;
					break;
				}
			}
		}
	}
	return  term1 / (4 * particles[index]->mass) + term2;
}

void substractAvgVel()
{
	*system_velocity = Point3(0, 0, 0);
	for (int i = 0; i < number_of_particles; i++)
	{
		*system_velocity += particles[i]->velocity;
	}
	*system_velocity /= number_of_particles;
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->velocity -= *system_velocity;
	}
}

void calcNewVelocityL()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->velocity = particles[i]->velocity + calcForce1(i) * time_step / 2;
	}
	substractAvgVel();
}

void calcNewVelocityO()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->velocity = particles[i]->velocity * exp(calcForce2(i) * time_step);
	}
	substractAvgVel();
}

void lod()
{
	calcNewVelocityL();
	calcNewPosition();
	calcNewVelocityO();
	calcNewPosition();
	calcNewVelocityL();
}