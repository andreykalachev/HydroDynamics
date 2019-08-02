#pragma once
#include <math.h>
#include "Math_operations.h"
#include "Tetrahedron.h"

/*
 * 1 amu	=	1.6605e-27 kg
 * 1 å		=	1e-10 m
 * 1 ps		=	1e-12 s
 * 1 pa		=	6.02e-8 reduced units
 * 1 amt	=	6.1e-3 reduced units
 * 1 m/s	=	1e-2 å/ps
 */

//current iteration
int iteration = 0;
//for triangulation
vector<Tet3*> delaunay_tets;
const int number_of_particles = 27;
//volume of the investigated box
const double volume = 1000 * number_of_particles;
//total volume of all particles (more than box volume because particles overlap)
double total_particles_volume = 0;
const double box_length = cbrt(volume);
double time_step = 0.05;
double elapsed_time = 0;
const double shear_viscosity = 5.47;
const double bulk_viscosity = 1.82;
//boltzmann constant
const double b_const = 0.8314;
//total system velocity
auto system_velocity = new Point3(0, 0, 0);
double system_density = 0;
vector<HParticle*> particles(number_of_particles);
//images of particles to create periodic conditions
vector<HParticle> particles_image(0);
//particles' density on the next step
vector<double> new_density(number_of_particles);
//random for thermal fluctuations
default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);

double calc_tet_volume(vector<Tet3*>::value_type& tet)
{
	double tet_volume = (tet->getCorner(0)->x() - tet->getCorner(1)->x()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z())
		+ (tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->z() - tet->getCorner(1)->z()) +
		(tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(0)->z() - tet->getCorner(1)->z()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z());

	return fabs(tet_volume / 6.0);
}

Point3 calculate_vector_b(vector<Tet3*>::value_type& tet, int corner, double tet_volume)
{
	auto _0 = *tet->getCorner(corner);
	auto _1 = *tet->getCorner(positive_mod((corner + 1), 4));
	auto _2 = *tet->getCorner(positive_mod((corner + 2), 4));
	auto _3 = *tet->getCorner(positive_mod((corner + 3), 4));

	auto vector_b = _0 - (_1 + _2 + _3) / 3.0;
	auto s = calculate_area(_0, _1, _2) + calculate_area(_0, _2, _3) + calculate_area(_0, _1, _3) + calculate_area(_1, _2, _3);
	return vector_b / (2 * s);
}

/*
 * analyze which tetrahedrons compose each particle
 * calculate all important parameters for these tetrahedrons
 * for each particle create list of its neighbors
 */
void analyze_tets()
{
	total_particles_volume = 0;
	for (int i = 0; i < number_of_particles; i++)
	{
		system_density += particles[i]->density / number_of_particles;
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
					auto new_tet = Tetrahedron(calc_tet_volume(tetrahedron), tetrahedron->getCircumcenter());
					new_tet.vectorB = calculate_vector_b(tetrahedron, corner, new_tet.volume);
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

								auto tt = Tetrahedron(new_tet.volume, new_tet.circumcenter, calculate_vector_b(tetrahedron, image_corner, new_tet.volume));
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
		total_particles_volume += particles[i]->volume;
		particles[i]->mass = particles[i]->volume * particles[i]->density;
	}
}

//https://doi.org/10.1063/1.4738763
double calc_pressure(double density)
{
	auto isothermal_compressibility = 1.7e4 / 6.1e-3;
	auto thermal_expansivity = 3.6e3;
	auto initial_density = 0.844;
	auto temperature = 86.5;
	auto result = density / initial_density / isothermal_compressibility + thermal_expansivity * temperature / isothermal_compressibility;

	// 1 amt = 6.1e-3 reduced units
	//in theory pressure should be 66atm

	//at the moment is set to constant
	result = 66 * 6.1e-3;

	return result;
}

//https://doi.org/10.1063/1.4738763
void calc_tempreture()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		//at the moment is set to constant
		particles[i]->temperature = 86.5;
	}
}

void calc_new_density()
{
	system_density = 0;
	for (int i = 0; i < number_of_particles; i++)
	{
		double sum = 0;
		for (const auto tet : particles[i]->tets)
		{
			sum += tet.volume * (tet.vectorB * tet.velocity) * tet.density;
		}
		new_density[i] = particles[i]->density + time_step * sum / particles[i]->volume;
		system_density += new_density[i] * particles[i]->volume / total_particles_volume;
	}
}

void calc_new_position()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->coordinates = (particles[i]->coordinates + particles[i]->velocity * time_step / 2.0) % box_length;
	}
}

Point3 calc_force1(int index)
{
	auto term1 = Point3(0, 0, 0), term2 = Point3(0, 0, 0), thermal_fluct = Point3(0, 0, 0), p1 = Point3(0, 0, 0), p2 = Point3(0, 0, 0);

	for (auto tet : particles[index]->tets)
	{
		term1 += tet.volume * (calc_pressure(tet.density) * tet.vectorB +
			particles[index]->density * tet.vectorB * tet.velocity * (4 * tet.velocity - particles[index]->velocity) / 4.0);

		thermal_fluct += tet.vectorB * (sqrt(4 * b_const*tet.temperature*shear_viscosity*tet.volume) *
			((tet.gaussMatrix + tet.gaussMatrix.getTransported()) / 2.0 - multiply_by_identity_matrix(tet.gaussMatrix.getDiagonalSum() / 3.0)) +
			sqrt(3 * b_const*tet.temperature*bulk_viscosity*tet.volume) * multiply_by_identity_matrix(tet.gaussMatrix.getDiagonalSum() / 3.0)) / (tet.volume * sqrt(time_step / 2));

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

double calc_force2(int index)
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

//calculate average system velocity and subtract it from every particle velocity (necessary because of the periodic conditions)
void subtract_avg_velocity()
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

void calc_new_velocity_l()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->velocity = particles[i]->velocity + calc_force1(i) * time_step / 2;
	}
	subtract_avg_velocity();
}

void calc_new_velocity_o()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->velocity = particles[i]->velocity * exp(calc_force2(i) * time_step);
	}
	subtract_avg_velocity();
}

//LOD integration scheme (BAOAB scheme)
//https://doi.org/10.1063/1.5030034
void lod()
{
	calc_new_velocity_l();
	calc_new_position();
	calc_new_velocity_o();
	calc_new_position();
	calc_new_velocity_l();
}