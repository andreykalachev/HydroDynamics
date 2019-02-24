#include <math.h>

//(initial) number of particles, new_particles
const int n = 10;
const int box_size = 100;
const double h = 0.1;
std::vector<HParticle*> particles(n);
vector<HParticle*> new_particles(n);

bool compareCoordinates(vector<Tet3*>::value_type& tet, int tet_corner, HParticle* particle)
{
	return tet->getCorner(tet_corner)->x() == particle->coordinates.x() && tet->getCorner(tet_corner)->y() == particle->coordinates.y() &&
		tet->getCorner(tet_corner)->z() == particle->coordinates.z();
}

double calcTetVolume(vector<Tet3*>::value_type& tet)
{
	auto tet_volume = (tet->getCorner(0)->x() - tet->getCorner(1)->x()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z())
		+ (tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->z() - tet->getCorner(1)->z()) +
		(tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(0)->z() - tet->getCorner(1)->z()) * (tet->getCorner(2)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->z() - tet->getCorner(1)->z()) * (tet->getCorner(3)->y() - tet->getCorner(1)->y()) * (tet->getCorner(0)->x() - tet->getCorner(1)->x()) -
		(tet->getCorner(2)->x() - tet->getCorner(1)->x()) * (tet->getCorner(0)->y() - tet->getCorner(1)->y()) * (tet->getCorner(3)->z() - tet->getCorner(1)->z());

	return abs(tet_volume);
}

Point3 calculateVectorB(vector<Tet3*>::value_type& tet, int corner, double tet_volume)
{
	// i don't know why
	auto sign = corner % 2 == 0 ? 1 : -1;

	auto _1 = abs(corner + sign * 1) % 4;
	auto _2 = abs(corner + sign * 2) % 4;
	auto _3 = abs(corner + sign * 3) % 4;

	auto x = -(tet->getCorner(_2)->y() - tet->getCorner(_1)->y()) *
		(tet->getCorner(_3)->z() - tet->getCorner(_1)->z()) +
		(tet->getCorner(_3)->y() - tet->getCorner(_1)->y()) *
		(tet->getCorner(_2)->z() - tet->getCorner(_1)->z());

	auto y = -(tet->getCorner(_2)->z() - tet->getCorner(_1)->z()) *
		(tet->getCorner(_3)->x() - tet->getCorner(_1)->x()) +
		(tet->getCorner(_3)->z() - tet->getCorner(_1)->z()) *
		(tet->getCorner(_2)->x() - tet->getCorner(_1)->x());

	auto z = -(tet->getCorner(_2)->x() - tet->getCorner(_1)->x()) *
		(tet->getCorner(_3)->y() - tet->getCorner(_1)->y()) +
		(tet->getCorner(_3)->x() - tet->getCorner(_1)->x()) *
		(tet->getCorner(_2)->y() - tet->getCorner(_1)->y());

	return Point3(x / tet_volume, y / tet_volume, z / tet_volume);
}

//calculate tetrahedrons' neighbors, volume, density, velocity and vector B
void tetsCount(int index)
{
	//loop for each tetrahedron
	for (auto& tet : tets)
	{
		double tet_density = 0;
		double tet_velocity_x = 0, tet_velocity_y = 0, tet_velocity_z = 0;
		double tet_volume = 0;

		for (int corner = 0; corner < 4; corner++)
		{
			//if index particle lie in the nth corner of the tetrahedron
			if (compareCoordinates(tet, corner, particles[index]))
			{
				//check all other corners of the tetrahedrons and add if there are any particles in these corners add them to the neighbors
				for (int j = 0; j < particles.size(); j++)
				{
					if (compareCoordinates(tet, (corner + 1) % 4, particles[j]) || compareCoordinates(tet, (corner + 2) % 4, particles[j]) ||
						compareCoordinates(tet, (corner + 3) % 4, particles[j]))
					{
						particles[index]->neighbours_points.push_back(particles[j]);
					}
				}

				//calculate volume of the tetrahedron and take it absolute value
				tet_volume = calcTetVolume(tet);
				particles[index]->tetsVolumes.push_back(tet_volume);

				//calculate density and velocity
				for (int k = 0; k < particles.size(); k++)
				{
					//loop for each corner of tet
					for (int corner_index = 0; corner_index < 4; corner_index++)
					{
						if (compareCoordinates(tet, corner_index, particles[k]))
						{
							tet_density += particles[k]->density / 4;
							tet_velocity_x += particles[k]->velocity.x() / 4;
							tet_velocity_y += particles[k]->velocity.y() / 4;
							tet_velocity_z += particles[k]->velocity.z() / 4;
						}
					}
				}
				particles[index]->tetsVelocities.push_back(Point3(tet_velocity_x, tet_velocity_y, tet_velocity_z));
				particles[index]->tetsDensities.push_back(tet_density);

				// calculate vector B
				auto vector_b = calculateVectorB(tet, corner, tet_volume);
				particles[index]->vectorsB.push_back(vector_b);
			}
		}
	}
}

void calcVolume(int index)
{
	double volume = 0;

	for (int i = 0; i < particles[index]->tetsVolumes.size(); i++)
	{
		volume += particles[index]->tetsVolumes[i];
	}
	particles[index]->volume = volume;
}

void calcNewDensity(int index)
{
	double sum = 0;
	for (int i = 0; i < particles[index]->tetsVolumes.size(); i++)
	{
		sum += particles[index]->tetsVolumes[i] * particles[index]->tetsDensities[i]
			* (particles[index]->vectorsB[i].x() * particles[index]->tetsVelocities[i].x()
				+ particles[index]->vectorsB[i].y() * particles[index]->tetsVelocities[i].y()
				+ particles[index]->vectorsB[i].z() * particles[index]->tetsVelocities[i].z());
	}
	new_particles[index]->density = particles[index]->density + h * sum / particles[index]->volume;
}

double calcPressure(int particle_index, int tet_index)
{
	static double a, b, c, d, e, f;
	a = -1.86789e-004;
	b = 1.9878e-001;
	c = -9.6160e001;
	d = 8.2666e004;
	e = -1.5356e006;
	f = 1.1584e-007;

	double x = particles[particle_index]->tetsDensities[tet_index];
	double x_2 = x * x;

	return f * x_2*x_2*x + a * x_2*x_2 + b * x_2*x + c * x_2 + d * x + e;
}

Point3 calcForce(int index)
{
	double term1_x = 0, term1_y = 0, term1_z = 0, term2_x = 0, term2_y = 0, term2_z = 0;

	particles[index]->mass = particles[index]->volume * particles[index]->density;

	for (int i = 0; i < particles[index]->tetsDensities.size(); i++)
	{
		auto pressure = calcPressure(index, i);
		auto tets_velocity = particles[index]->tetsVelocities[i];
		auto vector_b = particles[index]->vectorsB[i];
		auto tets_density = particles[index]->tetsDensities[i];

		term1_x += particles[index]->tetsVolumes[i] * (vector_b.x()*(pressure + pow(tets_velocity.x(), 2)*tets_density)
			+ vector_b.y()*tets_velocity.y()*tets_velocity.x()*tets_density + vector_b.z()*tets_velocity.z()*tets_velocity.x()*tets_density);

		term1_y += particles[index]->tetsVolumes[i] * (vector_b.x()*tets_velocity.x()*tets_velocity.y()*tets_density
			+ vector_b.y()*(pressure + pow(tets_velocity.y(), 2)*tets_density) + vector_b.z()*tets_velocity.y()*tets_velocity.z()*tets_density);

		term1_z += particles[index]->tetsVolumes[i] * (vector_b.x()*tets_velocity.z()*tets_velocity.x()*tets_density
			+ vector_b.y()*tets_velocity.z()*tets_velocity.y()*tets_density + vector_b.z()*(pressure + pow(tets_velocity.z(), 2)*tets_density));
	}

	for (int i = 0; i < particles[index]->neighbours_points.size(); i++)
	{
		auto neighbor = particles[index]->neighbours_points[i];

		for (int j = 0; j < neighbor->tetsDensities.size(); j++)
		{
			for (int k = 0; k < particles[index]->tetsDensities.size(); k++)
			{
				//[???Maybe???] if particle and its neighbor have common tetrahedron
				if (neighbor->tetsDensities[j] == particles[index]->tetsDensities[k] &&
					neighbor->tetsVelocities[j].x() == particles[index]->tetsVelocities[k].x() &&
					neighbor->tetsVelocities[j].y() == particles[index]->tetsVelocities[k].y() &&
					neighbor->tetsVelocities[j].z() == particles[index]->tetsVelocities[k].z() &&
					neighbor->tetsVolumes[j] == particles[index]->tetsVolumes[k])
				{
					auto scalar_product1 = (particles[index]->tetsVelocities[k].x() - neighbor->velocity.x()) * neighbor->vectorsB[j].x() +
						(particles[index]->tetsVelocities[k].y() - neighbor->velocity.y()) * neighbor->vectorsB[j].y() +
						(particles[index]->tetsVelocities[k].z() - neighbor->velocity.z()) * neighbor->vectorsB[j].z();

					auto scalar_product2 = particles[index]->vectorsB[k].x() * neighbor->vectorsB[j].x() +
						particles[index]->vectorsB[k].y() * neighbor->vectorsB[j].y() +
						particles[index]->vectorsB[k].z() * neighbor->vectorsB[j].z();

					auto scalar_product3 = (particles[index]->tetsVelocities[k].x() - neighbor->velocity.x()) * particles[index]->vectorsB[k].x() +
						(particles[index]->tetsVelocities[k].y() - neighbor->velocity.y()) * particles[index]->vectorsB[k].y() +
						(particles[index]->tetsVelocities[k].z() - neighbor->velocity.z()) * particles[index]->vectorsB[k].z();

					term2_x += particles[index]->tetsVolumes[k] * (scalar_product1 * particles[index]->vectorsB[k].x() +
						(particles[index]->tetsVelocities[k].x() - neighbor->velocity.x()) * scalar_product2 +
						scalar_product3 * neighbor->vectorsB[j].x() - 2.0 / 3.0 * scalar_product1 * particles[index]->vectorsB[k].x());

					term2_y += particles[index]->tetsVolumes[k] * (scalar_product1 * particles[index]->vectorsB[k].y() +
						(particles[index]->tetsVelocities[k].y() - neighbor->velocity.y()) * scalar_product2 +
						scalar_product3 * neighbor->vectorsB[j].y() - 2.0 / 3.0 * scalar_product1 * particles[index]->vectorsB[k].y());

					term2_z += particles[index]->tetsVolumes[k] * (scalar_product1 * particles[index]->vectorsB[k].z() +
						(particles[index]->tetsVelocities[k].z() - neighbor->velocity.z()) * scalar_product2 +
						scalar_product3 * neighbor->vectorsB[j].z() - 2.0 / 3.0 * scalar_product1 * particles[index]->vectorsB[k].z());
				}
			}
		}
	}

	double result_x = (term1_x + term2_x) / particles[index]->volume;
	double result_y = (term1_y + term2_y) / particles[index]->volume;
	double result_z = (term1_z + term2_z) / particles[index]->volume;

	return Point3(result_x, result_y, result_z);
}

void calcNewVelocity(int index)
{
	Point3 force = calcForce(index);

	//Why ???
	/*
	new_particles[index]->velocity = Point3(particles[index]->velocity.x() + h * force.x() / particles[index]->mass,
	particles[index]->velocity.y() + h * force.y() / particles[index]->mass,
	particles[index]->velocity.z() + h * force.z() / particles[index]->mass);
	*/

	new_particles[index]->velocity = Point3(particles[index]->velocity.x() + h * force.x() / 2 / particles[index]->mass,
		particles[index]->velocity.y() + h * force.y() / 2 / particles[index]->mass,
		particles[index]->velocity.z() + h * force.z() / 2 / particles[index]->mass);
}

void calcNewPosition(int index)
{
	new_particles[index]->coordinates = Point3(particles[index]->coordinates.x() + h * new_particles[index]->velocity.x(),
		particles[index]->coordinates.y() + h * new_particles[index]->velocity.y(),
		particles[index]->coordinates.z() + h * new_particles[index]->velocity.z());
}
