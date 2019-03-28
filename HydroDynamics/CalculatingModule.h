#include <math.h>

//number of particles, new_particles
const int n = 10;
const int box_size = 100;
const double h = 0.1;
const double bolzmana = 1;
std::vector<HParticle*> particles(n);
vector<HParticle*> new_particles(n);

static bool compare_coords(vector<Tet3*>::value_type& tet, int tet_corner, HParticle* particle)
{
	return tet->getCorner(tet_corner)->x() == particle->coordinates.x() && tet->getCorner(tet_corner)->y() == particle->coordinates.y() &&
		tet->getCorner(tet_corner)->z() == particle->coordinates.z();
}

static int positive_mod(int x, int y) {
	return (x % y + y) % y;
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
	for (auto& tet : tets)
	{
		double tet_density = 0, tet_temperature = 0, tet_velocity_x = 0, tet_velocity_y = 0, tet_velocity_z = 0, tet_volume = 0;

		for (int corner = 0; corner < 4; corner++)
		{
			//if index particle lie in the nth corner of the tetrahedron
			if (compare_coords(tet, corner, particles[index]))
			{
				//check all other corners of the tetrahedrons and add if there are any particles in these corners add them to the neighbors
				for (int j = 0; j < particles.size(); j++)
				{
					if (compare_coords(tet, (corner + 1) % 4, particles[j]) || compare_coords(tet, (corner + 2) % 4, particles[j]) ||
						compare_coords(tet, (corner + 3) % 4, particles[j]))
					{
						particles[index]->neighbours_points.push_back(particles[j]);
					}
					//loop for each corner of tet
					for (int corner_index = 0; corner_index < 4; corner_index++)
					{
						if (compare_coords(tet, corner_index, particles[j]))
						{
							tet_density += particles[j]->density / 4;
							tet_velocity_x += particles[j]->velocity.x() / 4;
							tet_velocity_y += particles[j]->velocity.y() / 4;
							tet_velocity_z += particles[j]->velocity.z() / 4;
							tet_temperature += particles[j]->temperature / 4;
						}
					}
				}
				tet_volume = calcTetVolume(tet);
				particles[index]->tetsVolumes.push_back(tet_volume);
				particles[index]->tetsVelocities.push_back(Point3(tet_velocity_x, tet_velocity_y, tet_velocity_z));
				particles[index]->tetsTemperature.push_back(tet_temperature);
				particles[index]->tetsDensities.push_back(tet_density);
				particles[index]->vectorsB.push_back(calculateVectorB(tet, corner, tet_volume));
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

void calcTempreture(int index)
{
	auto particle = particles[index];
	double velocity_squared = pow(particle->velocity.x(), 2) + pow(particle->velocity.y(), 2) + pow(particle->velocity.z(), 2);
	particle->temperature = particle->density * particle->volume * velocity_squared / bolzmana / 3;

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
	particle->mass = particle->volume * particle->density;

	for (int i = 0; i < particle->tetsDensities.size(); i++)
	{
		double pressure = calcPressure(particle->tetsDensities[i]);
		Point3 tets_velocity = particle->tetsVelocities[i];
		Point3 vector_b = particle->vectorsB[i];
		double tets_density = particle->tetsDensities[i];

		term1_x += particle->tetsVolumes[i] * (vector_b.x()*(pressure + pow(tets_velocity.x(), 2)*tets_density)
			+ vector_b.y()*tets_velocity.y()*tets_velocity.x()*tets_density + vector_b.z()*tets_velocity.z()*tets_velocity.x()*tets_density);

		term1_y += particle->tetsVolumes[i] * (vector_b.x()*tets_velocity.x()*tets_velocity.y()*tets_density
			+ vector_b.y()*(pressure + pow(tets_velocity.y(), 2)*tets_density) + vector_b.z()*tets_velocity.y()*tets_velocity.z()*tets_density);

		term1_z += particle->tetsVolumes[i] * (vector_b.x()*tets_velocity.z()*tets_velocity.x()*tets_density
			+ vector_b.y()*tets_velocity.z()*tets_velocity.y()*tets_density + vector_b.z()*(pressure + pow(tets_velocity.z(), 2)*tets_density));
	}

	for (int i = 0; i < particle->neighbours_points.size(); i++)
	{
		HParticle* neighbor = particle->neighbours_points[i];

		for (int j = 0; j < neighbor->tetsDensities.size(); j++)
		{
			for (int k = 0; k < particle->tetsDensities.size(); k++)
			{
				//if particle and its neighbor have common tetrahedron
				if (neighbor->tetsDensities[j] == particle->tetsDensities[k] &&
					neighbor->tetsVelocities[j].x() == particle->tetsVelocities[k].x() &&
					neighbor->tetsVelocities[j].y() == particle->tetsVelocities[k].y() &&
					neighbor->tetsVelocities[j].z() == particle->tetsVelocities[k].z() &&
					neighbor->tetsVolumes[j] == particle->tetsVolumes[k])
				{
					double scalar_product1 = (particle->tetsVelocities[k].x() - neighbor->velocity.x()) * neighbor->vectorsB[j].x() +
						(particle->tetsVelocities[k].y() - neighbor->velocity.y()) * neighbor->vectorsB[j].y() +
						(particle->tetsVelocities[k].z() - neighbor->velocity.z()) * neighbor->vectorsB[j].z();

					double scalar_product2 = particle->vectorsB[k].x() * neighbor->vectorsB[j].x() +
						particle->vectorsB[k].y() * neighbor->vectorsB[j].y() +
						particle->vectorsB[k].z() * neighbor->vectorsB[j].z();

					double scalar_product3 = (particle->tetsVelocities[k].x() - neighbor->velocity.x()) * particle->vectorsB[k].x() +
						(particle->tetsVelocities[k].y() - neighbor->velocity.y()) * particle->vectorsB[k].y() +
						(particle->tetsVelocities[k].z() - neighbor->velocity.z()) * particle->vectorsB[k].z();

					term2_x += particle->tetsTemperature[k] * particle->tetsVolumes[k] * (scalar_product1 * particle->vectorsB[k].x() +
						(particle->tetsVelocities[k].x() - neighbor->velocity.x()) * scalar_product2 +
						scalar_product3 * neighbor->vectorsB[j].x() - 2.0 / 3.0 * scalar_product1 * particle->vectorsB[k].x()) / neighbor->temperature;

					term2_y += particle->tetsTemperature[k] * particle->tetsVolumes[k] * (scalar_product1 * particle->vectorsB[k].y() +
						(particle->tetsVelocities[k].y() - neighbor->velocity.y()) * scalar_product2 +
						scalar_product3 * neighbor->vectorsB[j].y() - 2.0 / 3.0 * scalar_product1 * particle->vectorsB[k].y()) / neighbor->temperature;

					term2_z += particle->tetsTemperature[k] *  particle->tetsVolumes[k] * (scalar_product1 * particle->vectorsB[k].z() +
						(particle->tetsVelocities[k].z() - neighbor->velocity.z()) * scalar_product2 +
						scalar_product3 * neighbor->vectorsB[j].z() - 2.0 / 3.0 * scalar_product1 * particle->vectorsB[k].z()) / neighbor->temperature;
				}
			}
		}
	}

	double result_x = (term1_x + term2_x);
	double result_y = (term1_y + term2_y);
	double result_z = (term1_z + term2_z);

	return Point3(result_x, result_y, result_z);
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
	new_particles[index]->coordinates = Point3(particles[index]->coordinates.x() + h * new_particles[index]->velocity.x(),
		particles[index]->coordinates.y() + h * new_particles[index]->velocity.y(),
		particles[index]->coordinates.z() + h * new_particles[index]->velocity.z());
}