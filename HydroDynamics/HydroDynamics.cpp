#include "stdafx.h"
#include "fade3d/include_fade3d/Fade_3D.h"
#include "HParticle.h"
#include "CalculatingModule.h"
#include <iostream>
#include <ctime>
#include <iomanip>
#include <chrono>


using namespace std;
using namespace FADE3D;


void calculate_next_step()
{
	calcTempreture();
	if (iteration > 0) particles_image.erase(particles_image.begin(), particles_image.end());
	for (int i = 0, n = number_of_particles; i < n; i++)
	{
		auto x = particles[i]->coordinates.x();
		auto y = particles[i]->coordinates.y();
		auto z = particles[i]->coordinates.z();

		for (int p = -1; p <= 1; p++)
		{
			for (int j = -1; j <= 1; j++)
			{
				for (int k = -1; k <= 1; k++)
				{
					particles_image.push_back(particles[i]->get_copy(x + box_size * p, y + box_size * j, z + box_size * k));
				}
			}
		}
	}
	analyzeTets();
	calcNewDensity();
	lod();
}

/*
* clear all particle fields for every particle; take density's and velocity's value from new_particles
* vor every new_particles->coordinates take their value [% box_size] and assign to particles->coordinates
*/
void swap_and_clear()
{
	delaunay_tets.clear();
	elapsed_time += time_step;
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->clear();
		particles[i]->copy(new_particles[i]);
	}
}
/*
* add additional points to points vector (+26 points in the the vicinity of each initial point); the number of iterations = initial number of points
* all added points save to points vector in addition to initial ones
* then build tetrahedrons for all these points and leave only those which have one of the particles (from new_particles) in its corner
*/
void triangulation(Fade_3D *fade_3d)
{
	vector<Point3> points;
	//for every point in points we add neighbor points (26 for each) and push them to the initial points vector
	for (int i = 0, n = number_of_particles; i < n; i++)
	{
		auto x = particles[i]->coordinates.x();
		auto y = particles[i]->coordinates.y();
		auto z = particles[i]->coordinates.z();

		//create images of particles for periodic conditions
		for (int p = -1; p <= 1; p++)
		{
			for (int j = -1; j <= 1; j++)
			{
				for (int k = -1; k <= 1; k++)
				{
					points.push_back(Point3(x + box_size * p, y + box_size * j, z + box_size * k));
				}
			}
		}
	}

	//get vector of tetrahedrons from our vector of points
	fade_3d->insert(points);
	fade_3d->getTetrahedra(delaunay_tets);

	//loop for every tetrahedron to check if there is a particle (from new_particles) which lay in one of its corners. If there is not we erase this tetrahedron
	for (int i = 0; i < delaunay_tets.size(); i++)
	{
		bool b = false;
		for (int j = 0; j < number_of_particles; j++)
		{
			if (delaunay_tets[i]->hasVertex(particles[j]->coordinates))
			{
				b = true; break;
			}
		}
		if (!b) delaunay_tets.erase(delaunay_tets.begin() + i);
	}
	points.clear();
}

void display_for_plot(ofstream &velocity_file, ofstream &density_file, ofstream &velocity_fluctuations_file, ofstream &density_fluctuations_file, ofstream &coords_file, ofstream &volume_file)
{
	velocity_file << fixed << elapsed_time << " \t";
	density_file << fixed << elapsed_time << " \t";
	velocity_fluctuations_file << fixed << elapsed_time << " \t";
	density_fluctuations_file << fixed << elapsed_time << " \t";
	coords_file << fixed << elapsed_time << " \t";
	double avg_particlesvolume = 0;
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->display_velocity(velocity_file);
		particles[i]->display_density(density_file);
		particles[i]->display_velocity_fluctuations(velocity_fluctuations_file, system_velocity);
		particles[i]->display_density_fluctuations(density_fluctuations_file, system_density);
		particles[i]->display_coords(coords_file);
		avg_particlesvolume += particles[i]->volume;
		volume_file << fixed << particles[i]->volume << " \t";
	}
	velocity_file << fixed << system_velocity->x() * 1e3 << " \t" << system_velocity->y() * 1e3 << " \t" << system_velocity->z() * 1e3 << " \t" << endl;
	density_file << fixed << system_density * 1.66e3 << " \t" << endl;
	velocity_fluctuations_file << fixed << system_velocity->x() * 1e3 << " \t" << system_velocity->y() * 1e3 << " \t" << system_velocity->z() * 1e3 << " \t" << endl;
	density_fluctuations_file << fixed << system_density * 1.66e3 << " \t" << endl;
	coords_file << endl;
	volume_file << avg_particlesvolume / number_of_particles << "\t" << endl;
}

int main(int argc, char * argv[])
{
	/*set output to the file*/
	ofstream velocity_file, density_file, velocity_fluctuations_file, density_fluctuations_file, coords_file, volume_file;
	velocity_file.precision(6);
	density_file.precision(6);
	velocity_fluctuations_file.precision(6);
	density_fluctuations_file.precision(6);
	coords_file.precision(6);
	volume_file.precision(6);
	velocity_file.open("velocity.txt");
	density_file.open("density.txt");
	velocity_fluctuations_file.open("velocity_fluctuations.txt");
	density_fluctuations_file.open("density_fluctuations.txt");
	coords_file.open("coordinates.txt");
	volume_file.open("volume.txt");

	//coords and velocity initialization 
	for (int i = 0; i < number_of_particles; i++)
	{
		new_particles[i] = new HParticle();
		particles[i] = createRandomParticle(0, box_size);
		//particles[i]->velocity = createRandomPoint(-5e-3, 5e-3);
		//particles[i]->velocity = *new Point3(velocity_distribution(generator), velocity_distribution(generator), velocity_distribution(generator));
		particles[i]->velocity = *new Point3(0, 0, 0);
	}

	//subtract avg system velocity from particles velocities
	substractAvgVel();

	//for each particle start mainCount function
	for (iteration = 0; iteration < 100000; iteration++)
	{
		const auto fade_3d = new Fade_3D();
		triangulation(fade_3d);
		calculate_next_step();
		swap_and_clear();
		display_for_plot(velocity_file, density_file, velocity_fluctuations_file, density_fluctuations_file, coords_file, volume_file);
		delete fade_3d;
	}

	velocity_file.close();
	density_file.close();
	velocity_fluctuations_file.close();
	density_fluctuations_file.close();
	coords_file.close();
	volume_file.close();
}