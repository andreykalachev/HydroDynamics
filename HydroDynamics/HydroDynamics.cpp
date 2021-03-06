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
	system_velocity->init(Point3(0,0,0));
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
	for (int i = 0; i < number_of_particles; i++)
	{
		calcNewVelocity(i);
		calcNewDensity(i);
	}
	calcNewPosition();
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

void display_results(ofstream &file)
{
	file << "__________________ Series " << iteration << "__________________" << endl;
	double den_avg = 0, mass_avg = 0, temp_avg = 0, vel_avg = 0, mom_avg = 0;
	for (int i = 0; i < number_of_particles; i++)
	{
		mass_avg += particles[i]->mass / number_of_particles;
		den_avg += particles[i]->density / number_of_particles;
		temp_avg += particles[i]->temperature / number_of_particles;
		particles[i]->display(file);
	}
	file << endl << "density: " << den_avg << "\t mass: " << mass_avg << "\t temp: " << temp_avg << "\t vel: " << vel_avg << "\t mom: " << mom_avg
		<< "\t syst_vel: (" << system_velocity->x() << ", " << system_velocity->y() << ", " << system_velocity->z() << ")" << endl;
}

void display_for_plot(ofstream &file)
{
	file << fixed << elapsed_time << " \t";
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->display_for_plot(file);
	}
	file << scientific << system_velocity->x() << " \t" << system_velocity->y() << " \t" << system_velocity->z() << " \t" << endl;
}

void simulate_iteration(ofstream &file)
{
	const auto fade_3d = new Fade_3D();
	triangulation(fade_3d);
	calculate_next_step();
	swap_and_clear();
	display_for_plot(file);
	delete fade_3d;
}


int main(int argc, char * argv[])
{
	/*set output to the file*/
	ofstream file;
	file.precision(5);
	file.open("Results.txt");

	//coords and velocity initialization 
	for (int i = 0; i < number_of_particles; i++)
	{ 
		new_particles[i] = new HParticle();
		particles[i] = createRandomParticle(0 , box_size);
		particles[i]->velocity = createRandomPoint(-5e-1, 5e-1);
		*system_velocity += particles[i]->velocity;
	}

	//subtract avg system velocity from particles velocities
	*system_velocity /= number_of_particles;
	for (int i = 0; i < number_of_particles; i++) 
	{
		particles[i]->velocity -= *system_velocity;
	}

	//for each particle start mainCount function
	for (int k = 0; k < 10000; k++)
	{
		iteration = k;
		simulate_iteration(file);
	}

	file.close();
}
