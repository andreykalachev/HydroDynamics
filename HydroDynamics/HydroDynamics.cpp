#include "stdafx.h"
#include "GL/glut.h"
#include "fade3d/include_fade3d/Fade_3D.h"
#include "HParticle.h"
#include "DrawModule.h"
#include "CalculatingModule.h"
#include <iostream>
#include <ctime>
#include <iomanip>


using namespace std;
using namespace FADE3D;

/*
* clear all particle fields for every particle; take density's and velocity's value from new_particles
* vor every new_particles->coordinates take their value [% box_size] and assign to particles->coordinates
*/
void swap_and_clear()
{
	tets.clear();
	system_velocity = divide(system_velocity, number_of_particles);
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->clear();
		particles[i]->density = new_particles[i]->density;
		particles[i]->velocity = subtract(new_particles[i]->velocity, system_velocity);
		//particles[i]->velocity = new_particles[i]->velocity;
		particles[i]->momentum = new_particles[i]->momentum;
		particles[i]->velocity_absolute = new_particles[i]->velocity_absolute;
		particles[i]->momentum_absolute = new_particles[i]->momentum_absolute;
		particles[i]->temperature = new_particles[i]->temperature;

		auto x = positive_mod(new_particles[i]->coordinates.x(), box_size);
		auto y = positive_mod(new_particles[i]->coordinates.y(), box_size);
		auto z = positive_mod(new_particles[i]->coordinates.z(), box_size);

		particles[i]->coordinates = Point3(x, y, z);
	}
}
/*
* add additional points to points vector (+26 points in the the vicinity of each initial point); the number of iterations = initial number of points
* all added points save to points vector in addition to initial ones
* then build tetrahedrons for all these points and leave only those which have one of the particles (from new_particles) in its corner
*/
void triangulation(Fade_3D &fade_3d)
{
	vector<Point3> points;
	//for every point in points we add neighbor points (26 for each) and push them to the initial points vector
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
					points.push_back(Point3(x + box_size * p, y + box_size * j, z + box_size * k));
					if (iteration > 0) particles_image.erase(particles_image.begin());
					particles_image.push_back(particles[i]->copy(x + box_size * p, y + box_size * j, z + box_size * k));
				}
			}
		}
	}

	//get vector of tetrahedrons from our vector of points
	fade_3d.insert(points);
	fade_3d.getTetrahedra(tets);

	//loop for every tetrahedron to check if there is a particle (from new_particles) which lay in one of its corners. If there is not we erase this tetrahedron
	for (int i = 0; i < tets.size(); i++)
	{
		bool b = false;
		for (int j = 0; j < number_of_particles; j++)
		{
			if (tets[i]->hasVertex(particles[j]->coordinates))
			{
				b = true; break;
			}
		}
		if (!b) tets.erase(tets.begin() + i);
	}
}

void display_results(ofstream &file)
{
	file << "__________________ Series " << iteration << "__________________" << endl;
	double den_avg = 0, mass_avg = 0, temp_avg = 0, vel_avg = 0, mom_avg = 0;
	for (int i = 0; i < number_of_particles; i++)
	{
		mass_avg += particles[i]->mass / number_of_particles;
		den_avg += particles[i]->density / number_of_particles;
		vel_avg += particles[i]->velocity_absolute / number_of_particles;
		mom_avg += particles[i]->momentum_absolute / number_of_particles;
		temp_avg += particles[i]->temperature / number_of_particles;
		particles[i]->display(file);
	}
	file << endl << "density: " << den_avg << "\t mass: " << mass_avg << "\t temp: " << temp_avg << "\t vel: " << vel_avg << "\t mom: " << mom_avg
		<< "\t syst_vel: (" << system_velocity.x() << ", " << system_velocity.y() << ", " << system_velocity.z() << ")" << endl;
}

void display_for_plot(ofstream &file)
{
	file << fixed << iteration * time_step << " \t";
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->display_for_plot(file);
	}
	file << scientific << system_velocity.x() << " \t" << system_velocity.y() << " \t" << system_velocity.z() << " \t" << endl;
}

void simulate_iteration(ofstream &file)
{
	Fade_3D fade_3d;
	triangulation(fade_3d);
	tetsCount();
	for (int i = 0; i < number_of_particles; i++)
	{
		calcNewVelocity(i);
		calcNewDensity(i);
		calcTempreture(i);
		calcNewPosition(i);
	}
	swap_and_clear();
	display_for_plot(file);
}


int main(int argc, char * argv[])
{
	ofstream file;
	file.precision(5);
	file.open("Results.txt");

	//Create random particles, save their values to new_particles and save their coordinates to points vector
	for (int i = 0; i < number_of_particles; i++)
	{
		auto x = fRand(0, box_size);
		auto y = fRand(0, box_size);
		auto z = fRand(0, box_size);

		//file << x << " " << y << " " << z << endl;
		particles[i] = new HParticle(x, y, z);
		particles[i]->velocity = Point3(fRand(-5e-3, 5e-3), fRand(-5e-3, 5e-3), fRand(-5e-3, 5e-3));
		particles[i]->velocity_absolute = calculate_absolute_value(particles[i]->velocity);
		new_particles[i] = new HParticle();
		system_velocity = calculate_sum(system_velocity, particles[i]->velocity);
	}

	system_velocity = divide(system_velocity, number_of_particles);
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->velocity = subtract(particles[i]->velocity, system_velocity);
	}

	//for each particle start mainCount function
	for (int k = 0; k < 1000; k++)
	{
		if (iteration == 100) time_step *= 10;
		system_velocity = Point3(0, 0, 0);
		iteration = k;
		simulate_iteration(file);
	}

	file.close();
}