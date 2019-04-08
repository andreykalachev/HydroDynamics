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

int iteration = 0;
/*
* clear all particle fields for every particle; take density's and velocity's value from new_particles
* vor every new_particles->coordinates take their value [% box_size] and assign to particles->coordinates
*/
void swap_and_clear()
{
	tets.clear();
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->clear();
		particles[i]->density = new_particles[i]->density;
		particles[i]->velocity = new_particles[i]->velocity;
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

void display_results()
{
	double den_sum = 0, mass_sum = 0, temp_sum = 0, vel_sum = 0;
	for (int i = 0; i < number_of_particles; i++)
	{
		mass_sum += particles[i]->mass;
		den_sum += particles[i]->density;
		vel_sum += pow(particles[i]->velocity.x(), 2) + pow(particles[i]->velocity.y(), 2) + pow(particles[i]->velocity.z(), 2);
		temp_sum += particles[i]->temperature;
		//cout.precision(4);
		//cout << fixed;
		cout << "coords : ( " << particles[i]->coordinates.x() << ", " << particles[i]->coordinates.y() << ", " << particles[i]->coordinates.z() << " ) \t";
		//cout << scientific;
		cout << " velocity : ( " << particles[i]->velocity.x() << ", " << particles[i]->velocity.y() << ", " << particles[i]->velocity.z() << " )  \t"
			<< " mass :" << particles[i]->mass << "   \t" << " temp:" << particles[i]->temperature << "    \t" << "density:" << particles[i]->density << endl;
		cout << "---------" << endl;
	}
	cout << "density: " << den_sum / number_of_particles << "\t" << "mass: " << mass_sum / number_of_particles << "\t" << "temp: " << temp_sum / number_of_particles << "\t" << "vel: " << vel_sum / number_of_particles << endl;
	cout << "---------" << endl;
}

void display_coords()
{
	for (int i = 0; i < number_of_particles; i++)
	{
		//cout.precision(4);
		//cout << fixed;
		cout << "x :" << particles[i]->coordinates.x() << " " << particles[i]->coordinates.y() << " " << particles[i]->coordinates.z() << endl;
	}
}

void simulate_iteration()
{
	Fade_3D fade_3d;
	triangulation(fade_3d);
	for (int i = 0; i < number_of_particles; i++)
	{
		tetsCount(i);
		calcVolume(i);
		calcNewVelocity(i);
		calcNewDensity(i);
		calcTempreture(i);
		calcNewPosition(i);
	}
	swap_and_clear();
	display_results();
}


int main(int argc, char * argv[])
{
	//redirect output to Results.txt
	std::ofstream out("Results.txt");
	std::streambuf* coutbuf = std::cout.rdbuf();
	std::cout.rdbuf(out.rdbuf());

	//Create random particles, save their values to new_particles and save their coordinates to points vector
	for (int i = 0; i < number_of_particles; i++)
	{
		auto x = fRand(box_size * 0.1, box_size * 0.9);
		auto y = fRand(box_size * 0.1, box_size * 0.9);
		auto z = fRand(box_size * 0.1, box_size * 0.9);

		cout << "x :" << x << " y :" << y << " z :" << z << endl;
		particles[i] = new HParticle(x, y, z);
		particles[i]->velocity = Point3(fRand(-5, 5), fRand(-5, 5), fRand(-5, 5));
		new_particles[i] = new HParticle();
	}

	//for each particle start mainCount function
	for (int k = 0; k < 1000; k++)
	{
		iteration = k;
		cout << "__________________ Series " << k << "__________________" << endl;
		simulate_iteration();
	}
}