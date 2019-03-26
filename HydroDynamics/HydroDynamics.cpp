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
void swap_and_clear(vector<Point3>* points)
{
	points->clear();
	tets.clear();
	for (int i = 0; i < particles.size(); i++)
	{
		double x = new_particles[i]->coordinates.x();
		double y = new_particles[i]->coordinates.y();
		double z = new_particles[i]->coordinates.z();

		particles[i]->clear();
		particles[i]->density = new_particles[i]->density;
		particles[i]->velocity = new_particles[i]->velocity;

		if (x > box_size) x -= box_size;
		else if (x < 0.0) x += box_size;

		if (y > box_size) y -= box_size;
		else if (y < 0.0) y += box_size;

		if (z > box_size) z -= box_size;
		else if (z < 0.0) z += box_size;

		particles[i]->coordinates = Point3(x, y, z);
		points->push_back(particles[i]->coordinates);
	}
}
/*
* add additional points to points vector (+26 points in the the vicinity of each initial point); the number of iterations = initial number of points
* all added points save to points vector in addition to initial ones
* then build tetrahedrons for all these points and leave only those which have one of the particles (from new_particles) in its corner
*/
void pereodic(vector<Point3>* points, Fade_3D& dt)
{
	//for every point in points we add neighbor points (26 for each) and push them to the initial points vector
	for (int i = 0, n = points->size(); i < n; i++)
	{
		auto x = (*points)[i].x();
		auto y = (*points)[i].y();
		auto z = (*points)[i].z();

		points->push_back(Point3(x - box_size, y - box_size, z - box_size));
		points->push_back(Point3(x - box_size, y - box_size, z));
		points->push_back(Point3(x - box_size, y - box_size, z + box_size));
		points->push_back(Point3(x - box_size, y, z - box_size));
		points->push_back(Point3(x - box_size, y, z));
		points->push_back(Point3(x - box_size, y, z + box_size));
		points->push_back(Point3(x - box_size, y + box_size, z - box_size));
		points->push_back(Point3(x - box_size, y + box_size, z));
		points->push_back(Point3(x - box_size, y + box_size, z + box_size));
		points->push_back(Point3(x, y - box_size, z - box_size));
		points->push_back(Point3(x, y - box_size, z));
		points->push_back(Point3(x, y - box_size, z + box_size));
		points->push_back(Point3(x, y, z - box_size));
		points->push_back(Point3(x, y, z + box_size));
		points->push_back(Point3(x, y + box_size, z - box_size));
		points->push_back(Point3(x, y + box_size, z));
		points->push_back(Point3(x, y + box_size, z + box_size));
		points->push_back(Point3(x + box_size, y - box_size, z - box_size));
		points->push_back(Point3(x + box_size, y - box_size, z));
		points->push_back(Point3(x + box_size, y - box_size, z + box_size));
		points->push_back(Point3(x + box_size, y, z - box_size));
		points->push_back(Point3(x + box_size, y, z));
		points->push_back(Point3(x + box_size, y, z + box_size));
		points->push_back(Point3(x + box_size, y + box_size, z - box_size));
		points->push_back(Point3(x + box_size, y + box_size, z));
		points->push_back(Point3(x + box_size, y + box_size, z + box_size));
	}

	//get vector of tetrahedrons from our vector of points
	dt.insert(*points);
	dt.getTetrahedra(tets);

	//loop for every tetrahedron to check if there is a particle (from new_particles) which lay in one of its corners. If there is not we erase this tetrahedron
	for (int i = 0; i < tets.size(); i++)
	{
		bool a = false;
		//loop for every corner of i-th tetrahedron
		for (int j = 0; j < 4; j++)
		{
			//check if there is any particle which lay in the j corner of the i-th tetrahedron
			for (int k = 0; k < new_particles.size(); k++)
			{
				if (new_particles[k]->coordinates.x() == tets[i]->getCorner(j)->x() && new_particles[k]->coordinates.y() == tets[i]->getCorner(j)->y() && new_particles[k]->coordinates.z() == tets[i]->getCorner(j)->z())
				{
					a = true;
					break;
				}
			}
			if (a)
			{
				break;
			}
		}
		if (!a)
		{
			tets.erase(tets.begin() + i);
		}
	}
}

void display_results(vector<Point3>* points)
{
	for (int i = 0; i < particles.size(); i++)
	{
		cout.precision(4);
		cout << fixed;
		cout << "coords : ( " << (*points)[i].x() << ", " << (*points)[i].y() << ", " << (*points)[i].z() << " ) \t";
		cout.precision(4);
		cout << scientific;
		cout << " velocity : ( " << particles[i]->velocity.x() << ", " << particles[i]->velocity.y() << ", " <<
			particles[i]->velocity.z() << " )  \t"
			<< " mass :" << particles[i]->density * particles[i]->volume << endl;
		cout << "---------" << endl;
	}
}

void display_coords(vector<Point3>* points)
{
	for (int i = 0; i < particles.size(); i++)
	{
		cout.precision(4);
		cout << fixed;
		cout << "x :" << (*points)[i].x() << " " << (*points)[i].y() << " " << (*points)[i].z() << endl;
	}
}

void mainCount(vector<Point3>* points)
{
	//variable to interact with Fade_3D methods
	Fade_3D dt;

	pereodic(points, dt);
	for (int i = 0; i < particles.size(); i++)
	{
		tetsCount(i);
		calcVolume(i);
		calcNewDensity(i);
		calcNewVelocity(i);
		calcNewPosition(i);
	}

	swap_and_clear(points);
	display_results(points);
}

int main(int argc, char * argv[])
{
	vector<Point3> points;

	//redirect output to Results.txt
	std::ofstream out("Results.txt");
	std::streambuf* coutbuf = std::cout.rdbuf();
	std::cout.rdbuf(out.rdbuf());

	//Create random particles, save their values to new_particles and save their coordinates to points vector
	for (int i = 0; i < particles.size(); i++)
	{
		double x, y, z;

		x = rand() % box_size + 1;
		y = rand() % box_size + 1;
		z = rand() % box_size + 1;

		cout << "x :" << x << " y :" << y << " z :" << z << endl;
		particles[i] = new HParticle(x, y, z);
		new_particles[i] = new HParticle(x, y, z);
		points.push_back(particles[i]->coordinates);
	}

	//for each particle start mainCount function
	for (int k = 0; k < 1; k++)
	{
		cout << "__________________ Series " << k << "__________________" << endl;
		mainCount(&points);
	}

	for (int i = 0; i < particles.size(); i++)
	{
		delete particles[i];
		delete new_particles[i];
	}
}