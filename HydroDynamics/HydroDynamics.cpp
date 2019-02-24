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
	static const double eps = 0.01;
	for (int i = 0; i < particles.size(); i++)
	{
		double x = new_particles[i]->coordinates.x();
		double y = new_particles[i]->coordinates.y();
		double z = new_particles[i]->coordinates.z();


		particles[i]->vectorsB.clear();
		particles[i]->neighbours_points.clear();
		particles[i]->tetsVelocities.clear();
		particles[i]->tetsDensities.clear();
		particles[i]->tetsVolumes.clear();
		particles[i]->density = new_particles[i]->density;
		particles[i]->velocity = new_particles[i]->velocity;

		if (x - box_size > 0.0)
			x -= box_size;
		if (x < 0.0)
			x += box_size;
		if (y - box_size > 0.0)
			y -= box_size;
		if (y < 0.0)
			y += box_size;
		if (z - box_size > 0.0)
			z -= box_size;
		if (z < 0.0)
			z += box_size;

		particles[i]->coordinates = Point3(x, y, z);
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
		Point3 point1 = Point3((*points)[i].x() - 50, (*points)[i].y() - 50, (*points)[i].z() - 50);
		Point3 point2 = Point3((*points)[i].x() - 50, (*points)[i].y() - 50, (*points)[i].z());
		Point3 point3 = Point3((*points)[i].x() - 50, (*points)[i].y() - 50, (*points)[i].z() + 50);
		Point3 point4 = Point3((*points)[i].x() - 50, (*points)[i].y(), (*points)[i].z() - 50);
		Point3 point5 = Point3((*points)[i].x() - 50, (*points)[i].y(), (*points)[i].z());
		Point3 point6 = Point3((*points)[i].x() - 50, (*points)[i].y(), (*points)[i].z() + 50);
		Point3 point7 = Point3((*points)[i].x() - 50, (*points)[i].y() + 50, (*points)[i].z() - 50);
		Point3 point8 = Point3((*points)[i].x() - 50, (*points)[i].y() + 50, (*points)[i].z());
		Point3 point9 = Point3((*points)[i].x() - 50, (*points)[i].y() + 50, (*points)[i].z() + 50);
		Point3 point10 = Point3((*points)[i].x(), (*points)[i].y() - 50, (*points)[i].z() - 50);
		Point3 point11 = Point3((*points)[i].x(), (*points)[i].y() - 50, (*points)[i].z());
		Point3 point12 = Point3((*points)[i].x(), (*points)[i].y() - 50, (*points)[i].z() + 50);
		Point3 point13 = Point3((*points)[i].x(), (*points)[i].y(), (*points)[i].z() - 50);
		//Point3 point14 = Point3(points[i].x(), points[i].y(), points[i].z());
		Point3 point15 = Point3((*points)[i].x(), (*points)[i].y(), (*points)[i].z() + 50);
		Point3 point16 = Point3((*points)[i].x(), (*points)[i].y() + 50, (*points)[i].z() - 50);
		Point3 point17 = Point3((*points)[i].x(), (*points)[i].y() + 50, (*points)[i].z());
		Point3 point18 = Point3((*points)[i].x(), (*points)[i].y() + 50, (*points)[i].z() + 50);
		Point3 point19 = Point3((*points)[i].x() + 50, (*points)[i].y() - 50, (*points)[i].z() - 50);
		Point3 point20 = Point3((*points)[i].x() + 50, (*points)[i].y() - 50, (*points)[i].z());
		Point3 point21 = Point3((*points)[i].x() + 50, (*points)[i].y() - 50, (*points)[i].z() + 50);
		Point3 point22 = Point3((*points)[i].x() + 50, (*points)[i].y(), (*points)[i].z() - 50);
		Point3 point23 = Point3((*points)[i].x() + 50, (*points)[i].y(), (*points)[i].z());
		Point3 point24 = Point3((*points)[i].x() + 50, (*points)[i].y(), (*points)[i].z() + 50);
		Point3 point25 = Point3((*points)[i].x() + 50, (*points)[i].y() + 50, (*points)[i].z() - 50);
		Point3 point26 = Point3((*points)[i].x() + 50, (*points)[i].y() + 50, (*points)[i].z());
		Point3 point27 = Point3((*points)[i].x() + 50, (*points)[i].y() + 50, (*points)[i].z() + 50);

		points->push_back(point1);
		points->push_back(point2);
		points->push_back(point3);
		points->push_back(point4);
		points->push_back(point5);
		points->push_back(point6);
		points->push_back(point7);
		points->push_back(point8);
		points->push_back(point9);
		points->push_back(point10);
		points->push_back(point11);
		points->push_back(point12);
		points->push_back(point13);
		//points.push_back(point14);
		points->push_back(point15);
		points->push_back(point16);
		points->push_back(point17);
		points->push_back(point18);
		points->push_back(point19);
		points->push_back(point20);
		points->push_back(point21);
		points->push_back(point22);
		points->push_back(point23);
		points->push_back(point24);
		points->push_back(point25);
		points->push_back(point26);
		points->push_back(point27);
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
	swap_and_clear();
	points->clear();


	for (int i = 0; i < particles.size(); i++)
	{
		points->push_back(particles[i]->coordinates);
	}

	tets.clear();

	pereodic(points, dt);
	for (int i = 0; i < particles.size(); i++)
	{
		tetsCount(i);
		calcVolume(i);
		calcNewDensity(i);
		calcNewVelocity(i);
	}
	swap_and_clear();
	points->clear();


	//console log of particles
	for (int i = 0; i < particles.size(); i++)
	{
		points->push_back(particles[i]->coordinates);
		cout.precision(4);
		cout << fixed;
		cout << "coordinates : ( " << (*points)[points->size() - 1].x() << ", " << (*points)[points->size() - 1].y() << ", " << (*points)[points->size() - 1].z() << " )\t";
		cout.precision(4);
		cout << scientific;
		cout << " velocity : ( " << particles[i]->velocity.x() << ", " << particles[i]->velocity.y() << ", " << particles[i]->velocity.z() << " )\t"
			<< " mass :" << particles[i]->density *particles[i]->volume << endl;
		cout << "---------" << endl;
	}

	tets.clear();

}

int main(int argc, char * argv[])
{
	vector<Point3> points;

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

	//OpenGL drawing
	//initGL(argc, argv);

	//for each particle start mainCount function
	for (int k = 0; k < 10; k++)
	{
		cout << "__________________ Series " << k << "__________________" << endl;
		mainCount(&points);
	}

	for (int i = 0; i < particles.size(); i++)
	{
		delete particles[i];
		delete new_particles[i];
	}
	return 0;

}