#include "stdafx.h"
#include "fade3d/include_fade3d/Fade_3D.h"
#include "HParticle.h"
#include "CalculatingModule.h"
#include <ctime>
#include <iomanip>
#include <chrono>


using namespace std;
using namespace FADE3D;

//analyze particles on current step and calculate velocity and density on next step
void calculate_next_step()
{
	calc_tempreture();
	if (iteration > 0) particles_image.erase(particles_image.begin(), particles_image.end());
	//create images of particles for periodic conditions
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
					particles_image.push_back(particles[i]->get_copy(x + box_length * p, y + box_length * j, z + box_length * k));
				}
			}
		}
	}
	analyze_tets();
	calc_new_density();
	lod();
}

//update density, clear neighbors and tetrahedrons for each particle
void update_params()
{
	elapsed_time += time_step;
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->clear();
		particles[i]->density = new_density[i];
	}
}

//do triangulation for particles and its images and delete all tetrahedrons that don't belong to any particle
void triangulation(Fade_3D *fade_3d)
{
	vector<Point3> points;
	//clear results from previous triangulation
	delaunay_tets.clear();
	//create images of particles for periodic conditions
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
					points.push_back(Point3(x + box_length * p, y + box_length * j, z + box_length * k));
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

//results output
void display_results(ofstream &velocity_file, ofstream &density_file, ofstream &velocity_fluctuations_file, ofstream &density_fluctuations_file, 
	ofstream &coords_file, ofstream &volume_file, ofstream &temperature_file)
{
	velocity_file << fixed << elapsed_time << " \t";
	density_file << fixed << elapsed_time << " \t";
	velocity_fluctuations_file << fixed << elapsed_time << " \t";
	density_fluctuations_file << fixed << elapsed_time << " \t";
	coords_file << fixed << elapsed_time << " \t";
	for (int i = 0; i < number_of_particles; i++)
	{
		particles[i]->display_velocity(velocity_file);
		particles[i]->display_density(density_file);
		particles[i]->display_velocity_fluctuations(velocity_fluctuations_file, system_velocity);
		particles[i]->display_density_fluctuations(density_fluctuations_file, system_density);
		particles[i]->display_coords(coords_file);
		volume_file << fixed << particles[i]->volume << " \t";
		temperature_file << fixed << temperature[i] << " \t";
	}
	velocity_file << fixed << system_velocity->x() * 1e2 << " \t" << system_velocity->y() * 1e2 << " \t" << system_velocity->z() * 1e2 << " \t" << endl;
	density_file << fixed << system_density * 1.66e3 << " \t" << endl;
	velocity_fluctuations_file << fixed << system_velocity->x() * 1e2 << " \t" << system_velocity->y() * 1e2 << " \t" << system_velocity->z() * 1e2 << " \t" << endl;
	density_fluctuations_file << fixed << system_density * 1.66e3 << " \t" << endl;
	coords_file << endl;
	volume_file << total_particles_volume / number_of_particles << "\t" << endl;
	temperature_file << system_temperature << "\t" << endl;
}

int main(int argc, char * argv[])
{
	//set output to the files
	ofstream velocity_file, density_file, velocity_fluctuations_file, density_fluctuations_file, coords_file, volume_file, temperature_file;
	velocity_file.precision(6);
	density_file.precision(6);
	velocity_fluctuations_file.precision(6);
	density_fluctuations_file.precision(6);
	coords_file.precision(6);
	volume_file.precision(6);
	temperature_file.precision(6);
	velocity_file.open("velocity.txt");
	density_file.open("density.txt");
	velocity_fluctuations_file.open("velocity_fluctuations.txt");
	density_fluctuations_file.open("density_fluctuations.txt");
	coords_file.open("coordinates.txt");
	volume_file.open("volume.txt");
	temperature_file.open("temperature.txt");

	auto counter = 0;

	double ppr = cbrt(number_of_particles);
	double step = box_length / ppr;

	for (double j = step / 2; j < box_length; j += step)
	{
		for (double k = step / 2; k < box_length; k += step)
		{
			for (double p = step / 2; p < box_length; p += step)
			{
				particles[counter] = new HParticle(j, k, p);
				particles[counter]->velocity = *new Point3(0, 0, 0);
				counter++;
			}
		}
	}
	

	//subtract avg system velocity from particles velocities
	subtract_avg_velocity();

	//simulate system changes during n steps
	for (iteration = 0; iteration < 30000; iteration++)
	{
		const auto fade_3d = new Fade_3D();
		triangulation(fade_3d);
		calculate_next_step();
		update_params();
		display_results(velocity_file, density_file, velocity_fluctuations_file, density_fluctuations_file, coords_file, volume_file, temperature_file);
		delete fade_3d;
	}

	//close output files
	velocity_file.close();
	density_file.close();
	velocity_fluctuations_file.close();
	density_fluctuations_file.close();
	coords_file.close();
	volume_file.close();
	temperature_file.close();
}