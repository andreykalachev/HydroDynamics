#pragma once
#include <cmath>

static int positive_mod(int x, int y) {
	return (x % y + y) % y;
}

static double positive_mod(double x, double y) {
	return fmod(fmod(x, y) + y, y);
}

static double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

static double calculate_absolute_value(Point3 point)
{
	return sqrt(pow(point.x(), 2) + pow(point.y(), 2) + pow(point.z(), 2));
}

static Point3 calculate_sum(Point3 point1, Point3 point2)
{
	auto x = point1.x() + point2.x();
	auto y = point1.y() + point2.y();
	auto z = point1.z() + point2.z();

	return Point3(x, y, z);
}

static Point3 subtract(Point3 point1, Point3 point2)
{
	auto x = point1.x() - point2.x();
	auto y = point1.y() - point2.y();
	auto z = point1.z() - point2.z();

	return Point3(x, y, z);
}

static Point3 divide(Point3 point1, double divider)
{
	auto x = point1.x() / divider;
	auto y = point1.y() / divider;
	auto z = point1.z() / divider;

	return Point3(x, y, z);
}


static bool compare_coords(vector<Tet3*>::value_type& tet, int tet_corner, HParticle* particle)
{
	return tet->getCorner(tet_corner)->x() == particle->coordinates.x() && tet->getCorner(tet_corner)->y() == particle->coordinates.y() &&
		tet->getCorner(tet_corner)->z() == particle->coordinates.z();
}