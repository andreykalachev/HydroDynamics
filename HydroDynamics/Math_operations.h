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

static bool compare_coords(vector<Tet3*>::value_type& tet, int tet_corner, HParticle* particle)
{
	return  *tet->getCorner(tet_corner) == particle->coordinates;
}

static bool has_corner(vector<Tet3*>::value_type& tet, HParticle* particle, int* out_param_corner)
{
	for (int corner = 0; corner < 4; corner++)
	{
		if (compare_coords(tet, corner, particle)) {
			*out_param_corner = corner;
			return true;
		}
	}
	return false;
}
/*Operations with Point3*/
static Point3 operator-(Point3 point1, Point3 point2)
{
	return Point3(point1.x() - point2.x(), point1.y() - point2.y(), point1.z() - point2.z());
}

static Point3 operator+(Point3 point1, Point3 point2)
{
	return Point3(point1.x() + point2.x(), point1.y() + point2.y(), point1.z() + point2.z());
}

static double operator*(Point3 point1, Point3 point2)
{
	return point1.x() * point2.x() + point1.y() * point2.y() + point1.z() * point2.z();
}

static Point3 operator/(Point3 point1, double divider)
{
	return Point3(point1.x() / divider, point1.y() / divider, point1.z() / divider);
}

static Point3 operator*(Point3 point1, double factor)
{
	return Point3(point1.x() * factor, point1.y() * factor, point1.z() * factor);
}

static Point3 operator*(double factor, Point3 point1)
{
	return Point3(point1.x() * factor, point1.y() * factor, point1.z() * factor);
}

static Point3 operator%(Point3 point1, double divider)
{
	return Point3(positive_mod(point1.x(), divider), positive_mod(point1.y(), divider), positive_mod(point1.z(), divider));
}

static void operator-=(Point3 &point1, Point3 point2)
{
	point1.init(Point3(point1.x() - point2.x(), point1.y() - point2.y(), point1.z() - point2.z()));
}

static void operator+=(Point3 &point1, Point3 point2)
{
	point1.init(Point3(point1.x() + point2.x(), point1.y() + point2.y(), point1.z() + point2.z()));
}

static void operator/=(Point3 &point1, double divider)
{
	point1.init(Point3(point1.x() / divider, point1.y() / divider, point1.z() / divider));
}

static void operator*=(Point3 &point1, double factor)
{
	point1.init(Point3(point1.x() * factor, point1.y() * factor, point1.z() * factor));
}

static void operator%=(Point3 &point1, double divider)
{
	point1.init(Point3(positive_mod(point1.x(), divider), positive_mod(point1.y(), divider), positive_mod(point1.z(), divider)));
}

static bool operator!=(Point3 point1, Point3 point2)
{
	return point1.x() != point2.x() || point1.y() != point2.y() || point1.z() != point2.z();
}

static double distance(const Point3& p0, const Point3& p1)
{
	return sqrt(pow(p1.x() - p0.x(), 2) + pow(p1.y() - p0.y(), 2) + pow(p1.z() - p0.z(), 2));
}

static Point3 mod(Point3 point)
{
	return Point3(fabs(point.x()), fabs(point.y()), fabs(point.z()));
}


//additional operations
static HParticle* createRandomParticle(double min, double max)
{
	auto x = fRand(min, max);
	auto y = fRand(min, max);
	auto z = fRand(min, max);

	return new HParticle(fRand(min, max), fRand(min, max), fRand(min, max));
}

static Point3 createRandomPoint(double min, double max)
{
	auto x = fRand(min, max);
	auto y = fRand(min, max);
	auto z = fRand(min, max);

	return Point3(fRand(min, max), fRand(min, max), fRand(min, max));
}


static double calculate_area(Point3 _1, Point3 _2, Point3 _3)
{
	auto a = calculate_absolute_value(_1 - _2);
	auto b = calculate_absolute_value(_2 - _3);
	auto c = calculate_absolute_value(_3 - _1);

	auto p = (a + b + c) / 2.0;

	return sqrt(p * (p - a) * (p - b) * (p - c));
}