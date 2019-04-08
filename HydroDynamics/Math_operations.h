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