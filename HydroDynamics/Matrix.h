#pragma once
#include <random>
#include <time.h>
#include <chrono>

struct Matrix
{
	double matrix[3][3];

	Matrix() {}

	Matrix(double m[3][3])
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				matrix[i][j] = m[i][j];
			}
		}
	}

	void generateRandom(std::default_random_engine generator, std::normal_distribution<double> distribution)
	{
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		generator.seed(seed);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				matrix[i][j] = distribution(generator);
			}
		}
	}

	Matrix getTransported()
	{
		Matrix result;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result.matrix[i][j] = matrix[j][i];
			}
		}
		return result;
	}

	double getDiagonalSum()
	{
		double sum = 0;
		for (int i = 0; i < 3; i++)
		{
			sum += matrix[i][i];
		}
		return sum;
	}
};

static Matrix operator+(Matrix matrix1, Matrix matrix2)
{
	Matrix result;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.matrix[i][j] = matrix1.matrix[i][j] + matrix2.matrix[i][j];
		}
	}
	return result;
}

static Matrix operator-(Matrix matrix1, Matrix matrix2)
{
	Matrix result;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.matrix[i][j] = matrix1.matrix[i][j] - matrix2.matrix[i][j];
		}
	}
	return result;
}

static Matrix operator*(Matrix matrix, double factor)
{
	Matrix result;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.matrix[i][j] = factor * matrix.matrix[i][j];
		}
	}
	return result;
}

static Matrix operator*(double factor, Matrix matrix)
{
	return matrix * factor;
}

static Matrix operator/(Matrix matrix, double factor)
{
	Matrix result;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.matrix[i][j] = matrix.matrix[i][j] / factor;
		}
	}
	return result;
}

static Matrix multiply_by_identity_matrix(double factor)
{
	double result[3][3] = { {factor, 0, 0}, {0, factor, 0}, {0, 0, factor} };
	return Matrix(result);
}

static FADE3D::Point3 operator*(FADE3D::Point3 vector, Matrix matrix)
{
	return FADE3D::Point3
	(
		matrix.matrix[0][0] + matrix.matrix[1][0] + matrix.matrix[2][0],
		matrix.matrix[0][1] + matrix.matrix[1][1] + matrix.matrix[2][1],
		matrix.matrix[0][2] + matrix.matrix[1][2] + matrix.matrix[2][2]
	);
}