/*
    Copyright Vladimir Kolmogorov vnk@ist.ac.at 2014

    This file is part of SVM.

    SVM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SVM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SVM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef HGALHALSJGHASASFASFG
#define HGALHALSJGHASASFASFG

#include <stdlib.h>

inline int RandomInteger(int N) // returns random integer in [0,N-1]
{
	while ( 1 )
	{
		int x = (int) ( (((double)N) * rand()) / RAND_MAX );
		if (x<N) return x;
	}
}

inline void generate_permutation(int *buf, int n)
{
	int i, j;

	for (i=0; i<n; i++) buf[i] = i;
	for (i=0; i<n-1; i++)
	{
		j = i + RandomInteger(n-i);
		int tmp = buf[i]; buf[i] = buf[j]; buf[j] = tmp;
	}
}


inline void SetZero(double* w, int d)
{
	int k;
	for (k=0; k<d; k++) w[k] = 0;
}

inline double Norm(double* w, int d)
{
	int k;
	double val = 0;
	for (k=0; k<d; k++) val += w[k]*w[k];
	return val;
}

inline double DotProduct(double* u, double* w, int d)
{
	int k;
	double val = 0;
	for (k=0; k<d; k++) val += u[k]*w[k];
	return val;
}

inline double DotProduct(double* u, double* w, int* w_mapping, int di)
{
	if (!w_mapping) return DotProduct(u, w, di);
	int k;
	double val = 0;
	for (k=0; k<di; k++) val += u[k]*w[w_mapping[k]];
	return val;
}

inline void Add(double* u, double* v, int d)
{
	int k;
	for (k=0; k<d; k++) u[k] += v[k];
}

inline void Multiply(double* w, double c, int d)
{
	int k;
	for (k=0; k<d; k++) w[k] *= c;
}

inline void Multiply(double* u, double* v, double c, int d)
{
	int k;
	for (k=0; k<d; k++) u[k] = v[k] * c;
}

// returns <u-v,u>
inline double Op1(double* u, double* v, int d)
{
	int k;
	double val = 0;
	for (k=0; k<d; k++) val += (u[k]-v[k])*u[k];
	return val;
}

// returns <u-v,w>
inline double Op1(double* u, double* v, double* w, int d)
{
	int k;
	double val = 0;
	for (k=0; k<d; k++) val += (u[k]-v[k])*w[k];
	return val;
}

// returns <u-v,w>
inline double Op1(double* u, double* v, double* w, int* w_mapping, int d)
{
	int k;
	double val = 0;
	if (!w_mapping)
	{
		for (k=0; k<d; k++) val += (u[k]-v[k])*w[k];
	}
	else
	{
		for (k=0; k<d; k++) val += (u[k]-v[k])*w[w_mapping[k]];
	}
	return val;
}

// returns ||u-v||^2
inline double Op2(double* u, double* v, int d)
{
	int k;
	double val = 0;
	for (k=0; k<d; k++) val += (u[k]-v[k])*(u[k]-v[k]);
	return val;
}

// set u := (1-gamma)*u + gamma*v
inline void Interpolate(double* u, double* v, double gamma, int d)
{
	int k;
	double beta = 1-gamma;
	for (k=0; k<d; k++) u[k] = beta*u[k] + gamma*v[k];
}

//////////////////////////////////////////




#endif