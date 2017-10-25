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


class MovingAverage
{
public:
	// x is a pointer to an array of size d+1 that is modified externally. 
	// 
	//MovingAverage(double* x, int d)=0; // set avg = x
	//~MovingAverage()=0;

	virtual void StartUpdate(int di, int* mapping)=0; // notification that x is going to change soon
	virtual void FinishUpdate(double* gamma)=0; // set avg := (1-gamma)*avg + gamma*x

	virtual double* GetAverage(int r)=0;
};
	
//////////////////////////////////////////

class MovingAverageNaive : public MovingAverage
{
public:
	// x is a pointer to an array of size d+1 that is modified externally. 
	// 
	MovingAverageNaive(double* x, int d, int avg_num); // set avg = x
	~MovingAverageNaive();

	void StartUpdate(int di, int* mapping) {} // notification that x is going to change soon
	void FinishUpdate(double* gamma); // set avg := (1-gamma)*avg + gamma*x

	double* GetAverage(int r) { return avg[r]; }

private:
	int d, avg_num;
	double* x;
	double** avg;
};

inline MovingAverageNaive::MovingAverageNaive(double* _x, int _d, int _avg_num)
	: x(_x), d(_d), avg_num(_avg_num)
{
	int r;
	avg = new double*[avg_num];
	avg[0] = new double[(avg_num)*(d+1)];
	for (r=0; r<avg_num; r++)
	{
		if (r) avg[r] = avg[r-1] + d+1;
		memcpy(avg[r], x, (d+1)*sizeof(double));
	}
}

inline MovingAverageNaive::~MovingAverageNaive()
{
	delete [] avg[0];
	delete [] avg;
}

inline void MovingAverageNaive::FinishUpdate(double* gamma)
{
	int r, k;
	for (r=0; r<avg_num; r++)
	{
		if (gamma[r] == 0) continue;

		for (k=0; k<=d; k++)
		{
			avg[r][k] = (1-gamma[r])*avg[r][k] + gamma[r]*x[k];
		}
	}
}

//////////////////////////////////////////

#undef MOVING_AVERAGE_TEST

class MovingAverageLowDim : public MovingAverage
{
public:
	// x is a pointer to an array of size d+1 that is modified externally. 
	// 
	MovingAverageLowDim(double* x, int d, int avg_num); // set avg = x
	~MovingAverageLowDim();

	void StartUpdate(int di, int* mapping); // notification that x is going to change soon
	void FinishUpdate(double* gamma); // set avg := (1-gamma)*avg + gamma*x

	double* GetAverage(int r);

private:
	// avg[r] = x + alpha[r]*z[r]
	int d, avg_num;
	double* x;
	double* x_old;
	double** avg;
	double** z;
	double* alpha;

	int di;
	int* mapping;
#ifdef MOVING_AVERAGE_TEST
	MovingAverageNaive test;
#endif
};

inline MovingAverageLowDim::MovingAverageLowDim(double* _x, int _d, int _avg_num)
	: x(_x), d(_d), avg_num(_avg_num)
#ifdef MOVING_AVERAGE_TEST
	, test(_x, _d, _avg_num)
#endif
{
	int r;
	avg = new double*[2*avg_num];
	z = avg + avg_num;
	alpha = new double[2*(avg_num)*(d+1) + (d+1) + avg_num];
	x_old = alpha + avg_num;
	for (r=0; r<avg_num; r++)
	{
		avg[r] = ((r) ? z[r-1] : x_old) + (d+1);
		z[r] = avg[r] + d+1;
		memset(z[r], 0, (d+1)*sizeof(double));
		alpha[r] = 1;
	}
}

inline MovingAverageLowDim::~MovingAverageLowDim()
{
	delete [] alpha;
	delete [] avg;
}

inline void MovingAverageLowDim::StartUpdate(int _di, int* _mapping)
{
#ifdef MOVING_AVERAGE_TEST
	test.StartUpdate(_di, _mapping);
#endif

	di = _di;
	mapping = _mapping;
	int k;
	if (!mapping)
	{
		for (k=0; k<=di; k++) x_old[k] = x[k];
	}
	else
	{
		for (k=0; k<=di; k++) x_old[k] = x[mapping[k]];
	}
}

inline void MovingAverageLowDim::FinishUpdate(double* gamma)
{
#ifdef MOVING_AVERAGE_TEST
	test.FinishUpdate(gamma);
#endif

	int r, k;
	if (!mapping)
	{
		for (r=0; r<avg_num; r++)
		{
			double alpha_inv = 1/alpha[r];
			alpha[r] *= 1-gamma[r];
			for (k=0; k<=d; k++) z[r][k] -= alpha_inv*(x[k] - x_old[k]);
		}
	}
	else
	{
		for (r=0; r<avg_num; r++)
		{
			double alpha_inv = 1/alpha[r];
			alpha[r] *= 1-gamma[r];
			for (k=0; k<=di; k++) z[r][mapping[k]] -= alpha_inv*(x[mapping[k]] - x_old[k]);
		}
	}
	for (r=0; r<avg_num; r++)
	{
		if (alpha[r] < 1e-6)
		{
			for (k=0; k<=d; k++) z[r][k] *= alpha[r];
			alpha[r] = 1;
		}
	}

#ifdef MOVING_AVERAGE_TEST
	for (r=0; r<avg_num; r++)
	{
		double* avg_test = test.GetAverage(r);
		for (k=0; k<=d; k++)
		{
			double v = x[k] + alpha[r]*z[r][k] - avg_test[k];
			if (v < 0) v = -v;
			if (v > 1e-10)
			{
				printf("MovingAverage error\n");
				exit(1);
			}
		}
	}
#endif
}

inline double* MovingAverageLowDim::GetAverage(int r)
{
	int k;
	for (k=0; k<=d; k++)
	{
		z[r][k] *= alpha[r];
		avg[r][k] = x[k] + z[r][k];
	}
	alpha[r] = 1;
	return avg[r];
}


#endif