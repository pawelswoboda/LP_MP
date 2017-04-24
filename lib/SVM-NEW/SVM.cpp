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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SVM.h"
#include "SVMutils.h"



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

SVM::SVM(int _d, int _n, MaxFn _max_fn, CopyFn _copy_fn, CompareFn _compare_fn, DotProductFn _dot_product_fn, DotProductKernelFn _dot_product_kernel_fn, bool _zero_lower_bound) :
	d(_d), n(_n), 
	max_fn(_max_fn), copy_fn(_copy_fn), compare_fn(_compare_fn), dot_product_fn(_dot_product_fn), dot_product_kernel_fn(_dot_product_kernel_fn),
	zero_lower_bound(_zero_lower_bound),
	neg_phi_sum_norm(NOT_YET_COMPUTED),
	buf(1024)
{
	int i;
	w = (double*) buf.Alloc(d*sizeof(double));
	wi_buf = (double*) buf.Alloc(d*sizeof(double));
	terms = new Term*[n];
	for (i=0; i<n; i++) terms[i] = NULL;
	neg_phi_sum = (double*) buf.Alloc((d+1)*sizeof(double));
}

SVM::~SVM()
{
	if (terms)
	{
		int i;
		for (i=0; i<n; i++)
		{
			if (terms[i]) delete terms[i];
		}
		delete [] terms;
	}
}

void SVM::SetParams(double _lambda, double _mu, double _kappa)
{
	mu = _mu;
	kappa = _kappa;

	lambda_mu = _lambda/mu;
	lambda_mu_inv = mu/_lambda;
}

double SVM::GetNegPhiSumNorm()
{
	if (neg_phi_sum_norm == NOT_YET_COMPUTED) neg_phi_sum_norm = Norm(neg_phi_sum, d);
	return neg_phi_sum_norm;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double SVM::Evaluate(double* _w)
{
	int i;

	double* w_bak = new double[d+1];
	memcpy(w_bak, w, d*sizeof(double));
	if (w != _w) memcpy(w, _w, d*sizeof(double));
	YPtr y = new char[y_size_in_bytes_max];

	double v0 = Norm(w, d), v1 = 0;

	for (i=0; i<n; i++)
	{
		double* wi = terms[i]->ComputeRestriction(w);
		v1 += (*max_fn)(wi, y, kappa, terms[i]->term_data)*kappa;
		v1 += (*dot_product_fn)(wi, y, terms[i]->term_data);
	}

	memcpy(w, w_bak, d*sizeof(double));

	delete [] w_bak;
	delete [] (char*)y;

	return mu * (v0*lambda_mu/2 + v1);
}


void SVM::GetBounds(double& lower_bound, double& upper_bound)
{
	int i;

	YPtr y = new char[y_size_in_bytes_max];

	double norm = Norm(neg_phi_sum, d)*lambda_mu_inv/2;
	lower_bound = -norm - neg_phi_sum[d]*kappa;
	upper_bound = 0;
	for (i=0; i<n; i++)
	{
		double* wi_scaled = terms[i]->ComputeRestriction(neg_phi_sum); // = wi * lambda_mu
		upper_bound += (*max_fn)(wi_scaled, y, kappa*lambda_mu, terms[i]->term_data)*kappa;
		upper_bound += (*dot_product_fn)(wi_scaled, y, terms[i]->term_data)*lambda_mu_inv;
	}
	upper_bound += norm;

	lower_bound *= mu;
	upper_bound *= mu;

	delete [] (char*)y;
}








/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
////////////////////// Implementation of 'Term' /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////





SVM::Term::Term(int _d, TermData _term_data, int* _mapping, int _y_size_in_bytes, SVM* _svm, Buffer* _buf, bool maintain_products)
	: di(_d), term_data(_term_data), mapping(_mapping), svm(_svm), buf(_buf), num(0), num_max(0)
{
	y_size_in_bytes_plus = _y_size_in_bytes + sizeof(double);
	phi = (double*) buf->Alloc((di+1)*sizeof(double));
	y_arr = NULL;
	last_accessed = NULL;
	products = NULL;
	my_buf = NULL;
	//Allocate(svm->options.cp_max, maintain_products);
	Allocate(4, maintain_products); // start with up to 4 planes per term, then allocate more if necessary
}

void SVM::Term::Allocate(int num_max_new, bool maintain_products)
{
	int num_max_old = num_max;
	YPtr* y_arr_old = y_arr;
	float* last_accessed_old = last_accessed;
	double** products_old = products;
	char* my_buf_old = my_buf;
	num_max = num_max_new;

	int i, my_buf_size = num_max*sizeof(YPtr) + num_max*sizeof(float);
	if (maintain_products) my_buf_size += num_max*sizeof(double*) + num_max*num_max*sizeof(double);
	my_buf = new char[my_buf_size];

	y_arr = (YPtr*) my_buf;
	for (i=0; i<num_max_old; i++) y_arr[i] = y_arr_old[i];
	for ( ; i<num_max; i++) y_arr[i] = Y_NOT_ALLOCATED;

	last_accessed = (float*) (y_arr + num_max);
	memcpy(last_accessed, last_accessed_old, num_max_old*sizeof(int));

	if (maintain_products)
	{
		products = (double**)(last_accessed+num_max);
		for (i=0; i<num_max; i++)
		{
			products[i] = (i==0) ? ((double*)(products+num_max)) : (products[i-1] + num_max);
		}
		int t1, t2;
		for (t1=0; t1<num; t1++)
		for (t2=0; t2<num; t2++)
		{
			products[t1][t2] = products_old[t1][t2];
		}
	}

	if (my_buf_old) delete [] my_buf_old;
}


SVM::Term::~Term()
{
	if (my_buf) delete [] my_buf;
}

bool SVM::Term::isDuplicate(YPtr y)
{
	int t;

	for (t=0; t<num; t++)
	{
		if ((*svm->compare_fn)(y, y_arr[t], term_data))
		{
			last_accessed[t] = svm->timestamp;
			return true;
		}
	}
	return false;
}

int SVM::Term::AddPlane(YPtr y, int cp_max)
{
	int t, t2;

	if (num >= cp_max)
	{
		for (t=0, t2=1; t2<num; t2++)
		{
			if (last_accessed[t] > last_accessed[t2]) t = t2;
		}
	}
	else
	{
		if (num >= num_max)
		{
			int num_max_new = 2*num_max+1; if (num_max_new > cp_max) num_max_new = cp_max;
			Allocate(num_max_new, (products) ? true : false);
		}
		t = num ++;
		if (y_arr[t] == Y_NOT_ALLOCATED) y_arr[t] = (YPtr) buf->Alloc(y_size_in_bytes_plus);
		svm->total_plane_num ++;
	}
	memcpy(y_arr[t], y, y_size_in_bytes_plus);
	last_accessed[t] = svm->timestamp;

	if (products)
	{
		for (t2=0; t2<num; t2++)
		{
			products[t][t2] = products[t2][t] = NOT_YET_COMPUTED; // DotProduct(a[t], a[t2], d);
		}
	}

	return t;
}

void SVM::Term::DeletePlane(int t)
{
	num --;
	svm->total_plane_num --;
	if (t == num) return;
	YPtr tmp = y_arr[t]; y_arr[t] = y_arr[num]; y_arr[num] = tmp;
	last_accessed[t] = last_accessed[num];

	if (products)
	{
		int t2;
		for (t2=0; t2<num; t2++)
		{
			products[t][t2] = products[t2][t] = products[num][t2];
		}
		products[t][t] = products[num][num];
	}
}

int SVM::Term::Maximize(double* wi, double kappa)
{
	int t_best, t;
	double v_best;
	for (t=0; t<num; t++)
	{
		double v = (*svm->dot_product_fn)(wi, y_arr[t], term_data) + (*GetFreeTermPtr(y_arr[t]))*kappa;
		if (t == 0 || v_best <= v) { v_best = v; t_best = t; }
	}
	return t_best;
}

void SVM::Term::UpdateStats(int t_best)
{
	last_accessed[t_best] = svm->timestamp;
}

void SVM::Term::RemoveUnusedPlanes()
{
	if (svm->timestamp_threshold < 0) return;

	int t;
	for (t=0; t<num; t++)
	{
		if (last_accessed[t] < svm->timestamp_threshold && num > 1)
		{
			DeletePlane(t --);
		}
	}
}
