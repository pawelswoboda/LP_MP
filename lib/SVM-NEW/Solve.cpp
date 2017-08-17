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
#include "timer.h"


void SVM::AddCuttingPlane(int i, YPtr y)
{
	if (options.cp_max <= 0) return;
	if (terms[i]->isDuplicate(y)) return;
	terms[i]->AddPlane(y, options.cp_max);
}


void SVM::SetTerm(int i, TermData term_data, int di, int y_size_in_bytes, int* mapping)
{
	if (terms[i]) { printf("Error: SetTerm() cannot be called twice for the same term\n"); exit(1); }

	bool maintain_products = (options.kernel_max > 1) ? true : false;
	terms[i] = new Term(di, term_data, mapping, y_size_in_bytes, this, &buf, maintain_products);
}

void SVM::InitSolver()
{
	int i, k;

	SetZero(neg_phi_sum, d+1);
	SetZero(w, d);

	total_plane_num = 0;
	timestamp = 0;
	timestamp_threshold = -1;

	y_size_in_bytes_max = 0;
	for (i=0; i<n; i++)
	{
		if (y_size_in_bytes_max < terms[i]->y_size_in_bytes_plus) y_size_in_bytes_max = terms[i]->y_size_in_bytes_plus;
	}
	YPtr y_buf = (YPtr) new char[y_size_in_bytes_max];

	for (i=0; i<n; i++)
	{
		double* phi = terms[i]->phi;
		if (zero_lower_bound)
		{
			SetZero(phi, d+1);
		}
		else
		{
			double* wi = terms[i]->ComputeRestriction(w);
			YPtr y = (copy_fn) ? y_buf : phi;
			phi[terms[i]->di] = *terms[i]->GetFreeTermPtr(y) = (*max_fn)(wi, y, kappa, terms[i]->term_data);
			if (copy_fn) (*copy_fn)(phi, y, terms[i]->term_data);
			if (!terms[i]->mapping)
			{
				for (k=0; k<=d; k++) neg_phi_sum[k] -= phi[k];
			}
			else
			{
				for (k=0; k<=terms[i]->di; k++) neg_phi_sum[terms[i]->mapping[k]] -= phi[k];
			}
			AddCuttingPlane(i, y);
		}
	}

	delete [] (char*) y_buf;
}


double* SVM::Solve()
{
	time_start = get_time();

	InitSolver();
	int _i, i, k;
	
	// Throughtout the algorithm, we must have w = lambda_mu_inv*neg_phi_sum.
	// For efficiency don't recompute w every time when neg_phi_sum is updated;
	// instead, restore this only when calling external functions and upon termination.

	//Multiply(w, neg_phi_sum, lambda_mu_inv, d);
	neg_phi_sum_norm = NOT_YET_COMPUTED;
	lower_bound_last = GetCurrentLowerBound();

	int di_sum = 0;
	for (i=0; i<n; i++) di_sum += terms[i]->di;

	int approx_max = (options.cp_max <= 0) ? 0 : options.approx_max;
	int avg_num; // maintain 'avg_num' averaged vectors
	switch (options.avg_flag)
	{
		case 0: avg_num = 0; break;
		case 1:
		case 2: avg_num = 1; break;
		default: avg_num = 2; break;
	}

	int vec_size = (d+1)*sizeof(double);
	int alloc_size = (1+avg_num)*vec_size + y_size_in_bytes_max;
	if (options.randomize_method >= 1 || options.randomize_method <= 3)	alloc_size += n*sizeof(int);
	char* _buf = new char[alloc_size];
	char* _buf0 = _buf;

	YPtr y_new_buf = (YPtr) _buf; _buf += y_size_in_bytes_max;
	double* phi_new = (double*) _buf; _buf += vec_size;
	MovingAverage* neg_avg;
	if (di_sum > d*n / 2) neg_avg = new MovingAverageNaive(neg_phi_sum, d, avg_num);
	else                  neg_avg = new MovingAverageLowDim(neg_phi_sum, d, avg_num);
	int* permutation = NULL;
	if (options.randomize_method >= 1 && options.randomize_method <= 3)	{ permutation = (int*) _buf; _buf += n*sizeof(int); }
	double* neg_avg_buf = (double*) _buf;


	int k_avg[2] = { 0, 0 };
	callback_time = 0;

	if (options.randomize_method == 1) generate_permutation(permutation, n);

	for (iter=total_pass=0; iter<options.iter_max; iter++)
	{
		// recompute current_sum every 10 iterations for numerical stability.
		// Just in case, the extra runtime should be negligible.
		// For experiments in the paper this was not used (numerical stability was not an issue)
		if (iter > 0 && (iter % 10) == 0)
		{
			SetZero(neg_phi_sum, d+1);
			for (i=0; i<n; i++)
			{
				if (!terms[i]->mapping)
				{
					for (k=0; k<=terms[i]->di; k++) neg_phi_sum[k] -= terms[i]->phi[k];
				}
				else
				{
					for (k=0; k<=terms[i]->di; k++) neg_phi_sum[terms[i]->mapping[k]] -= terms[i]->phi[k];
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////

		timestamp = (float)(((int)timestamp) + 1); // When a plane is accessed, it is marked with 'timestamp'.
		                                           // Throughout the outer iteration, this counter will be gradually
		                                           // increased from 'iter+1' to 'iter+1.5', so that we
		                                           // (1) we can distinguish between planes added in the same iteration (when removing the oldest plane), and
		                                           // (2) we can easily determine whether a plane has been active during the last 'cp_inactive_iter_max' iterations
		if (options.cp_inactive_iter_max > 0) timestamp_threshold = timestamp - options.cp_inactive_iter_max;

		if (options.randomize_method == 2) generate_permutation(permutation, n);

		double _t[2];           // index 0: before calling real oracle
		double _lower_bound[2]; // index 1: after calling real oracle

		_t[0] = get_time();
		_lower_bound[0] = GetCurrentLowerBound();

		for (approx_pass=-1; approx_pass<approx_max; approx_pass++, total_pass++)
		{
			timestamp += (float) ( 0.5 / (approx_max+1) );

			if (options.randomize_method == 3) generate_permutation(permutation, n);

			for (_i=0; _i<n; _i++)
			{
				if (permutation)                        i = permutation[_i];
				else if (options.randomize_method == 0) i = _i;
				else                                    i = RandomInteger(n);
			
				double* phi = terms[i]->phi;
				double* wi_scaled = terms[i]->ComputeRestriction(neg_phi_sum); // = wi * lambda_mu
				int di = terms[i]->di;
				int* mapping = terms[i]->mapping;
				YPtr y_new = (copy_fn) ? y_new_buf : phi_new;

				// averaging
				double avg_gamma[2] = { 0, 0 };
				neg_avg->StartUpdate(di, mapping);
				int p = (approx_pass < 0 || options.avg_flag == 1) ? 0 : 1;
				if (p < avg_num)
				{
					avg_gamma[p] = 2.0 / (2 + (k_avg[p] ++));
				}

				if (approx_pass < 0) // call real oracle
				{
					*terms[i]->GetFreeTermPtr(y_new) = (*max_fn)(wi_scaled, y_new, kappa*lambda_mu, terms[i]->term_data);
					AddCuttingPlane(i, y_new);
				}
				else  // call approximate oracle
				{
					if (options.kernel_max > 1)
					{
						SolveWithKernel(i, options.kernel_max);
						//Multiply(w, neg_phi_sum, lambda_mu_inv, d);
						terms[i]->RemoveUnusedPlanes();

						neg_avg->FinishUpdate(avg_gamma);

						continue;
					}

					int t = terms[i]->Maximize(wi_scaled, kappa*lambda_mu);
					terms[i]->UpdateStats(t);
					memcpy(y_new, terms[i]->y_arr[t], terms[i]->y_size_in_bytes_plus);
					terms[i]->RemoveUnusedPlanes();
				}
				if (copy_fn) 
            {
               (*copy_fn)(phi_new, y_new, terms[i]->term_data);
               phi_new[di] = *terms[i]->GetFreeTermPtr(y_new); 
            }

				// min_{gamma \in [0,1]} B*gamma*gamma - 2*A*gamma
				double A = -Op1(phi, phi_new, neg_phi_sum, mapping, di) + (phi_new[di] - phi[di])*kappa*lambda_mu; // <phi-phi_new,phi_sum> + (b_new - b) * kappa * lambda
				double B = Op2(phi, phi_new, di); // ||current-current_new||^2
				double gamma;
				if (B<=0) gamma = (A <= 0) ? 0 : 1;
				else
				{
					gamma = A/B;
					if (gamma < 0) gamma = 0;
					if (gamma > 1) gamma = 1;
				}

				if (!mapping)
				{
					for (k=0; k<=di; k++)
					{
						double old = phi[k];
						phi[k] = (1-gamma)*phi[k] + gamma*phi_new[k];
						neg_phi_sum[k] -= phi[k] - old;
					}
				}
				else
				{
					for (k=0; k<=di; k++)
					{
						double old = phi[k];
						phi[k] = (1-gamma)*phi[k] + gamma*phi_new[k];
						neg_phi_sum[mapping[k]] -= phi[k] - old;
					}
				}
				//Multiply(w, neg_phi_sum, lambda_mu_inv, d);

				// averaging
				neg_avg->FinishUpdate(avg_gamma);
			}
			neg_phi_sum_norm = NOT_YET_COMPUTED;

			double t = get_time();
			lower_bound_last = GetCurrentLowerBound();

			if (approx_pass >= 0)
			{
				if ( (lower_bound_last - _lower_bound[1]) * (_t[1]-_t[0]) * options.approx_limit_ratio
				   < (_lower_bound[1]  - _lower_bound[0]) * (t-_t[1])      ) { approx_pass ++; break; }
			}

			_t[1] = t;
			_lower_bound[1] = lower_bound_last;
		}

		time_from_start = get_time() - time_start;
		if (options.callback_fn && (iter % options.callback_freq) == 0)
		{
			double t0 = get_time();

			bool res;
			double* w_tmp;
			double* neg_phi_sum_tmp;

			if (avg_num == 1)
			{
				neg_phi_sum_tmp = neg_phi_sum; w_tmp = w;
				neg_phi_sum = neg_avg->GetAverage(0); w = neg_avg_buf;
				Multiply(w, neg_phi_sum, lambda_mu_inv, d);
				res = (*options.callback_fn)(this);
				neg_phi_sum = neg_phi_sum_tmp; w = w_tmp;
			}
			else if (avg_num == 2)
			{
				neg_phi_sum_tmp = neg_phi_sum; w_tmp = w;
				neg_phi_sum = neg_avg_buf; w = neg_avg_buf + d+1;
				InterpolateBest(neg_avg->GetAverage(0), neg_avg->GetAverage(1), neg_phi_sum);
				Multiply(w, neg_phi_sum, lambda_mu_inv, d);
				res = (*options.callback_fn)(this);
				neg_phi_sum = neg_phi_sum_tmp; w = w_tmp;
			}
			else
			{
				Multiply(w, neg_phi_sum, lambda_mu_inv, d);
				res = (*options.callback_fn)(this);
			}

			callback_time += get_time() - t0;

			if (!res) break;
		}
	}

	if (avg_num == 1)
	{
		memcpy(neg_phi_sum, neg_avg->GetAverage(0), vec_size);
		Multiply(w, neg_phi_sum, lambda_mu_inv, d);
	}
	else if (avg_num == 2)
	{
		InterpolateBest(neg_avg->GetAverage(0), neg_avg->GetAverage(1), neg_phi_sum);
		Multiply(w, neg_phi_sum, lambda_mu_inv, d);
	}

	delete [] _buf0;
	delete neg_avg;
	return w;
}

void SVM::SolveWithKernel(int _i, int iter_max)
{
	Term* T = terms[_i];
	int num = T->num, di = T->di, i, t, iter;
	double** kk = T->products;
	double* ck = (double*) rbuf_SolveWithKernel.Alloc(3*num*sizeof(double)); // ck[i] = DotProduct(phi, T->a[i], d)
	double* sk = ck + num; // sk[i] = DotProduct(phi_sum, T->a[i], d)
	double cc = DotProduct(T->phi, T->phi, di);
	double cs = -DotProduct(T->phi, neg_phi_sum, T->mapping, di);
	double c_d = T->phi[di]*kappa;

	double* x = sk + num;
	double cx = 1;

	double gamma;

	double* neg_phi_sum_short = T->ComputeRestriction(neg_phi_sum);
	for (i=0; i<num; i++)
	{
		ck[i] = (*dot_product_fn)(T->phi,             T->y_arr[i], T->term_data);
		sk[i] = -(*dot_product_fn)(neg_phi_sum_short, T->y_arr[i], T->term_data);
		x[i] = 0;
	}

	for (iter=0; iter<iter_max; iter++)
	{
		if (iter > 0)
		{
			c_d += gamma*(*T->GetFreeTermPtr(T->y_arr[t])*kappa - c_d);
			double cc_new = cc + 2*gamma*(ck[t] - cc) + gamma*gamma*(kk[t][t] - 2*ck[t] + cc);
			double cs_new = cs+ gamma*(ck[t] + sk[t] - cc - cs) + gamma*gamma*(kk[t][t] - 2*ck[t] + cc);
			cc = cc_new;
			cs = cs_new;

			for (i=0; i<num; i++)
			{
				if (kk[i][t] == NOT_YET_COMPUTED) kk[i][t] = (*dot_product_kernel_fn)(T->y_arr[i], T->y_arr[t], T->term_data);
				double delta = gamma*(kk[i][t] - ck[i]);
				ck[i] += delta;
				sk[i] += delta;
			}
		}

		t = 0;
		double v_best;
		for (i=0; i<num; i++)
		{
			double v = -sk[i] * lambda_mu_inv + *T->GetFreeTermPtr(T->y_arr[i])*kappa;
			if (i==0 || v_best < v) { t = i; v_best = v; }
		}
		T->UpdateStats(t);

		if (kk[t][t] == NOT_YET_COMPUTED) kk[t][t] = (*dot_product_kernel_fn)(T->y_arr[t], T->y_arr[t], T->term_data);

		// min_{gamma \in [0,1]} B*gamma*gamma - 2*A*gamma
		double A = cs - sk[t] + (*T->GetFreeTermPtr(T->y_arr[t])*kappa - c_d)*lambda_mu;
		double B = cc + kk[t][t] - 2*ck[t];

		if (B<=0) gamma = (A <= 0) ? 0 : 1;
		else
		{
			gamma = A/B;
			if (gamma < 0) gamma = 0;
			if (gamma > 1) gamma = 1;
		}

		cx *= 1-gamma;
		for (i=0; i<num; i++) x[i] *= 1-gamma;
		x[t] += gamma;
	}

	if (!T->mapping)
	{
		for (i=0; i<=di; i++) { neg_phi_sum[i] += T->phi[i]; T->phi[i] *= cx; }
	}
	else
	{
		for (i=0; i<=di; i++) { neg_phi_sum[T->mapping[i]] += T->phi[i]; T->phi[i] *= cx; }
	}

	double* a;
	if (copy_fn) a = new double[di+1];
	for (t=0; t<num; t++)
	{
		if (x[t] == 0) continue;
		if (copy_fn) { (*copy_fn)(a, T->y_arr[t], T->term_data); a[di] = *T->GetFreeTermPtr(T->y_arr[t]); }
		else a = (double*) T->y_arr[t];
		for (i=0; i<=di; i++) T->phi[i] += x[t]*a[i];
	}
	if (copy_fn) delete [] a;

	if (!T->mapping)
	{
		for (i=0; i<=di; i++) neg_phi_sum[i] -= T->phi[i];
	}
	else
	{
		for (i=0; i<=di; i++) neg_phi_sum[T->mapping[i]] -= T->phi[i];
	}
}

void SVM::InterpolateBest(double* s1, double* s2, double* s_best)
{
	double A = Op1(s1, s2, d) - (s2[d] - s1[d])*kappa*lambda_mu;
	double B = Op2(s1, s2, d);
	double gamma;
	if (B<=0) gamma = (A <= 0) ? 0 : 1;
	else
	{
		gamma = A/B;
		if (gamma < 0) gamma = 0;
		if (gamma > 1) gamma = 1;
	}

	int k;
	for (k=0; k<=d; k++)
	{
		s_best[k] = (1-gamma)*s1[k] + gamma*s2[k];
	}
}

