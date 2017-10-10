//   SVM, version 1.02

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

/*
	Implements the MP-BCFW algorithm for training structural SVMs described in

		Neel Shah, Vladimir Kolmogorov, Christoph H. Lampert
		A Multi-Plane Block-Coordinate Frank-Wolfe Algorithm for Structural SVMs with a Costly max-Oracle
		Technical report arXiv:1408.6804, August 2014

	With parameters cp_max = approx_max = 0 reduces to the BCFW algorithm described in

		S. Lacoste-Julien, M. Jaggi, M. Schmidt, P. Pletscher
		Block-Coordinate Frank-Wolfe Optimization for Structural SVMs
		ICML 2013, Atlanta, USA, June 2013

	If you use this software for research purposes you should cite the aforementioned paper(s) in any resulting publication.

	///////////////////////////////////////////////////////////////////////////////////////////////

	Solves the following problem:
		min  1/2 \lambda ||w||^2 + \mu \sum_{i=1}^n H_i(w)
	where \lambda, \mu > 0 and
		H_i(w) = \max_y <a^{iy},[w \kappa]>
	for a constant \kappa > 0.

	It is assumed H_i(w) can be evaluated efficiently, i.e. the problem
		\max_y <a^{iy},[w \kappa]>
	can be solved efficiently for a given i and w.
*/

#ifndef OAISJNHFOASFASFASFASFNVASF
#define OAISJNHFOASFASFASFASFNVASF

#include "block.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SVMutils.h"
#include "timer.h"


//////////////////////////////
// NOTATION AND CONVENTIONS //
//////////////////////////////
// w: vector of size d
// a or a^{iy}: vector of size d+1
// di: "effective dimension" of term i, i.e. the number of non-zero components of vectors a^{iy} (excluding the free term a^{iy}[d]).
//     Note, di <= d.
// wi: vector of size di    (for term i)
// ai: vector of size di+1  (for term i)
// mapping: determines the correspondence between components of w and wi (and also a and ai).
//          It is an array of size di+1, with mapping[k] \in [0,d-1] for k\in[0,di-1] and mapping[di] = d
//
// In some applications labeling y that determines vector a^{iy} can be stored compactly (with less space than a^{iy}).
// The code assumes that y takes 'y_size_in_bytes' amount of memory. (This excludes the free term a^{iy}[d]).
// Type 'YPtr' is a pointer to such structure; labelings of these type are denoted as y, y1, y2, etc.
// The user must implement a function of type 'CopyFn' transforming labeling y to vector ai of size di (excluding the free term a^{iy}[d]).

typedef void *TermData; // pointer provided by the user in SetTerm()
typedef void *YPtr;


// The user must implement the following four functions. In all of them all variables (wi, ai, y, ...) are assumed to be already allocated.
// Below PAD(wi) denotes vector of size d obtained from we by padding with zeros in appropriate positions.

typedef double (*MaxFn)(double* wi, YPtr y, double kappa, TermData term_data); // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) kappa]> to y, and return the free term a^{iy}[d].
                                                                               // Note, kappa may be different from the value specified in the constructor
typedef void (*CopyFn)(double* ai, YPtr y, TermData term_data); // copies non-zero components of a^{iy} to 'ai' (excluding the free term; note, in this case 'ai' is of size di, not di+1).
                                                                // 'copy_fn' in the constructor can be NULL, then 'y' is exacly the same as ai[0:di-1] and y_size_in_bytes=di*sizeof(double).
typedef bool (*CompareFn)(YPtr y1, YPtr y2, TermData term_data); // returns true if y1==y2
typedef double (*DotProductFn)(double* wi, YPtr y, TermData term_data); // returns <PAD(wi), a^{iy}>. Note, the free term is excluded.
typedef double (*DotProductKernelFn)(YPtr y1, YPtr y2, TermData term_data); // returns <a^{iy1}, a^{iy2}>. Can be NULL, if options.kernel_max == 1.

class SVM;

struct default_SVM_callback {
   bool operator()(SVM* svm);
};


class SVM
{
   friend class default_SVM_callback; // this is not so nice from a design point of view!
public:
	// d = dimension of w, n = # of terms.
	//
	// If zero_lower_bound is true then H_i(w)\ge 0 for all i and w
	// (this is the case for the SVM objective - this lower bound is given by the ground truth labeling y^i).
	// If this flag is on then the algorithm is initialized with the plane (0,...,0,0) for each term,
	// otherwise the initial lower bound plane is obtained by calling the oracle for some w.
	SVM(int d, int n, MaxFn max_fn, CopyFn copy_fn, CompareFn compare_fn, DotProductFn dot_product_fn, DotProductKernelFn dot_product_kernel_fn, bool zero_lower_bound);
	~SVM();


	// 'term_data' will be passed to all user-defined functions.
	// if mapping==NULL, then must have di==d.
	void SetTerm(int i, TermData term_data, int di, int y_size_in_bytes, int* mapping=NULL);
	void SetParams(double lambda, double mu, double kappa);

   template<typename VISITOR>
	double* Solve(VISITOR visitor); // returns a pointer to an array of size 'd' containing solution (vector w).
	                 // For options to Solve(), see SVM::options below
   double* Solve() { return Solve(default_SVM_callback()); }

	void GetBounds(double& lower_bound, double& upper_bound); // of internally stored solution. Expensive - calls n oracles!

	double Evaluate(double* w); // returns the value of the objective function for given w. Expensive - calls n oracles!
   double Evaluate() { return Evaluate(w); }

	// To get block-coordinate Frank-Wolfe, set cp_max = approx_max = 0.
	// ('cp' stands for 'cutting plane')
	struct Options
	{
		Options() :
			randomize_method(2),
			avg_flag(3),

			iter_max(300),
			approx_max(1000), // <--- probably will not be reached (due to the param below)
			approx_limit_ratio(1.0),
			kernel_max(1),

			cp_max(100), // <--- probably will not be reached (due to the param below). Can be decreased if memory is an issue
			cp_inactive_iter_max(10), // <--- PERHAPS THE MOST IMPORTANT PARAMETER:
			                          //      for how many iterations inactive planes are kept in memory

			callback_freq(5),
			callback_fn(default_callback_fn),
			gap_threshold(1e-10), 
			print_flag(2),
			exclude_callback_time(false)
		{
		};

		int randomize_method; // 0: use default order for every iteration (0,1,...,n-1)
		                      // 1: generate a random permutation, use it for every iteration
		                      // 2: generate a new random permutation at every iteration
		                      // 3: generate a new random permutation at every exact & approximate pass
		                      // 4: for every step sample example in {0,1,...,n-1} uniformly at random

		int avg_flag; // 0: don't use averaging
		              // 1: compute weighted average of vectors after each update, as described in (Lacoste-Julien et al. ICML'13)
		              // 2: compute weighted average of vectors after each exact update, as described in (Lacoste-Julien et al. ICML'13)
		              // 3: compute two vectors: avg_exact - weighted average of vectors after each exact update
		              //                         avg_approx - weighted average of vectors after each approx update
		              //    return their best interpolation

		int iter_max;
		int approx_max; // >= 0. Each iter first performs one pass with calls to the 'real' oracle,
		                //       and then up to 'approx_max' passes with calls to the 'approximate' oracle
		                //       It is recommended to set it to a large number and rely on the criterion below.
		double approx_limit_ratio; // extra stopping criterion: approx. pass is stopped if
		                           //    approx_limit_ratio * (increase of the lower bound during B) / (time of B) 
		                           //                       < (increase of the lower bound during A) / (time of A)
		                           // where B corresponds to the last approx. pass and A corresponds to the sequence of steps
		                           // from the beginning of the current iter (including the exact pass) until B

		int kernel_max; // >= 1. Could be helpful if 'd' is very large.
		                // During approximate passes each term is processed 'kernel_max' times.
		                // If >1 then a specialized implementation with kernels is used, where the inner products between different planes are cached.

		///////////////////////////////
		// cutting planes parameters //
		///////////////////////////////

		// If there are more than 'cp_max' planes then remove the plane that has been inactive the longest time.
		// (A plane is active when it is added or when it is returned by the approximate oracle.)
		// Also after each approximate oracle call remove a plane if it hasn't been active during the last 'cp_inactive_iter_max' outer iterations (including the current one)
		int cp_max; // >= 0. 
		int cp_inactive_iter_max;  // if == 0 then this option is not used (so 0 corresponds to +\infty)

		///////////////////////////////////////////////
		// stopping criteria and printed information //
		///////////////////////////////////////////////

		// if callback_fn != NULL then this function will be called after every 'callback_freq' iterations (=callback_freq*n calls to max_fn).
		// If this function returns false then Solve() will terminate.
		// The default function checks the duality gap and prints all bounds.
		bool (*callback_fn)(SVM* svm);

		int callback_freq;
		double gap_threshold;
		int print_flag; // 0: don't print anything
		                // 1: print bounds and gap
		                // 2: print bounds and gap + average # of planes per term + # approx passes in the last outer iter
		bool exclude_callback_time; // when printing time in default_callback_fn().

		static bool default_callback_fn(SVM* svm)
		{
			double lower_bound, upper_bound;
			svm->GetBounds(lower_bound, upper_bound); // this is expensive! (calls real oracles)
			if (lower_bound < svm->lower_bound_last) lower_bound = svm->lower_bound_last; // can happen if averaging is used
			double dual_gap_bound = upper_bound - lower_bound;

			double t = svm->time_from_start; if (svm->options.exclude_callback_time) t -= svm->callback_time;

			if (svm->options.print_flag>0) printf("iterations: %d, time: %.2f sec, lower bound: %f, upper bound: %f, gap %f", svm->iter, t, lower_bound, upper_bound, dual_gap_bound);
			if (svm->options.print_flag>1) printf("   %.1f cp, %d it.", (double)svm->total_plane_num / svm->n, svm->approx_pass);
			printf("\n");

			if (dual_gap_bound < svm->options.gap_threshold) return false;
			return true;
		}
	} options;











//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

private:
	class Term // represents term H(w) = \max_y <a,[w 1]> and the current linear lower bound 'phi'
	{
	public:
	
		Term(int d, TermData term_data, int* mapping, int y_size_in_bytes, SVM* svm, Buffer* buf, bool maintain_products);
		~Term();

		int di;
		int* mapping;
		TermData term_data;
		int y_size_in_bytes_plus;

		int num; // number of planes 'y'
		double* phi; // array of size di+1

#define Y_NOT_ALLOCATED ((YPtr)NULL)

		YPtr* y_arr; // y_arr[t] points to an array of size y_size_in_bytes, 0<=t<num.
		float* last_accessed;  // timestamps, of size num

#define NOT_YET_COMPUTED (3e103) // some random number
		double** products; // products[t1][t2]=dot product of vectors a[t1] and a[t2], ignoring the last (d+1)-th coordinate (0<=t1,t2<num).
						   // valid only if maintain_products was true in the constructor. Computed on demand.

		SVM* svm;

		///////////////////////////////////////////////////////////////////////////////////////

		bool isDuplicate(YPtr y);
		int AddPlane(YPtr y, int cp_max); // if num>=cp_max then the plane with the lowest 'counter' will be deleted
		                                     // and the new plane 'a' will be inserted instead.
		                                     // returns id of the added plane
		void DeletePlane(int t); // plane 'num-1' is moved to position 't'.

		int Maximize(double* wi, double kappa); // returns id of the cutting plane 'a' that maximizes <[wi kappa],a>.

		void UpdateStats(int t); // increases 'counter' for 't' and decreases it for other planes, with parameter 'cp_history_size' (see implementation)

		void RemoveUnusedPlanes();

		double* ComputeRestriction(double* w)
		{
			if (!mapping) return w;
			int k;
			double* wi = svm->wi_buf;
			for (k=0; k<di; k++) wi[k] = w[mapping[k]];
			return wi;
		}

		double* GetFreeTermPtr(YPtr y) { return (double*) ( ((char*)y) + y_size_in_bytes_plus - sizeof(double) ); }

		////////////////////////////////////
	private:
		int num_max;
		Buffer* buf;
		char* my_buf;

		void Allocate(int num_max_new, bool maintain_products);
	};

	int d, n;
	double lambda_mu; // = lambda / mu
	double lambda_mu_inv; // = mu / lambda
	double mu;
	double kappa;
	MaxFn max_fn;
	CopyFn copy_fn;
	CompareFn compare_fn;
	DotProductFn dot_product_fn;
	DotProductKernelFn dot_product_kernel_fn;

	int y_size_in_bytes_max;
	Buffer buf;
	ReusableBuffer rbuf_SolveWithKernel;

	double* w; // of size d
	double* wi_buf; // of size d
	double* neg_phi_sum; // of size d+1
	double neg_phi_sum_norm; // can be NOT_YET_COMPUTED
	double GetNegPhiSumNorm();
	double GetCurrentLowerBound() { return (-GetNegPhiSumNorm()*lambda_mu_inv/2 - neg_phi_sum[d]*kappa)*mu; }
	Term** terms; // of size n

	float timestamp, timestamp_threshold;
	int total_plane_num;
	int iter, approx_pass, total_pass;
	double lower_bound_last;
	double time_start, time_from_start, callback_time;

	bool zero_lower_bound;

	void InitSolver();
	void AddCuttingPlane(int i, YPtr y);
	void SolveWithKernel(int i, int iter_max);
	void InterpolateBest(double* current_sum1, double* current_sum2, double* current_sum_best);
};

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

template<typename VISITOR>
double* SVM::Solve(VISITOR visitor)
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
		if ((iter % options.callback_freq) == 0)
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
				//res = (*options.callback_fn)(this);
            res = visitor(this);
				neg_phi_sum = neg_phi_sum_tmp; w = w_tmp;
			}
			else if (avg_num == 2)
			{
				neg_phi_sum_tmp = neg_phi_sum; w_tmp = w;
				neg_phi_sum = neg_avg_buf; w = neg_avg_buf + d+1;
				InterpolateBest(neg_avg->GetAverage(0), neg_avg->GetAverage(1), neg_phi_sum);
				Multiply(w, neg_phi_sum, lambda_mu_inv, d);
				//res = (*options.callback_fn)(this);
            res = visitor(this);
				neg_phi_sum = neg_phi_sum_tmp; w = w_tmp;
			}
			else
			{
				Multiply(w, neg_phi_sum, lambda_mu_inv, d);
				//res = (*options.callback_fn)(this);
            res = visitor(this);
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

bool default_SVM_callback::operator()(SVM* svm)
{
      double lower_bound, upper_bound;
      svm->GetBounds(lower_bound, upper_bound); // this is expensive! (calls real oracles)
      if (lower_bound < svm->lower_bound_last) lower_bound = svm->lower_bound_last; // can happen if averaging is used
      double dual_gap_bound = upper_bound - lower_bound;

      double t = svm->time_from_start; if (svm->options.exclude_callback_time) t -= svm->callback_time;

      if (svm->options.print_flag>0) printf("iterations: %d, time: %.2f sec, lower bound: %f, upper bound: %f, gap %f", svm->iter, t, lower_bound, upper_bound, dual_gap_bound);
      if (svm->options.print_flag>1) printf("   %.1f cp, %d it.", (double)svm->total_plane_num / svm->n, svm->approx_pass);
      printf("\n");

      if (dual_gap_bound < svm->options.gap_threshold) return false;
      return true; 
   } 


#endif
