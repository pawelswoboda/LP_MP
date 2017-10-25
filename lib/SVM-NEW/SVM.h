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
		min  1/2 \sum_{i=1}^n H_i(w)
	where 
		H_i(w) = \max_y <a^{iy},[w 1]>

	It is assumed H_i(w) can be evaluated efficiently, i.e. the problem
		\max_y <a^{iy},[w 1]>
	can be solved efficiently for a given i and w.
*/

#ifndef OAISJNHFOASFASFASFASFNVASF
#define OAISJNHFOASFASFASFASFNVASF

#include "block.h"

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
// Below PAD(wi) denotes vector of size d obtained from wi by padding with zeros in appropriate positions.

typedef double (*MaxFn)(double* wi, YPtr y, TermData term_data); // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) 1]> to y, and return the free term a^{iy}[d].
typedef void (*CopyFn)(double* ai, YPtr y, TermData term_data); // copies non-zero components of a^{iy} to 'ai' (excluding the free term; note, in this case 'ai' is of size di, not di+1).
                                                                // 'copy_fn' in the constructor can be NULL, then 'y' is exacly the same as ai[0:di-1] and y_size_in_bytes=di*sizeof(double).
typedef bool (*CompareFn)(YPtr y1, YPtr y2, TermData term_data); // returns true if y1==y2
typedef double (*DotProductFn)(double* wi, YPtr y, TermData term_data); // returns <PAD(wi), a^{iy}>. Note, the free term is excluded.
typedef double (*DotProductKernelFn)(YPtr y1, YPtr y2, TermData term_data); // returns <a^{iy1}, a^{iy2}>. Can be NULL, if options.kernel_max == 1.

class SVM
{
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

	double* Solve(); // returns a pointer to an array of size 'd' containing solution (vector w).
	                 // For options to Solve(), see SVM::options below

	//void GetBounds(double& lower_bound, double& upper_bound); // of internally stored solution. Expensive - calls n oracles!

	double Evaluate(double* w); // returns the value of the objective function for given w. Expensive - calls n oracles!

	// To get block-coordinate Frank-Wolfe, set cp_max = approx_max = 0.
	// ('cp' stands for 'cutting plane')
	struct Options
	{
		Options() :
			randomize_method(2),
			avg_flag(0),

			iter_max(100000),
			approx_max(1000), // <--- probably will not be reached (due to the param below)
			approx_limit_ratio(1.0),
			kernel_max(1),

			cp_max(100), // <--- probably will not be reached (due to the param below). Can be decreased if memory is an issue
			cp_inactive_iter_max(10), // <--- PERHAPS THE MOST IMPORTANT PARAMETER:
			                          //      for how many iterations inactive planes are kept in memory

			callback_freq(5),
			c_min(1e-5),
			c_max(1e5),
			c_init(1),
			c_increase_factor(2.0),
			c_decrease_factor(0.5),

			serious_iter_min(200),
			serious_iter_max(10),
			serious_iter_init(50),
			serious_iter_increase_factor(1.1),
			serious_iter_decrease_factor(0.9),

			gap_threshold(1e-10), 
			print_flag(2)
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
		              // (NOT IMPLEMENTED) 3: compute two vectors: avg_exact - weighted average of vectors after each exact update
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
		// If this function returns 0: Solve() will terminate.
		//                         -1: update w and continue
		//                          1: continue
		// TODO    The default function checks the duality gap and prints all bounds.

		int callback_freq;
		double c_min;
		double c_max;
		double c_init;
		double c_increase_factor;
		double c_decrease_factor;

		int serious_iter_min;
		int serious_iter_max;
		int serious_iter_init;
		double serious_iter_increase_factor;
		double serious_iter_decrease_factor;


		double gap_threshold;
		int print_flag; // 0: don't print anything
		                // 1: print bounds and gap
		                // 2: print bounds and gap + average # of planes per term + # approx passes in the last outer iter

		//static int default_callback_fn(SVM* svm)
		//{
		//	printf("%d, %.2f sec:   %f %f, ||phi*||=%f\n", svm->iter, svm->time_from_start, svm->GetCurrentLowerBound(), svm->Evaluate(svm->z), svm->GetPhiStarNorm());
		//	if (svm->iter_since_last_update > 200)
		//	{
		//		return -1;
		//	}
		//	return 1;
			/*
			double lower_bound, upper_bound;
			svm->GetBounds(lower_bound, upper_bound); // this is expensive! (calls real oracles)
			if (lower_bound < svm->lower_bound_last) lower_bound = svm->lower_bound_last; // can happen if averaging is used
			double dual_gap_bound = upper_bound - lower_bound;

			double t = svm->time_from_start; if (svm->options.exclude_callback_time) t -= svm->callback_time;

			if (svm->options.print_flag>0) printf("%d, %.2f sec:   %f %f, gap %f", svm->iter, t, lower_bound, upper_bound, dual_gap_bound);
			if (svm->options.print_flag>1) printf("   %.1f cp, %d it.", (double)svm->total_plane_num / svm->n, svm->approx_pass);
			printf("\n");

			if (dual_gap_bound < svm->options.gap_threshold) return false;
			return true;
			*/
		//}
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

		int Maximize(double* wi); // returns id of the cutting plane 'a' that maximizes <[wi 1],a>.

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
	double c;
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
	double* z; // of size d+1. = [w_prox - c*phi_sum, -c*phi_sum[d]]
	//double neg_phi_sum_norm; // can be NOT_YET_COMPUTED
	//double GetNegPhiSumNorm();
	double GetCurrentLowerBound()
	{
		int i;
		double sum = 0;
		for (i=0; i<d; i++) sum += (w[i] - z[i])*(w[i] + z[i]);
		return ( 0.5*sum - z[d] ) / c;
	}
	double GetPhiStarNorm()
	{
		int i;
		double sum = 0;
		for (i=0; i<d; i++) sum += (w[i]-z[i])*(w[i]-z[i]);
		return sum / (c*c);
	}
	Term** terms; // of size n

	float timestamp, timestamp_threshold;
	int total_plane_num;
	int iter, approx_pass, total_pass;
	double lower_bound_last;
	double time_start;

	bool zero_lower_bound;

	void InitSolver();
	void AddCuttingPlane(int i, YPtr y);
	void SolveWithKernel(int i, int iter_max);
	void InterpolateBest(double* current_sum1, double* current_sum2, double* current_sum_best);
};



#endif
