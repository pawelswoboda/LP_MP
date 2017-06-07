LP_MP
========


LP_MP is a C++ framework for developing scalable convergent message passing solvers for a wide range of LP-relaxations to discrete optimization problems.
For a theoretical introduction to the techniques used and the class of problems that can be optimized see [here](https://arxiv.org/abs/1612.05460).

## Solvers
We provide solvers for the following optimization problems:
* **Discrete graphical models**: message passing as done by TRWS, SRMP or MPLP, input in opengm (hdf5) or [uai](http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php) format (text). Tightening with [frustrated cycles](http://cs.nyu.edu/~dsontag/papers/sontag_uai12.pdf) is also supported.
* **Multicut**, (also known as correlation clustering) with cycle<!--, odd wheel and odd bicycle wheel --> and odd wheel inequality separation, input in opengm, a custom format used by [Andres' graph package](https://github.com/bjoern-andres/graph) (both hdf5) and a simple text format. Primal solutions are obtained with the efficient [Kernighan&Lin algorithm](https://github.com/bjoern-andres/graph).
* **Lifted multicut**, which additionally has lifted edges consituting soft connectivity priors. Input formats, inequalities and primal heuristic are similar to multicut.
<!---* **(Asymmetric) Multiway cut with input in the opengm format. -->
* **Graph matching**, input accepted in the format as used by the [dual decomposition graph matching solver of Vladimir Kolmogorov](http://pub.ist.ac.at/~vnk/software/GraphMatching-v1.02.src.zip) or in a custom uai format (both text).
* **Multi-label discrete tomography**, with input accpeted in a custom uai format (text).
<!---* **Tracking by detection** for some cell-tracking problems, with input in a custom text format.-->

<!---*
Additionally, interfaces to the (I)LP-solvers [gurobi](http://www.gurobi.com) and [cplex](http://www.ibm.com/software/integration/optimization/cplex-optimizer/) are available for solving the above optimization problems. In this case the message passing solvers can act as pre-solvers and initial bound providers, improving performance of the subsequent optimization performed by the LP-solvers.
-->

Parallel optimization can be enabled in cmake by setting `LP_MP_PARALLEL` to `ON`.

SAT-based rounding can be enabled for some problems by setting `WITH_SAT_BASED_ROUNDING` to `ON`.

A large number of datasets can be automatically downloaded for evaluating solvers.

## Installation
Type `git clone https://github.com/pawelswoboda/LP_MP.git` for downloading and `cmake` for building.

Prerequisites:
* Clang 3.8 or GCC 5.4 upwards
* HDF5 (for inputs in hdf5 format)
* Gurobi (for the gurobi interface. Compile libgurobi_c++ with clang if you use clang for compiling LP_MP)
* CPLEX (for the cplex interface)
* Sqlite3 (for evaluation)

## Contact
* [Paul Swoboda](https://github.com/pawelswoboda)
* [Jan Kuske](https://github.com/DerJFK)
