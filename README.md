LP_MP
========

LP_MP is a C++ framework for developing scalable convergent message passing solvers for a wide range of relaxations to discrete optimization problems.

## Solvers
We provide solvers for the following optimization problems:
* Discrete graphical models: message passing as done by TRWS, SRMP or MPLP, input in opengm (hdf5) or [uai](http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php) format (text).
* (Lifted) Multicut, with cycle and odd wheel inequality separation, input in opengm, a custom format used by [Andres' graph package](https://github.com/bjoern-andres/graph) (both hdf5) and a simple text format.
* Graph matching, input accepted in the format as used by the [dual decomposition graph matching solver of Vladimir Kolmogorov](http://pub.ist.ac.at/~vnk/software/GraphMatching-v1.02.src.zip) or in a custom uai format (both text).
* Multi-label discrete tomography, with input accpeted in a custom uai format (text).

Additionally, an interface to the (I)LP-solver [gurobi](http://www.gurobi.com) is available for solving the above optimization problems. In this case the message passing solvers can act as pre-solvers and initial bound providers, improving performance of the subsequent optimization performed by gurobi.

A large number of datasets can be automatically downloaded for evaluating solvers.

## Installation
type `git clone https://github.com/pawelswoboda/LP_MP.git` to download and `cmake` for building.

Prerequisites:
* Clang 3.8
* HDF5 (for inputs in hdf5 format)
* Gurobi (for the gurobi interface. Compile libgurobi_c++ with clang)
* Sqlite3 (for evaluation)

