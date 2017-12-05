LP_MP
========

LP_MP is a C++ framework for developing scalable dual (Lagrangean) decomposition based algorithms solvers for a wide range of LP-relaxations to discrete optimization problems.
For a theoretical introduction to the techniques used and the class of problems that can be optimized see [here](https://arxiv.org/abs/1612.05460).

## Solvers
Solvers are provided in separate projects and include
* **[Discrete graphical models](https://github.com/pawelswoboda/LP_MP-MRF)**,
* **[Multicut](https://github.com/pawelswoboda/LP_MP-Cut)**, 
* **[Graph matching](https://github.com/pawelswoboda/LP_MP-QAP)**, 
* **[Discrete tomography](https://github.com/pawelswoboda/LP_MP-Discrete-tomography)**.

Parallel optimization can be enabled in cmake by setting `LP_MP_PARALLEL` to `ON`.

An interface to external solvers is provided by [DD_ILP](https://github.com/pawelswoboda/DD_ILP).

## Installation
Type `git clone https://github.com/pawelswoboda/LP_MP.git` for downloading, then `cd LP_MP` and `git submodule update --init` for downloading dependencies` and finally `cmake` for building.

Prerequisites:
* Clang 5.0 or GCC 7.0 upwards for C++17 compatibility.

## Contact
* [Paul Swoboda](https://github.com/pawelswoboda)
* [Jan Kuske](https://github.com/DerJFK)
