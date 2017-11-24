LP_MP
========

LP_MP is a C++ framework for developing scalable monotonuously improving dual block coordinate ascent (also known as message passing) solvers for a wide range of LP-relaxations to discrete optimization problems.
For a theoretical introduction to the techniques used and the class of problems that can be optimized see [here](https://arxiv.org/abs/1612.05460).

## Solvers
Solvers are provided in other projects and include
* **[Discrete graphical models](https://github.com/pawelswoboda/LP_MP-MRF)**,
* **[Multicut](https://github.com/pawelswoboda/LP_MP-Cut)**, 
* **[Graph matching](https://github.com/pawelswoboda/LP_MP-QAP)**, 
* **[Discrete tomography](https://github.com/pawelswoboda/LP_MP-Discrete-tomography)**.

Parallel optimization can be enabled in cmake by setting `LP_MP_PARALLEL` to `ON`.

SAT-based rounding can be enabled for some problems by setting `WITH_SAT_BASED_ROUNDING` to `ON`.

## Installation
Type `git clone https://github.com/pawelswoboda/LP_MP.git` for downloading and `cmake` for building.

Prerequisites:
* Clang 3.8 or GCC 5.4 upwards
* Sqlite3 (for evaluation)

## Contact
* [Paul Swoboda](https://github.com/pawelswoboda)
* [Jan Kuske](https://github.com/DerJFK)
