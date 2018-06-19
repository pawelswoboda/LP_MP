LP_MP
========

[![Build Status](https://travis-ci.org/pawelswoboda/LP_MP.svg?branch=master)](https://travis-ci.org/pawelswoboda/LP_MP)

LP_MP is a C++ framework for developing scalable dual (Lagrangean) decomposition based algorithms solvers for a wide range of LP-relaxations to discrete optimization problems.
For a theoretical introduction to the techniques used and the class of problems that can be optimized see [1].

## Solvers
Solvers are provided in separate projects and include
* **[Discrete graphical models](https://github.com/pawelswoboda/LP_MP-MRF)**,
* **[Multicut](https://github.com/pawelswoboda/LP_MP-Cut)**, 
* **[Graph matching](https://github.com/pawelswoboda/LP_MP-QAP)**, 
* **[Discrete tomography](https://github.com/pawelswoboda/LP_MP-Discrete-tomography)**.

## Optimization techniques
Optimization techniques include
* **Messsage passing [1]**,
* **Subgradient ascent with a proximal bundle method based on the Frank-Wolfe algorithm [2]**, [Vladimir Kolmogorov's](http://http://pub.ist.ac.at/~vnk/) [original implementation](http://pub.ist.ac.at/~vnk/papers/FWMAP.html).
* An interface to **external solvers** is provided by [DD_ILP](https://github.com/pawelswoboda/DD_ILP).


## Installation
Type `git clone https://github.com/pawelswoboda/LP_MP.git` for downloading, then `cd LP_MP` and `git submodule update --init` for downloading dependencies` and finally `cmake` for building.

Prerequisites:
* Clang 5.0 or GCC 7.0 upwards for C++17 compatibility.

## References
* [1]: [`P. Swoboda, J. Kuske and B. Savchynskyy. A Dual Ascent Framework for Lagrangean Decomposition of Combinatorial Problems. In CVPR 2017.`](http://openaccess.thecvf.com/content_cvpr_2017/html/Swoboda_A_Dual_Ascent_CVPR_2017_paper.html)
* [2]: `P. Swoboda and V. Kolmogorov. MAP inference via Block-Coordinate Frank-Wolfe Algorithm. arXiv.`
