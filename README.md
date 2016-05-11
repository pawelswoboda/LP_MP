# LP_MP
Solving LPs arising from discrete optimization with convergent message passing

Installation: 
Required:
- Install Vc from https://github.com/VcDevel/Vc/releases, take version 1.2.0
- Edit CmakeLists file and set Vc-paths correctly

Optional:
- HDF5 for reading datasets, in debian: apt-get install libhdf5-dev. Also download opengm for reading graphical models then into lib
- for graph matching: Download lemon-1.3.1 into lib
- Gurobi: Compile libgurobi_c++ with clang, if clang compiles LP_MP
