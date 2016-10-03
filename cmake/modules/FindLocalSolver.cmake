FIND_LIBRARY(LOCALSOLVER_LIB
  NAMES liblocalsolver
  HINTS /opt
  PATHS /usr/lib
        ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH)
        
FIND_PATH(LOCALSOLVER_INCLUDE_DIR localsolver.h)