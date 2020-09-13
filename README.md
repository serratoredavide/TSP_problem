# TSP_problem

### Ricerca Operativa 2 project report
University of Padua | Computer Engineering Department


#### Authors
Serratore Davide - M.1207660 | Vesco Omar - M.1197699

---

## Repository content

	RO2_report.pdf                \\ Project report
	data/*.tsp                    \\ TSPLIB instances used in the report http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/index.html

	Makefile                      \\ Makefile to compile code
	src/compact_models_TSP.c      \\ Compact models for exact resolution
	src/heuristic_models_TSP.c    \\ Heuristic models CPLEX-based
	src/heuristic_solver_TSP.c    \\ Heuristic algorithms not based on CPLEX
	src/tsp.c                     \\ General structure of the problem and some useful function, loop and lazy approach
	src/utils.c                   \\ Utility functions
	include/*.h                   \\ Headers
	main.c                        \\ Main code 
	a280.tsp                      \\ An instance from TSPLIB (280 nodes) http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/



## Compile and Execution

To compile (ILOG CPLEX needed):

	make 

To execute code. Some parameters are needed:

	./exec_tsp -file data/xxxxx.tsp 

Parameters:

	-help | --help          \\ for help
	-file | -f | -input     \\ input instances
	-time_limit             \\ Time limit in seconds. For VNS : number of iterations
	-solver_type            \\ integer defining the type of model (use: ./exec_tsp -help )
	-init                   \\ initialization type for heuristic solver without cplex. 0 for Nearest Neighborhood, 1 for Random GRASP


VERBOSE is also defined in tsp.h. To plot tour: VERBOSE >= 80

`make clean` to clean up executable files
