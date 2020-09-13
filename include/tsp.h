#ifndef TSP_H_

#define TSP_H_

#include <cplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// #include <pthread.h>

#define VERBOSE 100  // printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//model type
#define TSP_LOOP 0
#define TSP_LAZY_CALLBACK 1
#define COMPACT_TSP_MTZ 2
#define COMPACT_TSP_LAZY_MTZ 3
#define COMPACT_TSP_FLOW1 4
#define COMPACT_TSP_LAZY_FLOW1 5
#define HARD_FIXING_TSP 6
#define SOFT_FIXING_TSP 7
#define OPT_2_TSP 8
#define VNS_TSP 9
#define SA_TSP 10

//data structures

/**
 * @brief Instance Structure
 */
typedef struct {
    //input data
    int nnodes;      // number of nodes
    double *xcoord;  // x coords of nodes
    double *ycoord;  // y coords of nodes

    // parameters
    int model_type;         // what model to use to solve problem
    int randomseed;         // random seed
    double timelimit;       // overall time limit, in sec.s
    int initialization;     // initialization for heuristic_solver: 0 for NN, 1 for randomGRASP
    char input_file[1000];  // input file

    //heuristic data
    int *nodes_tour;   // tour of heuristic model not based on cplex
    double *distance;  // distances between nodes used in heuristic_solver_TSP.c
} instance;

/**
 * @brief Record of the solution 
 * successor : next node
 * component : index of the subtour
 * ncomp : number of subtours in the solution
 */
typedef struct {
    int *successor;
    int *component;
    int ncomp;
} solution;

void free_instance(instance *inst);
void free_solution(solution *sol);
double dist(int i, int j, const instance *inst);
int TSPopt(const instance *inst);
int TSPopt_callback(instance *inst);
int CPXPUBLIC mycallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user);
int addcut(instance *inst, double *x_value, CPXCALLBACKCONTEXTptr context);
void build_model(const instance *inst, CPXENVptr env, CPXLPptr lp);
int triangle_xpos(int i, int j, const instance *inst);
void build_solution(const double *x_value, const instance *inst, solution *solution);

void init_env_callback(CPXENVptr env, CPXLPptr lp, instance *inst);
void init_cplex_env(CPXENVptr *env, CPXLPptr *lp, const instance *inst);
void init_heuristic_inst(instance *inst);

#endif /* TSP_H_ */
