#ifndef EURISTIC_SOLVER_TSP

#define EURISTIC_SOLVER_TSP
#include "tsp.h"

void nearest_neighborhood(instance *inst);
void random_GRASP(instance *inst);
double nearest_neighborhood_from_node(instance *inst, int start_node, int *nodes);
double random_GRASP_from_node(instance *inst, int start_node, int *index_neighborhood, int index_neigh_size, int *nodes);
void opt_2(instance *inst, int *tour, int plot);
void opt_2_swap(int *nodes, int start, int end);
void VNS(instance *inst, int random_start);
void make_4_opt_move(int *local_best_tour, int *new_tour, int *swapper, int num_nodes);
void make_5_opt_move(int *local_best_tour, int *new_tour, int *swapper, int num_nodes);
void make_6_opt_move(int *local_best_tour, int *new_tour, int *swapper, int num_nodes);
void simulated_annealing(instance *inst, int random_start);
void build_cplex_sol(instance *inst, double *x_value);

#endif /*EURISTIC_SOLVER_TSP*/