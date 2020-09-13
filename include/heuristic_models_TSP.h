#ifndef EURISTIC_MODELS_TSP

#define EURISTIC_MODELS_TSP
#include <cplex.h>
// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>

#include "tsp.h"

void set_seed(int seed);
double get_rand();
void fix_variable(const instance *inst, const double *x_value, const double fixing_prob, double *new_lb, int *index, char *lb_char, CPXENVptr env, CPXLPptr lp);
int hard_fixing_TSP(instance *inst);
int soft_fixing_TSP(instance *inst);

#endif /*EURISTIC_MODELS_TSP*/
