#ifndef COMPACT_MODELS_TSP

#define COMPACT_MODELS_TSP
#include <cplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tsp.h"

int compact_TSPopt(const instance *inst, const int model_type, const int enable_lazy);
void atsp_init_model(const instance *inst, CPXENVptr env, CPXLPptr lp);
void build_atsp_solution(const double *x_value, const instance *inst, solution *solution);
void MTZ_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp);
void lazy_MTZ_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp);
void FLOW1_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp);
void lazy_FLOW1_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp);
int matrix_xpos(int i, int j, const instance *inst);

#endif /*COMPACT_MODELS_TSP*/
