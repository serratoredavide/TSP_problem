#include "heuristic_models_TSP.h"

#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "heuristic_solver_TSP.h"
#include "mystring.h"
#include "tsp.h"
#include "utils.h"

/**
 * @brief Set the seed for random choice
 * 
 * @param seed 
 */
void set_seed(int seed) {
    if (seed == -1)
        srand(time(NULL));
    else
        srand(seed);
}

/**
 * @brief Get random number between 0 and 1
 * 
 * @return double value between 0 and 1
 */
double get_rand() {
    return (double)rand() / (double)RAND_MAX;
}

/**
 * @brief Fix some edges randomly
 * 
 * @param inst Pointer to instance struct 
 * @param x_value Array solution returned by CPLEX (they represent edges)
 * @param new_lb Array of lower bounds 
 * @param index Array that indicates edges to be modify
 * @param lb_char Array of 'L' needed for CPLEX
 * @param env CPLEX environment
 * @param lp CPLEX data problem
 */
void fix_variable(const instance *inst, const double *x_value, const double fixing_prob, double *new_lb, int *index, char *lb_char, CPXENVptr env, CPXLPptr lp) {
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;

    int modified = 0;
    for (int i = 0; i < num_x_columns; i++) {
        if (x_value[i] > 0.5)
            if (get_rand() > 1 - fixing_prob) {
                index[modified] = i;
                modified++;
            }
    }
    if (CPXtightenbds(env, lp, modified, index, lb_char, new_lb)) print_error(ERROR_LOWER_BOUND_FIX_VAR_STRING);
    CPXwriteprob(env, lp, MODEL_NAME, NULL);
}

/**
 * @brief Perform Hard Fixing for the instance 
 * 
 * @param inst Pointer to instance struct
 * @return int 0 if successful and nonzero if an error occurs
 */
int hard_fixing_TSP(instance *inst) {
    double local_time = 600.0;
    if (inst->timelimit == -1.0) print_error(TIME_LIMIT_NOT_FOUND_STRING);
    if (inst->timelimit < local_time) print_error("time_limit parameter too low. Required at least 600 secs");
    set_seed(inst->randomseed);  //set seed for fixing probability

    time_t start, end;
    int time_taken;

    //init cplex environment
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    init_cplex_env(&env, &lp, inst);
    init_env_callback(env, lp, inst);
    //local_time for CPLEX
    if (CPXsetdblparam(env, CPX_PARAM_TILIM, local_time)) print_error(ERROR_SETTING_TIMELIM_STRING);

    //useful value
    int mipopt_code, getx_code;
    int iterations = 0;
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;
    double best_bound = CPX_INFBOUND;
    double obj_value = CPX_INFBOUND;

    build_model(inst, env, lp);

    //useful data
    double *x_value = (double *)calloc(num_x_columns, sizeof(double));
    //default lower bound, upper bound value and char
    double *default_lb = (double *)calloc(num_x_columns, sizeof(double));
    double *default_ub = (double *)malloc(num_x_columns * sizeof(double));
    char *lb_char = (char *)malloc(num_x_columns * sizeof(char));
    char *ub_char = (char *)malloc(num_x_columns * sizeof(char));
    int *default_index = (int *)malloc(num_x_columns * sizeof(int));
    //index array used to fix variable
    int *new_index = (int *)malloc(num_x_columns * sizeof(int));
    //for initial solution
    int beg = 0;
    int *index = malloc(num_nodes * sizeof(int));
    double *values = malloc(num_nodes * sizeof(double));
    char **start_names = calloc(1, sizeof(char *));
    start_names[0] = "Grasp + 2-Opt initialization";

    //init of some arrays
    for (int i = 0; i < num_x_columns; i++) {
        default_ub[i] = 1.0;
        default_index[i] = i;
        lb_char[i] = LOWER_BOUND_CHAR;
        ub_char[i] = UPPER_BOUND_CHAR;
    }

    init_heuristic_inst(inst);
    time(&start);
    //generate initial solution
    random_GRASP(inst);
    opt_2(inst, inst->nodes_tour, 1);
    //solution to CPLEX
    for (int i = 0; i < num_nodes; i++) {
        index[i] = triangle_xpos(inst->nodes_tour[i], inst->nodes_tour[(i + 1) % num_nodes], inst);
        values[i] = 1.0;
    }
    CPXaddmipstarts(env, lp, 1, num_nodes, &beg, index, values, NULL, start_names);

    mipopt_code = CPXmipopt(env, lp);  //CPLEX with callback is used
    if (mipopt_code) return mipopt_code;
    time(&end);

    //print obj_value and bound
    if (CPXgetobjval(env, lp, &obj_value)) print_error(ERROR_OBJVALUE_STRING);
    printf("Solution value without fixing: %.5f\n", obj_value);
    if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
    printf("Best Bound without fixing : %.5f\n", best_bound);

    double fixing_prob = 0.90;
    double decr_value = 0.20;
    double last_value = obj_value;
    while (end - start < inst->timelimit && fixing_prob > 0.0) {
        getx_code = CPXgetx(env, lp, x_value, 0, num_x_columns - 1);
        if (getx_code) print_error(ERROR_GET_X_STRING);

        if (iterations > 0) {
            //reset bounds of variables x_ij
            if (CPXtightenbds(env, lp, num_x_columns, default_index, lb_char, default_lb)) print_error(ERROR_DEFAULT_LB_STRING);
            if (CPXtightenbds(env, lp, num_x_columns, default_index, ub_char, default_ub)) print_error(ERROR_DEFAULT_UP_STRING);
        }
        fix_variable(inst, x_value, fixing_prob, default_ub, new_index, lb_char, env, lp);
        mipopt_code = CPXmipopt(env, lp);  //CPLEX with callback is used
        if (mipopt_code) return mipopt_code;

        if (CPXgetobjval(env, lp, &obj_value)) print_error(ERROR_OBJVALUE_STRING);
        if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
        printf("Solution value with fixing prob %.2f: %.5f\n", fixing_prob, obj_value);
        printf("Best Bound with fixing prob %.2f: %.5f\n", fixing_prob, best_bound);
        //check for improvment
        if (last_value <= obj_value)
            fixing_prob -= decr_value;
        else
            last_value = obj_value;

        iterations++;
        time(&end);
    }

    //last call to CPLEX
    if (end - start < inst->timelimit) {
        //reset bounds of variables x_ij
        if (CPXtightenbds(env, lp, num_x_columns, default_index, lb_char, default_lb)) print_error(ERROR_DEFAULT_LB_STRING);
        if (CPXtightenbds(env, lp, num_x_columns, default_index, ub_char, default_ub)) print_error(ERROR_DEFAULT_UP_STRING);

        //set to CPLEX remaining time
        double remaining_time = inst->timelimit - (end - start);
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time)) print_error(ERROR_SETTING_TIMELIM_STRING);
        mipopt_code = CPXmipopt(env, lp);  //CPLEX with callback is used
        if (mipopt_code) return mipopt_code;

        if (CPXgetobjval(env, lp, &obj_value)) print_error(ERROR_OBJVALUE_STRING);
        if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
        printf("Solution value with fixing prob %.2f: %.5f\n", 0.0, obj_value);
        printf("Best Bound with fixing prob %.2f: %.5f\n", 0.0, best_bound);
        time(&end);
    }

    time_taken = end - start;
    printf("Computation time: %d seconds\n", time_taken);

    // plotting
    FILE *gnuplotpipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) {
        //get x solution
        int getx_code = CPXgetx(env, lp, x_value, 0, num_x_columns - 1);
        if (getx_code) print_error(ERROR_GET_X_STRING);

        solution *sol = (solution *)malloc(sizeof(solution));
        //malloc solution struct
        sol->successor = (int *)malloc(num_nodes * sizeof(int));
        sol->component = (int *)malloc(num_nodes * sizeof(int));

        //build solution to plot
        build_solution(x_value, inst, sol);

        gnuplotpipe = open_init_plot();
        temp = fopen("points", "w");
        createPointsFile(inst, sol, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotpipe, plot, 1);
        fclose(temp);
        fflush(gnuplotpipe);
        pclose(gnuplotpipe);

        //free memory
        free_solution(sol);
        free(sol);
    }

    //free memory
    free(x_value);
    free(default_lb);
    free(default_ub);
    free(default_index);
    free(new_index);
    free(lb_char);
    free(ub_char);
    free(index);
    free(values);
    free(start_names);
    //free and close cplex model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

int soft_fixing_TSP(instance *inst) {
    double local_time = 600.0;
    if (inst->timelimit == -1.0) print_error(TIME_LIMIT_NOT_FOUND_STRING);
    if (inst->timelimit < local_time) print_error("time_limit parameter too low. Required at least 600 secs");

    time_t start, end;
    int time_taken;

    // open cplex model
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    init_cplex_env(&env, &lp, inst);
    init_env_callback(env, lp, inst);
    //local_time for CPLEX
    if (CPXsetdblparam(env, CPX_PARAM_TILIM, local_time)) print_error(ERROR_SETTING_TIMELIM_STRING);

    //useful value
    int mipopt_code, getx_code;
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;
    double best_bound = CPX_INFBOUND;
    double obj_value = CPX_INFBOUND;

    build_model(inst, env, lp);

    //useful data
    double *x_value = (double *)calloc(num_x_columns, sizeof(double));
    int *index = (int *)malloc(num_nodes * sizeof(int));
    double *value = (double *)malloc(num_nodes * sizeof(double));
    char *name = "soft_constraint";
    //for initial solution
    int beg = 0;
    char **start_names = calloc(1, sizeof(char *));
    start_names[0] = "Grasp + 2-Opt initialization";

    init_heuristic_inst(inst);
    time(&start);
    //generate initial solution
    random_GRASP(inst);
    opt_2(inst, inst->nodes_tour, 1);
    //solution to CPLEX
    for (int i = 0; i < num_nodes; i++) {
        index[i] = triangle_xpos(inst->nodes_tour[i], inst->nodes_tour[(i + 1) % num_nodes], inst);
        value[i] = 1.0;
    }
    CPXaddmipstarts(env, lp, 1, num_nodes, &beg, index, value, NULL, start_names);
    mipopt_code = CPXmipopt(env, lp);

    if (mipopt_code) return mipopt_code;
    time(&end);

    if (CPXgetobjval(env, lp, &obj_value)) print_error(ERROR_OBJVALUE_STRING);
    printf("Solution value without fixing: %.5f\n", obj_value);
    if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
    printf("Best Bound without fixing : %.5f\n", best_bound);

    int iterations = 0;
    // int k[5] = {2,3,5,10,20};
    int k = 2;
    int incr_value = 2;
    double last_value = obj_value;
    double rhs = num_nodes - k;
    char sense = GREATER_EQUAL_CONSTRAINT_PARAMETER;
    // int beg = 0;
    int row_to_delete = CPXgetnumrows(env, lp) - 1;
    while (end - start < inst->timelimit && k < num_nodes) {
        //get object value
        if (CPXgetobjval(env, lp, &obj_value)) print_error(ERROR_OBJVALUE_STRING);
        //get solution
        getx_code = CPXgetx(env, lp, x_value, 0, num_x_columns - 1);
        if (getx_code) print_error(ERROR_GET_X_STRING);

        //delete last row
        if (iterations > 0)
            if (CPXdelrows(env, lp, row_to_delete, row_to_delete)) print_error(ERROR_DEL_ROW_STRING);

        int count = 0;
        for (int i = 0; i < num_x_columns; i++) {
            if (x_value[i] > 0.5) {
                index[count] = i;
                value[count] = 1.0;
                count++;
            }
        }
        rhs = num_nodes - k;
        //add new row
        CPXaddrows(env, lp, 0, 1, num_nodes, &rhs, &sense, &beg, index, value, NULL, &name);

        //solve problem
        mipopt_code = CPXmipopt(env, lp);
        if (mipopt_code) return mipopt_code;
        if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
        printf("Solution value with k %d: %.5f\n", k, obj_value);
        printf("Best Bound with k %d: %.5f\n", k, best_bound);
        //check for improvment
        if (last_value <= obj_value)
            k += incr_value;
        else
            last_value = obj_value;

        iterations++;
        time(&end);
    }

    //last call to CPLEX
    if (end - start < inst->timelimit) {
        //delete last row
        if (CPXdelrows(env, lp, row_to_delete, row_to_delete)) print_error(ERROR_DEL_ROW_STRING);

        //set to CPLEX remaining time
        double remaining_time = inst->timelimit - (end - start);
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time)) print_error(ERROR_SETTING_TIMELIM_STRING);

        //solve problem
        mipopt_code = CPXmipopt(env, lp);
        if (mipopt_code) return mipopt_code;
        if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
        printf("Solution value with k %d: %.5f\n", num_nodes, obj_value);
        printf("Best Bound with k %d: %.5f\n", num_nodes, best_bound);
    }

    time(&end);
    time_taken = end - start;
    printf("Computation time: %d seconds\n", time_taken);

    // plotting
    FILE *gnuplotpipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) {
        //get x solution
        int getx_code = CPXgetx(env, lp, x_value, 0, num_x_columns - 1);
        if (getx_code) print_error(ERROR_GET_X_STRING);

        solution *sol = (solution *)malloc(sizeof(solution));
        //malloc solution struct
        sol->successor = (int *)malloc(num_nodes * sizeof(int));
        sol->component = (int *)malloc(num_nodes * sizeof(int));

        //build solution to plot
        build_solution(x_value, inst, sol);

        gnuplotpipe = open_init_plot();
        temp = fopen("points", "w");
        createPointsFile(inst, sol, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotpipe, plot, 1);
        fclose(temp);
        fflush(gnuplotpipe);
        pclose(gnuplotpipe);

        //free memory
        free_solution(sol);
        free(sol);
    }

    free(x_value);
    free(index);
    free(value);
    free(start_names);
    //free and close cplex model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}