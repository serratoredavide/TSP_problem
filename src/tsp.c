#include "tsp.h"

#include <unistd.h>
// #include <time.h>

#include "mystring.h"
#include "utils.h"

//#define DEBUG

/**
 * @brief Free memory
 * 
 * @param inst Pointer to instance struct
 */
void free_instance(instance *inst) {
    free(inst->xcoord);
    free(inst->ycoord);
    if (inst->nodes_tour != NULL) free(inst->nodes_tour);
    if (inst->distance != NULL) free(inst->distance);
}
/**
 * @brief Free memory
 * 
 * @param sol Pointer to solution struct
 */
void free_solution(solution *sol) {
    free(sol->successor);
    free(sol->component);
}

/**
 * @brief Compute the distance using Euclidian distance
 * 
 * @param i First node index
 * @param j Second node index
 * @param inst Pointer to struct
 * @return double : Distance value
 */
double dist(int i, int j, const instance *inst) {
    double dx = inst->xcoord[i] - inst->xcoord[j];
    double dy = inst->ycoord[i] - inst->ycoord[j];
    return sqrt(dx * dx + dy * dy);
}

/**
 * @brief Compute the optimal solution using loop model
 * 
 * @param inst Instance struct
 * @param sol Solution struct
 * @return int 0 if successful and nonzero if an error occurs.
 */
int TSPopt(const instance *inst) {
    // open cplex model
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    init_cplex_env(&env, &lp, inst);

    //initialize model and useful data
    build_model(inst, env, lp);
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;
    double *x_value = (double *)calloc(num_x_columns, sizeof(double));
    char **rows_name = (char **)calloc(1, sizeof(char *));
    rows_name[0] = (char *)calloc(100, sizeof(char));
    solution *sol = (solution *)malloc(sizeof(solution));
    time_t start, end;
    int time_taken;

    //malloc solution struct
    sol->successor = (int *)malloc(num_nodes * sizeof(int));
    sol->component = (int *)malloc(num_nodes * sizeof(int));

    FILE *gnuplotpipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) gnuplotpipe = open_init_plot();

    //we want ncomp = 1
    int iteration = 1;
    int done = 0;
    int already_plotted = 0;
    int number_subtour_elimination = 0;
    double global_time = (inst->timelimit == -1.0) ? CPX_INFBOUND : (inst->timelimit);
    time(&start);
    end = start;
    while (!done && (end - start) < global_time) {
        done = 1;
        if (inst->timelimit != 1.0)
            if (CPXsetdblparam(env, CPX_PARAM_TILIM, global_time - (end - start))) print_error(ERROR_SETTING_TIMELIM_STRING);
        // Solve the problem
        int mipopt_code = CPXmipopt(env, lp);
        if (mipopt_code) return mipopt_code;

        //get x solution
        int getx_code = CPXgetx(env, lp, x_value, 0, num_x_columns - 1);
        if (getx_code) print_error(ERROR_GET_X_STRING);

        //build solution to plot
        build_solution(x_value, inst, sol);
        if (VERBOSE >= VERBOSE_LITTLE_OUTPUT) {
            printf("Problem solved: iteration %d\n", iteration);
            iteration++;
            printf("Subtour found %d\n", sol->ncomp);
        }

        //Plotting
        if (VERBOSE >= VERBOSE_TO_PLOT) {
            temp = fopen("points", "w");
            createPointsFile(inst, sol, temp);
            if (!already_plotted) {
                char *plot[] = {"plot 'points' with lines"};
                sendCommands(gnuplotpipe, plot, 1);
                already_plotted++;
            } else {
                char *replot[] = {"replot"};
                sendCommands(gnuplotpipe, replot, 1);
            }
            fclose(temp);
            fflush(gnuplotpipe);
            sleep(1);
        }

        if (sol->ncomp >= 2) {
            done = 0;
            //update constraints
            //initialize useful structure
            double *rhs = (double *)calloc((sol->ncomp), sizeof(double));
            char sense = LESS_EQUAL_CONSTRAINT_PARAMETER;  //'L';
            //last row before subtour
            int lastrow_before_st = CPXgetnumrows(env, lp);

            //count number of nodes for each component
            for (int i = 0; i < num_nodes; i++) {
                rhs[sol->component[i]] += 1;
            }
            //update rhs,set sense to <= and set name to the new row
            for (int i = 0; i < sol->ncomp; i++) {
                rhs[i] -= 1;
                number_subtour_elimination++;
                sprintf(rows_name[0], "subtour_constraint %d", number_subtour_elimination);
                CPXnewrows(env, lp, 1, &rhs[i], &sense, NULL, rows_name);
            }

            for (int i = 0; i < num_nodes; i++) {
                for (int j = i + 1; j < num_nodes; j++) {
                    if (sol->component[i] == sol->component[j]) {
                        CPXchgcoef(env, lp, lastrow_before_st + sol->component[i], triangle_xpos(sol->successor[i], sol->successor[j], inst), 1.0);
                    }
                }
            }

            free(rhs);
        }
        time(&end);
    }
    time(&end);
    time_taken = end - start;
    printf("Computation time: %d seconds\n", time_taken);

    if (VERBOSE >= VERBOSE_TO_PLOT) pclose(gnuplotpipe);
    if (VERBOSE >= 100) CPXwriteprob(env, lp, MODEL_NAME, NULL);
    if (VERBOSE >= 10) {
        double best_solution = CPX_INFBOUND;
        if (CPXgetobjval(env, lp, &best_solution)) print_error(ERROR_OBJVALUE_STRING);
        printf("Solution found\n");
        printf("Best Object value: %f\n", best_solution);
    }

    //free pointer and close cplex model
    free(x_value);
    free(rows_name[0]);
    free(rows_name);
    free_solution(sol);
    free(sol);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

/**
 * @brief Compute the optimal solution using lazy callback
 * 
 * @param inst Instance struct
 * @return int 0 if successful and nonzero if an error occurs
 */
int TSPopt_callback(instance *inst) {
    // open cplex model
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    init_cplex_env(&env, &lp, inst);
    init_env_callback(env, lp, inst);

    //initialize model and useful data
    build_model(inst, env, lp);
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;
    time_t start, end;
    int time_taken;

    time(&start);
    //solve problem
    int mipopt_code = CPXmipopt(env, lp);
    if (mipopt_code) return mipopt_code;
    time(&end);
    time_taken = end - start;
    printf("Computation time: %d seconds\n", time_taken);

    if (VERBOSE >= 10) {
        double best_solution = CPX_INFBOUND;
        double best_bound = CPX_INFBOUND;
        if (CPXgetobjval(env, lp, &best_solution)) print_error(ERROR_OBJVALUE_STRING);
        if (CPXgetbestobjval(env, lp, &best_bound)) print_error(ERROR_BEST_OBJVALUE_STRING);
        printf("Solution found\n");
        printf("Best Object value: %f\n", best_solution);
        printf("Best bound: %f\n", best_bound);

        FILE *gnuplotpipe = NULL;
        FILE *temp = NULL;
        //plotting
        if (VERBOSE >= VERBOSE_TO_PLOT) {
            //get x solution
            double *x_value = (double *)calloc(num_x_columns, sizeof(double));
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
            free(x_value);
            free_solution(sol);
            free(sol);
        }
    }

    //free memory and close cplex model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

/**
 * @brief Function called by cplex on lazy callback
 * 
 * @param context Callback context structure
 * @param contextid Id to identify  the contest
 * @param user User input. It should be Instance struct
 * @return int mycallback 0 if successful and nonzero if an error occurs
 */
int CPXPUBLIC mycallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user) {
    instance *inst = (instance *)user;
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;

    double *x_value = (double *)malloc(num_x_columns * sizeof(double));
    double objval = (double)CPX_INFBOUND;

    if (CPXcallbackgetcandidatepoint(context, x_value, 0, num_x_columns - 1, &objval)) return 1;

    addcut(inst, x_value, context);
    free(x_value);
    return 0;
}

/**
 * @brief Add constraints on lazy callback
 * 
 * @param inst Instance struct
 * @param x_value Solution value
 * @param context Callback context structure
 * @return int Return 0
 */
int addcut(instance *inst, double *x_value, CPXCALLBACKCONTEXTptr context) {
    solution *sol = (solution *)malloc(sizeof(solution));
    sol->successor = (int *)malloc(inst->nnodes * sizeof(int));
    sol->component = (int *)malloc(inst->nnodes * sizeof(int));
    build_solution(x_value, inst, sol);

    if (sol->ncomp >= 2) {
        //initialize useful structure
        double *nnodes_component = (double *)calloc((sol->ncomp), sizeof(double));
        double rhs;
        int nzcnt;
        int beg = 0;
        char sense = LESS_EQUAL_CONSTRAINT_PARAMETER;

        //count number of nodes for each component
        for (int i = 0; i < inst->nnodes; i++) {
            nnodes_component[sol->component[i]] += 1;
        }

        for (int k = 0; k < sol->ncomp; k++) {
            rhs = nnodes_component[k] - 1;
            nzcnt = nnodes_component[k] * (nnodes_component[k] - 1) / 2;
            int *index = (int *)malloc(nzcnt * sizeof(int));
            double *value = (double *)malloc(nzcnt * sizeof(double));

            int count = 0;
            for (int i = 0; i < inst->nnodes; i++) {
                if (sol->component[i] == k) {
                    for (int j = i + 1; j < inst->nnodes; j++) {
                        if (sol->component[j] == k) {
                            index[count] = triangle_xpos(sol->successor[i], sol->successor[j], inst);
                            value[count] = 1.0;
                            count++;
                        }
                    }
                }
            }
            CPXcallbackrejectcandidate(context, 1, nzcnt, &rhs, &sense, &beg, index, value);
            free(index);
            free(value);
        }

        free(nnodes_component);
    }
    //free memory
    free_solution(sol);
    free(sol);

    return 0;
}

/**
 * @brief Build the model for cplex (Symmetric TSP)
 * 
 * @param inst Pointer to struct
 * @param env Environment
 * @param lp Problem variable
 */
void build_model(const instance *inst, CPXENVptr env, CPXLPptr lp) {
    char binary = CPX_BINARY;

    //required by cplex
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    // add binary var.s x(i,j) for i < j
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = dist(i, j, inst);  // cost == distance
            double lb = 0.0;                //values' lower bound
            double ub = 1.0;                //values' upper bound
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(ERROR_CREATING_NEWCOLUMN_STRING);
            //Check the position
            if (CPXgetnumcols(env, lp) - 1 != triangle_xpos(i, j, inst)) print_error(ERROR_CHECKING_POSITION_STRING);
        }
    }

    // add the degree constraints
    for (int h = 0; h < inst->nnodes; h++)  // degree constraints
    {
        int newrow = CPXgetnumrows(env, lp);
        double rhs = 2.0;
        char sense = EQUALITY_CONSTRAINT_PARAMETER;
        sprintf(cname[0], "degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_STRING);
        for (int i = 0; i < inst->nnodes; i++) {
            if (i == h) continue;
            if (CPXchgcoef(env, lp, newrow, triangle_xpos(i, h, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_STRING);
        }
    }

    if (VERBOSE >= 100) CPXwriteprob(env, lp, MODEL_NAME, NULL);

    free(cname[0]);
    free(cname);
}

/**
 * @brief Compute edge's position in the model (Symmetric TSP)
 * 
 * @param i First node
 * @param j Second node
 * @param inst Pointer to struct
 * @return int : Position in the model
 */
int triangle_xpos(int i, int j, const instance *inst) {
    if (i == j) print_error(" i == j in xpos");
    if (i > j) return triangle_xpos(j, i, inst);
    return i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
}

/**
 * @brief Build solution to verify how many subtour cplex found
 * 
 * @param x_value Value assigned to arcs from cplex
 * @param inst Instance struct
 * @param solution Struct to store solution
 */
void build_solution(const double *x_value, const instance *inst, solution *solution) {
#ifdef DEBUG
    //check the degree of each node in the graph
    int *degree = (int *)calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            int k = triangle_xpos(i, j, inst);
            //Chek if value is near 0 or 1: the expected solution
            if (fabs(x_value[k]) > EPSILON && fabs(x_value[k] - 1.0) > EPSILON) print_error("Wrong x_value in build_solution()");
            if (x_value[k] > 0.5) {
                ++degree[i];
                ++degree[j];
            }
        }
    }
    for (int i = 0; i < inst->nnodes; i++) {
        if (degree[i] != 2) print_error("Wrong degree in build_solution()");
    }
    free(degree);
#endif
    solution->ncomp = -1;
    for (int i = 0; i < inst->nnodes; i++) {
        solution->successor[i] = -1;
        solution->component[i] = -1;
    }

    for (int start = 0; start < inst->nnodes; start++) {
        if (solution->component[start] == -1) {
            solution->ncomp += 1;

            int i = start;
            int done = 0;
            //done = 1 when a tour/subtour ends
            while (!done) {
                done = 1;
                solution->component[i] = solution->ncomp;
                for (int j = 0; j < inst->nnodes; j++) {
                    if (i != j && x_value[triangle_xpos(i, j, inst)] > 0.5 && solution->component[j] == -1) {
                        //Found an edge available in the solution with a node not visited yet
                        solution->successor[i] = j;
                        i = j;  //update next node
                        solution->component[j] = solution->ncomp;
                        done = 0;
                        break;
                    }
                }
            }
            solution->successor[i] = start;
        }
    }
    solution->ncomp += 1;  //index component start from 0. ncomp is the number of component
}

/**
 * @brief Set callback mode for CPLEX
 * 
 * @param env CPLEX environment pointer
 * @param lp CPLEX problem pointer
 */
void init_env_callback(CPXENVptr env, CPXLPptr lp, instance *inst) {
    if (CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF)) print_error(ERROR_SETTING_MIPCREDLP_STRING);
    if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, mycallback, inst)) print_error(ERROR_SETTING_CALLBACK_STRING);
}

/**
 * @brief Initialize CPLEX environment
 * 
 * @param env CPLEX environment pointer
 * @param lp CPLEX problem pointer
 */
void init_cplex_env(CPXENVptr *env, CPXLPptr *lp, const instance *inst) {
    int error;
    *env = CPXopenCPLEX(&error);
    if (error) print_error(ERROR_OPEN_CPLEX_STRING);
    *lp = CPXcreateprob(*env, &error, PROBLEM_NAME_STRING);
    if (error) print_error(ERROR_CREATE_PROBLEM_STRING);

    if (VERBOSE >= 100)
        if (CPXsetintparam(*env, CPX_PARAM_SCRIND, CPX_ON)) print_error(ERROR_SETTING_SCRIND_STRING);
    if (inst->timelimit != -1.0)
        if (CPXsetdblparam(*env, CPX_PARAM_TILIM, inst->timelimit)) print_error(ERROR_SETTING_TIMELIM_STRING);
    if (inst->randomseed != -1)
        if (CPXsetintparam(*env, CPX_PARAM_RANDOMSEED, inst->randomseed)) print_error(ERROR_SETTING_RANDSEED_STRING);
}

/**
 * @brief Initialize instance for heuristic problem
 * 
 * @param inst Pointer to instance structure
 */
void init_heuristic_inst(instance *inst) {
    int num_nodes = inst->nnodes;
    int num_x_columns = (num_nodes * num_nodes) - num_nodes * (num_nodes + 1) / 2;
    inst->nodes_tour = (int *)malloc(num_nodes * sizeof(int));
    inst->distance = (double *)malloc(num_x_columns * sizeof(double));

    //init distance between nodes
    for (int i = 0; i < num_x_columns; i++) {
        inst->distance[i] = -1.0;
    }
}
