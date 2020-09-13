#include "compact_models_TSP.h"

#include "mystring.h"
#include "utils.h"

/**
 * @brief Compute the optimal solution using cplex
 * 
 * @param inst Instance struct
 * @param sol Solution struct
 * @param model_type 0 for MTZ, 1 for Flow1
 * @param enable_lazy 1 to enable lazy constraints
 * @return int 0 if successful and nonzero if an error occurs
 */
int compact_TSPopt(const instance *inst, const int model_type, const int enable_lazy) {
    if (model_type > 1 || enable_lazy > 1)
        return -1;

    time_t start, end;
    int time_taken;

    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    init_cplex_env(&env, &lp, inst);
    atsp_init_model(inst, env, lp);

    //Choose model type
    if (model_type == 0) {  //mtz model
        if (enable_lazy == 1) {
            printf("lazy_MTZ_build_model\n");
            lazy_MTZ_build_model(inst, env, lp);
        } else {
            printf("MTZ_build_model\n");
            MTZ_build_model(inst, env, lp);
        }
    } else {  //flow1 model
        if (enable_lazy == 1) {
            printf("lazy_flow1_build_model\n");
            lazy_FLOW1_build_model(inst, env, lp);
        } else {
            printf("flow1_build_model\n");
            FLOW1_build_model(inst, env, lp);
        }
    }

    //save model
    if (VERBOSE >= 100) CPXwriteprob(env, lp, MODEL_NAME, NULL);

    time(&start);
    // Solve the problem
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
            int num_nodes = inst->nnodes;
            int num_x_columns = (inst->nnodes) * (inst->nnodes);
            //get x solution
            double *x_value = (double *)calloc(num_x_columns, sizeof(double));
            int getx_code = CPXgetx(env, lp, x_value, 0, num_x_columns - 1);
            if (getx_code) print_error(ERROR_GET_X_STRING);

            solution *sol = (solution *)malloc(sizeof(solution));
            //malloc solution struct
            sol->successor = (int *)malloc(num_nodes * sizeof(int));
            sol->component = (int *)malloc(num_nodes * sizeof(int));

            //build solution to plot
            build_atsp_solution(x_value, inst, sol);

            gnuplotpipe = open_init_plot();
            temp = fopen("points", "w");
            createPointsFile(inst, sol, temp);
            char *plot[] = {"plot 'points' with lines"};
            sendCommands(gnuplotpipe, plot, 1);
            fclose(temp);
            fflush(gnuplotpipe);
            pclose(gnuplotpipe);

            //free pointer
            free(x_value);
            free_solution(sol);
            free(sol);
        }
    }

    //free and close cplex model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    return 0;
}

/**
 * @brief Build model for asymmetric version
 * 
 * @param inst Instance struct
 * @param env Environment
 * @param lp Problem variable
 */
void atsp_init_model(const instance *inst, CPXENVptr env, CPXLPptr lp) {
    char binary = CPX_BINARY;
    double lb = 0.0;
    //required by cplex
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = dist(i, j, inst);
            double ub = (i == j) ? 0.0 : 1.0;
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(ERROR_CREATING_NEWCOLUMN_STRING);
            if (CPXgetnumcols(env, lp) - 1 != matrix_xpos(i, j, inst)) print_error(ERROR_CHECKING_POSITION_STRING);
        }
    }

    // add the degree constraints
    double rhs = 1.0;
    char sense = EQUALITY_CONSTRAINT_PARAMETER;

    for (int h = 0; h < inst->nnodes; h++) {  // degree constraints
        int newrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "degree_Out(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_DEGREE_OUT_STRING);
        sprintf(cname[0], "degree_In(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_DEGREE_IN_STRING);
        for (int i = 0; i < inst->nnodes; i++) {
            if (i == h) continue;
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(i, h, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_DEGREE_IN_STRING);
            if (CPXchgcoef(env, lp, newrow + 1, matrix_xpos(h, i, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_DEGREE_OUT_STRING);
        }
    }

    free(cname[0]);
    free(cname);
}

/**
 * @brief Build atsp solution to verify how many subtour CPLEX found
 * 
 * @param x_value Value assigned to arcs from CPLEX
 * @param inst Instance struct
 * @param solution Solution struct
 */
void build_atsp_solution(const double *x_value, const instance *inst, solution *solution) {
    //Init solution struct value
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
                    if (i != j && x_value[matrix_xpos(i, j, inst)] > 0.5 && solution->component[j] == -1) {
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
 * @brief Add "Miller, Tucker and Zemlin" constraints to the model
 * 
 * @param inst Instance struct
 * @param env Environment
 * @param lp Problem variable
 */
void MTZ_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp) {
    //required by cplex
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    //add "Miller, Tucker and Zemlin" Constraints
    int last_xvalue_pos = CPXgetnumcols(env, lp) - 1;
    char continuous = CPX_CONTINUOUS;
    double lb = 0.0;
    double ub = 0.0;
    lb = 0;
    ub = inst->nnodes - 2;
    for (int i = 1; i < inst->nnodes; i++) {
        sprintf(cname[0], "u(%d)", i + 1);
        if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &continuous, cname)) print_error(ERROR_CREATING_NEWCOLUMN_MTZ_STRING);
        if (CPXgetnumcols(env, lp) - 1 != (last_xvalue_pos + i)) print_error(ERROR_CHECKING_POSITION_MTZ_STRING);
    }

    double rhs = inst->nnodes - 1;
    char sense = LESS_EQUAL_CONSTRAINT_PARAMETER;
    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = 1; j < inst->nnodes; j++) {
            if (i == j) continue;
            int newrow = CPXgetnumrows(env, lp);
            sprintf(cname[0], "MTZ_Constraint(%d,%d)", i + 1, j + 1);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_MTZ_STRING);
            //update coef
            if (CPXchgcoef(env, lp, newrow, last_xvalue_pos + i, 1.0)) print_error(ERROR_MODIFYING_COEF_MTZ_STRING);
            if (CPXchgcoef(env, lp, newrow, last_xvalue_pos + j, -1.0)) print_error(ERROR_MODIFYING_COEF_MTZ_STRING);
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(i, j, inst), inst->nnodes)) print_error(ERROR_MODIFYING_COEF_MTZ_STRING);
        }
    }

    free(cname[0]);
    free(cname);
}

/**
 * @brief Add "Miller, Tucker and Zemlin" constraint to the model as lazy constraing
 * 
 * @param inst Instance struct
 * @param env Environment
 * @param lp Problem variable
 */
void lazy_MTZ_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp) {
    //required by cplex
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    //add "Miller, Tucker and Zemlin" Constraints
    int last_xvalue_pos = CPXgetnumcols(env, lp) - 1;
    char continuous = CPX_CONTINUOUS;
    double lb = 0.0;
    double ub = 0.0;
    lb = 0;
    ub = inst->nnodes - 2;
    for (int i = 1; i < inst->nnodes; i++) {
        sprintf(cname[0], "u(%d)", i + 1);
        if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &continuous, cname)) print_error(ERROR_CREATING_NEWCOLUMN_MTZ_STRING);
        if (CPXgetnumcols(env, lp) - 1 != (last_xvalue_pos + i)) print_error(ERROR_CHECKING_POSITION_MTZ_STRING);
    }

    double rhs = inst->nnodes - 1;
    char sense = LESS_EQUAL_CONSTRAINT_PARAMETER;
    int nnz = 3;
    int index[3];
    int beg = 0;
    double value[3] = {+1.0, -1.0, inst->nnodes};
    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = 1; j < inst->nnodes; j++) {
            if (i == j) continue;
            sprintf(cname[0], "MTZ_constraints(%d,%d)", i + 1, j + 1);
            //update index
            index[0] = last_xvalue_pos + i;
            index[1] = last_xvalue_pos + j;
            index[2] = matrix_xpos(i, j, inst);
            //add lazy constraint
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &beg, index, value, cname)) print_error(ERROR_ADDING_NEWLAZY_MTZ_STRING);
        }
    }

    //add more lazy constraints: x_ij + x_ji <=1
    rhs = 1.0;
    nnz = 2;
    value[0] = value[1] = 1.0;
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            sprintf(cname[0], "SEC on node pair (%d,%d)", i + 1, j + 1);
            index[0] = matrix_xpos(i, j, inst);
            index[1] = matrix_xpos(j, i, inst);
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &beg, index, value, cname)) print_error(ERROR_ADDING_NEWLAZY_2NODE_SEL_STRING);
        }
    }

    free(cname[0]);
    free(cname);
}

/**
 * @brief Build model for flow1
 * 
 * @param inst Pointer to struct
 * @param env Environment
 * @param lp Problem variable
 */
void FLOW1_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp) {
    char continuos = CPX_CONTINUOUS;
    double lb = 0.0;
    double ub;
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    //add y_ij variables
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
            ub = ((i == j) || (j == 0)) ? 0.0 : CPX_INFBOUND;
            if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &continuos, cname)) print_error(ERROR_CREATING_NEWCOLUMN_FLOW1_STRING);
            if (CPXgetnumcols(env, lp) - 1 != matrix_xpos(inst->nnodes + i, j, inst)) print_error(ERROR_CHECKING_POSITION_FLOW1_STRING);
        }
    }

    double rhs = inst->nnodes - 1;
    int newrow = CPXgetnumrows(env, lp);
    char sense = EQUALITY_CONSTRAINT_PARAMETER;
    sprintf(cname[0], "y_degree(%d)", 1);
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);

    //add ys constraints for node 1
    for (int j = 0; j < inst->nnodes; j++) {
        if (j == 0) continue;
        if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
    }

    //add sum y_hj = sum y_ih - 1 for each h in V\{1} constraints
    rhs = -1.0;
    for (int h = 1; h < inst->nnodes; h++) {
        int newrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "y_degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);
        //set coefficient for outgoing edges
        for (int j = 0; j < inst->nnodes; j++) {
            if (j == h) continue;
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + h, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
        }
        //set coefficient for incoming edges
        for (int i = 0; i < inst->nnodes; i++) {
            if (i == h) continue;
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + i, h, inst), -1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
        }
    }

    sense = LESS_EQUAL_CONSTRAINT_PARAMETER;
    rhs = 0.0;
    int i = 0;

    //add y_ij <= (n - 1) * x_ij constraints for i = 1, j != 1
    for (int j = 0; j < inst->nnodes; j++) {
        if (i == j) continue;
        int newrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "y_set1(%d,%d)", i + 1, j + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);
        //set coefficient for x(1,j)
        if (CPXchgcoef(env, lp, newrow, matrix_xpos(i, j, inst), -(inst->nnodes - 1))) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
        //set coefficient for y(1,j)
        if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + i, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
    }

    //add y_ij <= (n - 2) * x_ij constraints for y(i,j) i,j != 1
    for (i = 1; i < inst->nnodes; i++) {
        for (int j = 1; j < inst->nnodes; j++) {
            if (i == j) continue;
            int newrow = CPXgetnumrows(env, lp);
            sprintf(cname[0], "y_set2(%d,%d)", i + 1, j + 1);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);
            //set coefficient for x(1,j)
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(i, j, inst), -(inst->nnodes - 2))) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
            //set coefficient for y(1,j)
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + i, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
        }
    }

    free(cname[0]);
    free(cname);
}

/**
 * @brief Build model for flow 1 with lazy constraints
 * 
 * @param inst Pointer to struct
 * @param env Environment
 * @param lp Problem variable
 */

void lazy_FLOW1_build_model(const instance *inst, CPXENVptr env, CPXLPptr lp) {
    char continuous = CPX_CONTINUOUS;
    double lb = 0.0;
    double ub;
    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    //add y_ij variables
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            sprintf(cname[0], "y(%d, %d)", i + 1, j + 1);
            ub = ((i == j) || (j == 0)) ? 0.0 : CPX_INFBOUND;
            if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &continuous, cname)) print_error(ERROR_CREATING_NEWCOLUMN_FLOW1_STRING);
            if (CPXgetnumcols(env, lp) - 1 != matrix_xpos(inst->nnodes + i, j, inst)) print_error(ERROR_CHECKING_POSITION_FLOW1_STRING);
        }
    }

    double rhs = inst->nnodes - 1;
    int newrow = CPXgetnumrows(env, lp);
    char sense = EQUALITY_CONSTRAINT_PARAMETER;
    sprintf(cname[0], "y_degree(%d)", 1);
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);

    //add ys constraints for node 1
    for (int j = 0; j < inst->nnodes; j++) {
        if (j == 0) continue;
        if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
    }

    //add sum y_hj = sum y_ih - 1 constraints for each h in V\{1}
    rhs = -1.0;
    for (int h = 1; h < inst->nnodes; h++) {
        int newrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "y_degree(%d)", h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);
        //set coefficient for outgoing edges
        for (int j = 0; j < inst->nnodes; j++) {
            if (j == h) continue;
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + h, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
        }
        //set coefficient for incoming edges
        for (int i = 0; i < inst->nnodes; i++) {
            if (i == h) continue;
            if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + i, h, inst), -1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
        }
    }

    sense = LESS_EQUAL_CONSTRAINT_PARAMETER;
    rhs = 0.0;
    int i = 0;
    int izero = 0;

    //add y_ij <= (n - 1) * x_ij constraints for i = 1, j != 1
    for (int j = 0; j < inst->nnodes; j++) {
        if (i == j) continue;
        int newrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "y_set1(%d,%d)", i + 1, j + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(ERROR_CREATING_NEWROW_FLOW1_STRING);
        //set coefficient for x_1j
        if (CPXchgcoef(env, lp, newrow, matrix_xpos(i, j, inst), -(inst->nnodes - 1))) print_error(ERROR_MODIFYING_COEF_X_FLOW1_STRING);
        //set coefficient for y_1j
        if (CPXchgcoef(env, lp, newrow, matrix_xpos(inst->nnodes + i, j, inst), 1.0)) print_error(ERROR_MODIFYING_COEF_Y_FLOW1_STRING);
    }

    //add y_ij <= (n - 2) * x_ij lazy constraints for y(i,j) i,j != 1
    int index[2];
    double value[2];
    int nnz = 2;
    rhs = 0.0;

    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = 1; j < inst->nnodes; j++) {
            if (i == j) continue;
            //int newrow = CPXgetnumrows(env, lp);
            sprintf(cname[0], "y_set2(%d,%d)", i + 1, j + 1);
            index[0] = matrix_xpos(i, j, inst);
            value[0] = -(inst->nnodes - 2);
            index[1] = matrix_xpos(inst->nnodes + i, j, inst);
            value[1] = 1.0;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) print_error(ERROR_ADDING_NEWLAZY_FLOW1_STRING);
        }
    }

    //add 1.0 * x_ij + 1.0 * x_ji <= 1 lazy constraints for each arc (i,j) with i < j
    rhs = 1.0;
    nnz = 2;
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            sprintf(cname[0], "SEL on node pair (%d, %d)", i + 1, j + 1);
            index[0] = matrix_xpos(i, j, inst);
            value[0] = 1.0;
            index[1] = matrix_xpos(j, i, inst);
            value[1] = 1.0;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) print_error(ERROR_ADDING_NEWLAZY_2NODE_SEL_STRING);
        }
    }

    free(cname[0]);
    free(cname);
}

/**
 * @brief Compute edge's poition in the model
 * 
 * @param i First node
 * @param j Second node
 * @param inst Instance struct
 * @return int Position in the model
 */
int matrix_xpos(int i, int j, const instance *inst) {
    return i * inst->nnodes + j;
}
