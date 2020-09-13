#ifndef MY_STRING
#define MY_STRING

#include <string.h>

#define FILE_NOT_FOUND_STRING "Input file not found!"
#define REPEATED_DIMENSION_STRING "Repeated DIMENSION section in input file"
#define EDGE_WEIGHT_TYPE_ERROR_STRING "Format error: only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!"
#define DIMENSION_ERROR_STRING "... DIMENSION section should appear before NODE_COORD_SECTION section"
#define NODE_ERROR_STRING "...Unknown node in NODE_COORD_SECTION section"
#define WRONG_FORMAT_ERROR "...Wrong format for the current simplified parser!"
#define TIME_LIMIT_NOT_FOUND_STRING "Time limit not defined"
#define ERROR_INIT_VALUE_STRING "Wrong initialization value"

//Error
#define ERROR_OPEN_CPLEX_STRING "Error opening cplex variable"
#define ERROR_CREATE_PROBLEM_STRING "Error creating TSP problem variable"
#define ERROR_GET_X_STRING "Error getting x"

#define ERROR_SETTING_SCRIND_STRING "Error getting on CPX_PARAM_SCRIND"
#define ERROR_SETTING_MIPCREDLP_STRING "Error getting off CPX_PARAM_MIPCBREDLP"
#define ERROR_SETTING_CALLBACK_STRING "Error setting callbacksetfunc"
#define ERROR_SETTING_TIMELIM_STRING "Error setting time limit"
#define ERROR_SETTING_RANDSEED_STRING "Error setting randomseed"

#define ERROR_CREATING_NEWCOLUMN_STRING "wrong CPXnewcols on x var.s"
#define ERROR_CHECKING_POSITION_STRING "wrong position for x var.s"
#define ERROR_CREATING_NEWROW_STRING "wrong CPXnewrows [degree]"
#define ERROR_MODIFYING_COEF_STRING "wrong CPXchgcoef [degree]"

#define ERROR_CREATING_NEWROW_DEGREE_IN_STRING "Wrong CPXnewrows Degree_In"
#define ERROR_CREATING_NEWROW_DEGREE_OUT_STRING "Wrong CPXnewrows Degree_Out"
#define ERROR_MODIFYING_COEF_DEGREE_IN_STRING "Wrong CPXchgcoef Degree_In"
#define ERROR_MODIFYING_COEF_DEGREE_OUT_STRING "Wrong CPXchgcoef Degree_Out"

#define ERROR_CREATING_NEWROW_MTZ_STRING "wrong CPXnewrows MTZ_Constraint"
#define ERROR_CREATING_NEWCOLUMN_MTZ_STRING "Wrong CPXnewcols on MTZ u variable"
#define ERROR_CHECKING_POSITION_MTZ_STRING "Wrong position for MTZ u variable"
#define ERROR_MODIFYING_COEF_MTZ_STRING "Wrong CPXchgcoef MTZ_Constraint"
#define ERROR_ADDING_NEWLAZY_MTZ_STRING "Wrong adding lazy MTZ_Constraint"

#define ERROR_CREATING_NEWROW_FLOW1_STRING "Wrong CPXnewrows FLOW1_Constraint"
#define ERROR_CREATING_NEWCOLUMN_FLOW1_STRING "Wrong CPXnewcols on FLOW1 y variable"
#define ERROR_CHECKING_POSITION_FLOW1_STRING "Wrong position for FLOW1 y variable"
#define ERROR_MODIFYING_COEF_X_FLOW1_STRING "Wrong CPXchgcoef FLOW1_Constraint"
#define ERROR_MODIFYING_COEF_Y_FLOW1_STRING "Wrong CPXchgcoef FLOW1_Constraint"
#define ERROR_ADDING_NEWLAZY_FLOW1_STRING "Wrong adding lazy FLOW1_Constraint"
#define ERROR_ADDING_NEWLAZY_2NODE_SEL_STRING "Wrong CPXlazyconstraints on 2-node SECs"

#define ERROR_LOWER_BOUND_FIX_VAR_STRING "Error setting lower_bound in fix_variable"
#define ERROR_OBJVALUE_STRING "Error getting solution value"
#define ERROR_BEST_OBJVALUE_STRING "Error getting Best object value"
#define ERROR_DEL_ROW_STRING "Error deleting the last row"
#define ERROR_DEFAULT_LB_STRING "Error setting default LB"
#define ERROR_DEFAULT_UP_STRING "Error setting default UB"

#define ERROR_PLOT_VALUE_STRING "Wrong plot value"
#define ERROR_DISTANCE_STRING "Distance not computed!"

#define NN_INIT_STRING "Nearest Neighborhood initialization\n"
#define GRASP_INIT_STRING "Random GRASP initialization\n"

//TSP Utils
#define PROBLEM_NAME_STRING "TSP"
#define MODEL_NAME "model.lp"

// #define BINARY_VALUE_PARAMETER 'B'
// #define INTEGER_VALUE_PARAMETER 'I'
#define EQUALITY_CONSTRAINT_PARAMETER 'E'
#define LESS_EQUAL_CONSTRAINT_PARAMETER 'L'
#define GREATER_EQUAL_CONSTRAINT_PARAMETER 'G'

#define LOWER_BOUND_CHAR 'L'
#define UPPER_BOUND_CHAR 'U'

// #define EPSILON 1e-5

//GNUPlot
#define GNUPLOT_OPEN_COMMAND "gnuplot -persistent"

//VERBOSE
#define VERBOSE_LITTLE_OUTPUT 20
#define VERBOSE_TO_PLOT 80
#endif /* MY_STRING */
