#include "utils.h"

#include <stdio.h>

#include "mystring.h"

/**
 * @brief Print error and close the program
 * 
 * @param err String to show
 */
void print_error(const char *err) {
    printf("\n\n ERROR: %s \n\n", err);
    fflush(NULL);
    exit(1);
}
/**
 * @brief Creates a gnuplot pipe and returns it
 * 
 * @return FILE* 
 */
FILE *open_init_plot() {
    FILE *gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
    char *commands[] = {"set title \"Graph\"", "plot 'points' with lines"};
    sendCommands(gnuplotPipe, commands, 1);
    fflush(gnuplotPipe);
    return gnuplotPipe;
}

/**
 * @brief Sends commands to the pipe
 * 
 * @param gnuplotPipe 
 * @param commandsForGnuplot 
 * @param numCommands 
 */

void sendCommands(FILE *gnuplotPipe, char *commandsForGnuplot[], int numCommands) {
    for (int i = 0; i < numCommands; i++)
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
}

/**
 * @brief Creates the file with the points to plot
 * 
 * @param inst 
 * @param sol 
 * @param temp 
 */
void createPointsFile(const instance *inst, const solution *sol, FILE *temp) {
    int visited_comp = 0;
    for (int start = 0; start < inst->nnodes; start++) {
        if (sol->component[start] == visited_comp) {
            int j = sol->successor[start];  //prendo successivo
            fprintf(temp, "%lf %lf \n", inst->xcoord[start], inst->ycoord[start]);
            fprintf(temp, "%lf %lf \n", inst->xcoord[j], inst->ycoord[j]);
            while (j != start) {
                j = sol->successor[j];
                fprintf(temp, "%lf %lf \n", inst->xcoord[j], inst->ycoord[j]);
            }
            ++visited_comp;
            fprintf(temp, "\n");
        }
    }
}

/**
 * @brief Creates the file with the points to plot from heuristic solution
 * 
 * @param inst Instance structure
 * @param temp Temporary file to write coord
 */
void createHeurPointsFile(const instance *inst, FILE *temp) {
    for (int i = 0; i < inst->nnodes; i++) {
        fprintf(temp, "%lf %lf \n", inst->xcoord[inst->nodes_tour[i]], inst->ycoord[inst->nodes_tour[i]]);
    }
    fprintf(temp, "%lf %lf \n", inst->xcoord[inst->nodes_tour[0]], inst->ycoord[inst->nodes_tour[0]]);
}

/**
 * @brief Plot points and egdes of a tour using GNUPlot 
 * 
 * @param inst Pointer to struct
 * @param sol Solution struct
 */
void plotData(const instance *inst, const solution *sol) {
    char *commandsForGnuplot[] = {"set title \"Graph\"", "plot 'points' with lines"};
    FILE *temp = fopen("points", "w");
    FILE *gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
    //int numCommands = 2;
    // for (int i = 0; i < inst->nnodes; i++) {
    //     fprintf(temp, "%lf %lf \n", inst->xcoord[i], inst->ycoord[i]);  //Write the data to a temporary file
    // }

    //Plot
    for (int i = 0; i < 1; i++) {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);  //Send commands to gnuplot one by one.
    }

    int visited_comp = 0;
    for (int start = 0; start < inst->nnodes; start++) {
        if (sol->component[start] == visited_comp) {
            int j = sol->successor[start];  //prendo successivo
            fprintf(temp, "%lf %lf \n", inst->xcoord[start], inst->ycoord[start]);
            fprintf(temp, "%lf %lf \n", inst->xcoord[j], inst->ycoord[j]);
            while (j != start) {
                j = sol->successor[j];
                fprintf(temp, "%lf %lf \n", inst->xcoord[j], inst->ycoord[j]);
            }
            ++visited_comp;
            fprintf(temp, "\n");
        }
    }

    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[1]);

    fclose(temp);
    pclose(gnuplotPipe);
}

/**
 * @brief Parsing of command line to get more info
 * 
 * @param argc 
 * @param argv 
 * @param inst Pointer to struct to update info
 */
void parse_command_line(int argc, char **argv, instance *inst) {
    if (VERBOSE >= 100) printf(" running %s with %d parameters \n", argv[0], argc - 1);

    // default
    strcpy(inst->input_file, "NULL");
    inst->timelimit = -1.0;
    inst->randomseed = -1;
    inst->model_type = 0;
    inst->initialization = 0;

    int help = 0;
    if (argc < 1) help = 1;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-file") == 0) {
            strcpy(inst->input_file, argv[++i]);
            continue;
        }  // input file
        if (strcmp(argv[i], "-input") == 0) {
            strcpy(inst->input_file, argv[++i]);
            continue;
        }  // input file
        if (strcmp(argv[i], "-f") == 0) {
            strcpy(inst->input_file, argv[++i]);
            continue;
        }  // input file
        if (strcmp(argv[i], "-time_limit") == 0) {
            inst->timelimit = atof(argv[++i]);
            continue;
        }  // total time limit
        if (strcmp(argv[i], "-randomseed") == 0) {
            inst->randomseed = atoi(argv[++i]);
            continue;
        }  // randomseed
        if (strcmp(argv[i], "-model_type") == 0) {
            inst->model_type = atoi(argv[++i]);
            continue;
        }  // model_type
        if (strcmp(argv[i], "-init") == 0) {
            int value = atoi(argv[++i]);
            if (value > 1) print_error(ERROR_INIT_VALUE_STRING);
            inst->initialization = value;
            continue;
        }
        if (strcmp(argv[i], "-help") == 0) {
            help = 1;
            continue;
        }  // help
        if (strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }  // help
        help = 1;
    }

    // print current parameters
    if (VERBOSE >= 10) {
        printf("\n----- Available parameters -----\n");
        printf("-file %s\n", inst->input_file);
        printf("-model_type %d\n", inst->model_type);
        if (inst->timelimit != -1.0) printf("-time_limit %lf\n", inst->timelimit);
        if (inst->randomseed != -1) printf("-randomseed %d\n", inst->randomseed);
        if (inst->initialization != -1) printf("-init %d\n", inst->initialization);
        printf("\nenter -help or --help for help\n");
        printf("--------------------------------\n\n");
    }

    if (help) {
        char models[11][50] = {
            "0: Tsp with loop",
            "1: Tsp with lazy callback",
            "2: Compact MTZ",
            "3: Compact Lazy MTZ",
            "4: Compact Flow 1",
            "5: Compact Lazy Flow1"
            "6: Hard Fixing",
            "7: Soft Fixing",
            "8: 2-Opt",
            "9: Variable Neighborhood Search",
            "10: Simulated annealing"};
        printf("Parameters:\n");
        printf("-file | -f | -input : input file\n");
        printf("-time_limit : time limit in seconds.\n(N.B.VNS : number of iterations)\n");
        printf("-model_type : how to solve problem:\n");
        for (int i = 0; i < 11; i++) {
            printf("%s\n", models[i]);
        }
        printf("-init : initialization type for heuristic solver without cplex. 0 for NN, 1 for GRASP\n");

        exit(1);
    }
}

/**
 * @brief Read input file
 * 
 * @param inst Pointer to struct
 */
void read_input(instance *inst)  // simplified CVRP parser, not all SECTIONs detected
{
    FILE *fin = fopen(inst->input_file, "r");
    if (fin == NULL) print_error(FILE_NOT_FOUND_STRING);
    char line[180];
    char *par_name;
    char *token1;
    char *token2;

    inst->nnodes = -1;
    inst->nodes_tour = NULL;
    inst->distance = NULL;

    int active_section = 0;  // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION

    int do_print = (VERBOSE >= 1000);

    while (fgets(line, sizeof(line), fin) != NULL) {
        if (VERBOSE >= 2000) {
            printf("%s", line);
            fflush(NULL);
        }
        if (strlen(line) <= 1) continue;  // skip empty lines
        par_name = strtok(line, " :");
        if (VERBOSE >= 3000) {
            printf("parameter \"%s\" ", par_name);
            fflush(NULL);
        }

        if (strncmp(par_name, "NAME", 4) == 0) {
            active_section = 0;
            continue;
        }

        if (strncmp(par_name, "COMMENT", 7) == 0) {
            active_section = 0;
            token1 = strtok(NULL, "");
            //if (VERBOSE >= 10) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
            continue;
        }

        if (strncmp(par_name, "TYPE", 4) == 0) {
            token1 = strtok(NULL, " :");
            //if (strncmp(token1, "CVRP", 4) != 0) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!");
            active_section = 0;
            continue;
        }

        if (strncmp(par_name, "DIMENSION", 9) == 0) {
            if (inst->nnodes >= 0) print_error(REPEATED_DIMENSION_STRING);
            token1 = strtok(NULL, " :");
            inst->nnodes = atoi(token1);
            if (do_print) printf(" ... nnodes %d\n", inst->nnodes);
            inst->xcoord = (double *)calloc(inst->nnodes, sizeof(double));
            inst->ycoord = (double *)calloc(inst->nnodes, sizeof(double));
            active_section = 0;
            continue;
        }

        if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0) {
            token1 = strtok(NULL, " :");
            if (strncmp(token1, "EUC_2D", 6) != 0) print_error(EDGE_WEIGHT_TYPE_ERROR_STRING);
            active_section = 0;
            continue;
        }

        if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0) {
            if (inst->nnodes <= 0) print_error(DIMENSION_ERROR_STRING);
            active_section = 1;
            continue;
        }

        if (strncmp(par_name, "EOF", 3) == 0) {
            active_section = 0;
            break;
        }

        if (active_section == 1)  // within NODE_COORD_SECTION
        {
            int i = atoi(par_name) - 1;
            if (i < 0 || i >= inst->nnodes) print_error(NODE_ERROR_STRING);
            token1 = strtok(NULL, " :,");
            token2 = strtok(NULL, " :,");
            inst->xcoord[i] = atof(token1);
            inst->ycoord[i] = atof(token2);
            if (do_print) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i + 1, inst->xcoord[i], inst->ycoord[i]);
            continue;
        }

        print_error(WRONG_FORMAT_ERROR);
    }

    fclose(fin);
}