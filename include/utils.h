#include "tsp.h"

FILE *open_init_plot();
void sendCommands(FILE *gnuplotPipe, char *commandsForGnuplot[], int numCommands);
void createPointsFile(const instance *inst, const solution *sol, FILE *temp);
void createHeurPointsFile(const instance *inst, FILE *temp);
void closePlot(FILE *temp);
void plotData(const instance *inst, const solution *sol);
void parse_command_line(int argc, char **argv, instance *inst);
void print_error(const char *err);
void read_input(instance *inst);
