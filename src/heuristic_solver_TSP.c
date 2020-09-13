#include "heuristic_solver_TSP.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "mystring.h"
#include "utils.h"

/**
 * @brief Initialize solution with Nearest neighborhood
 * 
 * @param inst Instance structure
 */
void nearest_neighborhood(instance *inst) {
    int num_nodes = inst->nnodes;
    double best_obj_val = DBL_MAX;
    double tmp_obj_val = 0.0;
    int *tmp_nodes = (int *)malloc(num_nodes * sizeof(int));

    for (int start_node = 0; start_node < num_nodes; start_node++) {
        //init nodes
        for (int i = 0; i < num_nodes; i++) {
            tmp_nodes[i] = i;
        }

        tmp_obj_val = nearest_neighborhood_from_node(inst, start_node, tmp_nodes);

        //check for best solution
        if (tmp_obj_val < best_obj_val) {
            //update best solution
            best_obj_val = tmp_obj_val;
            for (int i = 0; i < num_nodes; i++) {
                inst->nodes_tour[i] = tmp_nodes[i];
            }
        }
    }

    printf("NN solution cost: %f\n", best_obj_val);

    // plot solution
    FILE *gnuplotPipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) {
        gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
        char *commands[] = {"set title \"Nearest Neighborhood\""};
        sendCommands(gnuplotPipe, commands, 1);
        //temporary file
        temp = fopen("points", "w");
        createHeurPointsFile(inst, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotPipe, plot, 1);
        fclose(temp);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }

    free(tmp_nodes);
}

/**
 * @brief Initialize solution with GRASP
 * 
 * @param inst Instance structure
 */
void random_GRASP(instance *inst) {
    if (inst->randomseed == -1)
        srand(time(NULL));  //random seed
    else
        srand(inst->randomseed);

    int num_nodes = inst->nnodes;
    double best_obj_val = DBL_MAX;
    double tmp_obj_val = 0.0;
    int *tmp_nodes = (int *)malloc(num_nodes * sizeof(int));

    int *index_neighborhood = (int *)malloc(3 * sizeof(int));

    //# of iterations?
    for (int iteration = 0; iteration < 100; iteration++) {
        //init index nodes
        for (int i = 0; i < num_nodes; i++) {
            tmp_nodes[i] = i;
        }

        int start_node = rand() % num_nodes;

        //random GRASP from start_node
        tmp_obj_val = random_GRASP_from_node(inst, start_node, index_neighborhood, 3, tmp_nodes);

        //check for best solution
        if (tmp_obj_val < best_obj_val) {
            //update best solution
            best_obj_val = tmp_obj_val;
            for (int i = 0; i < num_nodes; i++) {
                inst->nodes_tour[i] = tmp_nodes[i];
            }
        }
    }
    printf("GRASP solution cost: %f\n", best_obj_val);

    // plot solution
    FILE *gnuplotPipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) {
        gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
        char *commands[] = {"set title \"RandomGRASP\""};
        sendCommands(gnuplotPipe, commands, 1);
        //temporary file
        temp = fopen("points", "w");
        createHeurPointsFile(inst, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotPipe, plot, 1);
        fclose(temp);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }

    free(index_neighborhood);
    free(tmp_nodes);
}

/**
 * @brief Nearest Neighborhood starting from start_node
 * 
 * @param inst Instance structure
 * @param start_node Starting node
 * @param nodes Output route
 * @return double Route cost
 */
double nearest_neighborhood_from_node(instance *inst, int start_node, int *nodes) {
    //swap for index 0
    nodes[0] = start_node;
    nodes[start_node] = 0;
    double tmp_obj_val = 0;  //init tmp obj value
    //nearest neighborhood starting from start
    for (int i = 0; i < inst->nnodes - 1; i++) {
        int index_nearest_node = -1;
        double best_distance = DBL_MAX;
        for (int j = i + 1; j < inst->nnodes; j++) {
            //find nearest node from i
            int dist_position = triangle_xpos(nodes[i], nodes[j], inst);
            //compute distance
            if (inst->distance[dist_position] == -1.0)
                inst->distance[dist_position] = dist(nodes[i], nodes[j], inst);
            //update best node
            if (inst->distance[dist_position] < best_distance) {
                index_nearest_node = j;
                best_distance = inst->distance[dist_position];
            }
        }
        tmp_obj_val += best_distance;
        //swap value
        int tmp = nodes[i + 1];
        nodes[i + 1] = nodes[index_nearest_node];
        nodes[index_nearest_node] = tmp;
    }
    //update last edge
    tmp_obj_val += inst->distance[triangle_xpos(nodes[0], nodes[inst->nnodes - 1], inst)];

    return tmp_obj_val;
}

/**
 * @brief Random GRASP starting from start_node
 * 
 * @param inst Instance structure
 * @param start_node Starting node
 * @param index_neighborhood Array to store closest nodes 
 * @param index_neigh_size Size of index_neighborhood
 * @param nodes Output route
 * @return double Route cost
 */
double random_GRASP_from_node(instance *inst, int start_node, int *index_neighborhood, int index_neigh_size, int *nodes) {
    //swap for index 0
    nodes[0] = start_node;
    nodes[start_node] = 0;
    int neighborhood_size;
    double tmp_obj_val = 0;  //init tmp obj value
    for (int i = 0; i < inst->nnodes - 1; i++) {
        neighborhood_size = 0;
        double worst_distance = DBL_MAX;
        //find nearest node from i
        for (int j = i + 1; j < inst->nnodes; j++) {
            int dist_position = triangle_xpos(nodes[i], nodes[j], inst);
            //compute distance
            if (inst->distance[dist_position] == -1.0)
                inst->distance[dist_position] = dist(nodes[i], nodes[j], inst);

            //update best nodes
            //neighborhood is empty yet
            if (neighborhood_size < index_neigh_size) {
                int done = 0;
                for (int k = 0; k < neighborhood_size && !done; k++) {
                    double local_distance = inst->distance[triangle_xpos(nodes[i], nodes[index_neighborhood[k]], inst)];
                    if (inst->distance[dist_position] < local_distance) {
                        //position found
                        done = 1;
                        for (int h = neighborhood_size - 1; h >= k; h--) {
                            index_neighborhood[h + 1] = index_neighborhood[h];
                        }
                        index_neighborhood[k] = j;
                    }
                }
                if (done == 0) {
                    index_neighborhood[neighborhood_size] = j;       //insert in the last position
                    worst_distance = inst->distance[dist_position];  //update worst distance
                }
                neighborhood_size++;
            } else {
                //neighborhood contains index_neigh_size elements
                if (inst->distance[dist_position] < worst_distance) {
                    //new neighborhood to insert
                    int done = 0;
                    for (int k = 0; k < neighborhood_size && !done; k++) {
                        double local_distance = inst->distance[triangle_xpos(nodes[i], nodes[index_neighborhood[k]], inst)];
                        if (inst->distance[dist_position] < local_distance) {
                            //position found
                            done = 1;
                            for (int h = neighborhood_size - 2; h >= k; h--) {
                                index_neighborhood[h + 1] = index_neighborhood[h];
                            }
                            index_neighborhood[k] = j;
                            //update worst distance
                            worst_distance = inst->distance[triangle_xpos(nodes[i], nodes[index_neighborhood[neighborhood_size - 1]], inst)];
                        }
                    }
                }
            }
        }

        //0.90 prob greedy, 0.10 prob random
        int winner;
        double number = rand() / ((double)RAND_MAX);
        if (number > 0.10)
            winner = 0;
        else {
            //pick one node at random from neighborhood
            winner = rand() % neighborhood_size;
        }

        //update solution
        tmp_obj_val += inst->distance[triangle_xpos(nodes[i], nodes[index_neighborhood[winner]], inst)];
        //swap value
        int tmp = nodes[i + 1];
        nodes[i + 1] = nodes[index_neighborhood[winner]];
        nodes[index_neighborhood[winner]] = tmp;
    }
    //update last edge
    tmp_obj_val += inst->distance[triangle_xpos(nodes[0], nodes[inst->nnodes - 1], inst)];

    return tmp_obj_val;
}

/**
 * @brief 2-Opt optimization of a tour
 * 
 * @param inst Instance Structure
 * @param tour Tour to optimize
 * @param plot 1 to plot, 0 otherwise (VERBOSE needed too)
 */
void opt_2(instance *inst, int *tour, int plot) {
    if (plot > 1) print_error(ERROR_PLOT_VALUE_STRING);
    int num_nodes = inst->nnodes;
    int swap = 1;
    double c1, c2, c3, c4, delta;

    while (swap != 0) {
        swap = 0;
        double best_delta = 0.0;
        int best_swap_1, best_swap_2;
        for (int i = 1; i < num_nodes - 1; i++) {
            for (int k = i + 1; k < num_nodes; k++) {
                //check for distances
                if (inst->distance[triangle_xpos(tour[i - 1], tour[i], inst)] == -1)
                    inst->distance[triangle_xpos(tour[i - 1], tour[i], inst)] = dist(tour[i - 1], tour[i], inst);
                if (inst->distance[triangle_xpos(tour[k], tour[(k + 1) % num_nodes], inst)] == -1)
                    inst->distance[triangle_xpos(tour[k], tour[(k + 1) % num_nodes], inst)] = dist(tour[k], tour[(k + 1) % num_nodes], inst);
                if (inst->distance[triangle_xpos(tour[i - 1], tour[k], inst)] == -1)
                    inst->distance[triangle_xpos(tour[i - 1], tour[k], inst)] = dist(tour[i - 1], tour[k], inst);
                if (inst->distance[triangle_xpos(tour[i], tour[(k + 1) % num_nodes], inst)] == -1)
                    inst->distance[triangle_xpos(tour[i], tour[(k + 1) % num_nodes], inst)] = dist(tour[i], tour[(k + 1) % num_nodes], inst);

                c1 = inst->distance[triangle_xpos(tour[i - 1], tour[i], inst)];                //dist(tour[i - 1], tour[i], inst);
                c2 = inst->distance[triangle_xpos(tour[k], tour[(k + 1) % num_nodes], inst)];  //dist(tour[k], tour[(k + 1) % num_nodes], inst);
                c3 = inst->distance[triangle_xpos(tour[i - 1], tour[k], inst)];                //dist(tour[i - 1], tour[k], inst);
                c4 = inst->distance[triangle_xpos(tour[i], tour[(k + 1) % num_nodes], inst)];  //dist(tour[i], tour[(k + 1) % num_nodes], inst);
                if (c1 == -1 || c2 == -1 || c3 == -1 || c4 == -1) print_error(ERROR_DISTANCE_STRING);
                delta = (c3 + c4) - (c1 + c2);
                if (delta < best_delta) {
                    //update best swap data
                    best_delta = delta;
                    best_swap_1 = i;
                    best_swap_2 = k;
                    swap = 1;
                }
            }
        }
        if (swap == 1) {
            opt_2_swap(tour, best_swap_1, best_swap_2);
        }
    }

    double distance = 0.0;
    for (int i = 0; i < num_nodes; i++) {
        distance += dist(tour[i], tour[(i + 1) % num_nodes], inst);
    }
    printf("2-Opt Solution: %f\n", distance);

    // plot solution
    FILE *gnuplotPipe = NULL;
    FILE *temp = NULL;
    if (plot == 1 && VERBOSE >= VERBOSE_TO_PLOT) {
        gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
        char *commands[] = {"set title \"2-Opt\""};
        sendCommands(gnuplotPipe, commands, 1);
        //temporary file
        temp = fopen("points", "w");
        createHeurPointsFile(inst, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotPipe, plot, 1);
        fclose(temp);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }
}

/**
 * @brief 2-Opt swap
 * 
 * @param nodes Route
 * @param start Starting index node for swap
 * @param end Ending index node for swap
 */
void opt_2_swap(int *nodes, int start, int end) {
    int tmp;
    while (start < end) {
        tmp = nodes[start];
        nodes[start] = nodes[end];
        nodes[end] = tmp;
        start++;
        end--;
    }
}

/**
 * @brief Variable Neighborhood Search approach to optimize
 * 
 * @param inst Instance structure
 * @param random_start 
 */
void VNS(instance *inst, int random_start) {
    if (random_start > 1)
        print_error(ERROR_INIT_VALUE_STRING);
    if (inst->timelimit == -1.0) print_error(TIME_LIMIT_NOT_FOUND_STRING);
    //initialization
    if (random_start == 0) {
        nearest_neighborhood(inst);
        if (VERBOSE >= 20) printf(NN_INIT_STRING);
    } else {
        if (VERBOSE >= 20) printf(GRASP_INIT_STRING);
        random_GRASP(inst);
    }

    //local search
    opt_2(inst, inst->nodes_tour, 0);

    int iterations = 0;
    double new_distance, local_best_distance;
    int dim_swapper = 4;
    int swap_size;
    int count_failed_4_opt = 0;
    int count_failed_5_opt = 0;
    int num_iterations = inst->timelimit;
    int dim_continous_failed = num_iterations / 100;
    int *swapper = (int *)malloc(6 * sizeof(int));
    int *local_best_tour = inst->nodes_tour;
    int *new_tour = (int *)calloc(inst->nnodes, sizeof(int));

    local_best_distance = 0;
    for (int i = 0; i < inst->nnodes; i++) {
        local_best_distance += dist(local_best_tour[i], local_best_tour[(i + 1) % inst->nnodes], inst);
    }

    while (iterations <= num_iterations) {
        //  switch to bigger space
        if (count_failed_5_opt > dim_continous_failed)
            dim_swapper = 6;
        else if (count_failed_4_opt > dim_continous_failed)
            dim_swapper = 5;

        printf("\nITERATION NUMBER %d\n", iterations);
        swap_size = 0;
        for (int i = 0; i < dim_swapper; i++) {
            int done = 0;
            while (!done) {
                done = 1;
                swapper[i] = rand() % (inst->nnodes);
                for (int j = 0; j < swap_size; j++) {
                    if (abs(swapper[i] - swapper[j]) < 2 || swapper[i] > (inst->nnodes - 2)) {
                        done = 0;
                        break;
                    }
                }
            }
            swap_size++;
        }

        //sort swapper
        for (int i = 0; i < dim_swapper; i++)
            for (int j = 0; j < dim_swapper; j++)
                if (swapper[j] > swapper[i]) {
                    int tmp = swapper[i];
                    swapper[i] = swapper[j];
                    swapper[j] = tmp;
                }

        // x-opt move
        if (dim_swapper == 4) {
            make_4_opt_move(local_best_tour, new_tour, swapper, inst->nnodes);
        } else if (dim_swapper == 5) {
            make_5_opt_move(local_best_tour, new_tour, swapper, inst->nnodes);
        } else {
            make_6_opt_move(local_best_tour, new_tour, swapper, inst->nnodes);
        }

        opt_2(inst, new_tour, 0);

        new_distance = 0;
        for (int i = 0; i < inst->nnodes; i++) {
            new_distance += dist(new_tour[i], new_tour[(i + 1) % inst->nnodes], inst);
        }

        if (new_distance < local_best_distance) {
            if (dim_swapper == 4)
                count_failed_4_opt = 0;
            else if (dim_swapper == 5)
                count_failed_5_opt = 0;

            local_best_distance = new_distance;
            for (int i = 0; i < inst->nnodes; i++) {
                local_best_tour[i] = new_tour[i];
            }
        } else {
            if (dim_swapper == 4)
                count_failed_4_opt++;
            else if (dim_swapper == 5)
                count_failed_5_opt++;
        }
        iterations++;
    }

    // plot solution
    FILE *gnuplotPipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) {
        gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
        char *commands[] = {"set title \"Variable Neighborhood Search\""};
        sendCommands(gnuplotPipe, commands, 1);
        //temporary file
        temp = fopen("points", "w");
        createHeurPointsFile(inst, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotPipe, plot, 1);
        fclose(temp);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }

    //free memory
    free(swapper);
    free(new_tour);
}

/**
 * @brief 4-Opt move on a tour
 * Takes in input the current optimal tour, the nodes for the move and creates a new one through a 4-opt move
 * 
 * @param local_best_tour The current optimal tour
 * @param new_tour The new tour after the 4-opt move
 * @param swapper Random nodes used for the 4-opt move
 * @param num_nodes Total number of nodes
 */
void make_4_opt_move(int *local_best_tour, int *new_tour, int *swapper, int num_nodes) {
    int size = 0;
    for (int i = 0; i <= swapper[0]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[2] + 1; i <= swapper[3]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[1] + 1; i <= swapper[2]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[0] + 1; i <= swapper[1]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[3] + 1; i < num_nodes; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
}

/**
 * @brief 5-Opt move on a tour
 * Takes in input the current optimal tour, the nodes for the move and creates a new one through a 5-opt move
 * 
 * @param local_best_tour The current optimal tour
 * @param new_tour The new tour after the 5-opt move
 * @param swapper Random nodes used for the 5-opt move
 * @param num_nodes Total number of nodes
 */
void make_5_opt_move(int *local_best_tour, int *new_tour, int *swapper, int num_nodes) {
    int size = 0;
    for (int i = 0; i <= swapper[0]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[2] + 1; i <= swapper[3]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[0] + 1; i <= swapper[1]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[3] + 1; i <= swapper[4]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[1] + 1; i <= swapper[2]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[4] + 1; i < num_nodes; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
}

/**
 * @brief 6-Opt move on a tour
 * Takes in input the current optimal tour, the nodes for the move and creates a new one through a 6-opt move
 * 
 * @param local_best_tour The current optimal tour
 * @param new_tour The new tour after the 6-opt move
 * @param swapper Random nodes used for the 6-opt move
 * @param num_nodes Total number of nodes
 */
void make_6_opt_move(int *local_best_tour, int *new_tour, int *swapper, int num_nodes) {
    int size = 0;
    for (int i = 0; i <= swapper[0]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[1] + 1; i <= swapper[2]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[3] + 1; i <= swapper[4]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[2] + 1; i <= swapper[3]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[4] + 1; i <= swapper[5]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[0] + 1; i <= swapper[1]; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
    for (int i = swapper[5] + 1; i < num_nodes; i++) {
        new_tour[size] = local_best_tour[i];
        size++;
    }
}

/**
 * @brief Simulated Annealing Approach to optimize
 * 
 * @param inst Instance structure
 * @param random_start Initialization flag: 0 for NN, 1 for GRASP 
 */
void simulated_annealing(instance *inst, int random_start) {
    if (inst->timelimit == -1.0) print_error(TIME_LIMIT_NOT_FOUND_STRING);
    if (random_start > 1)
        print_error(ERROR_INIT_VALUE_STRING);
    //initialization
    if (random_start == 0) {
        if (VERBOSE >= 20) printf(NN_INIT_STRING);
        nearest_neighborhood(inst);
    } else {
        if (VERBOSE >= 20) printf(GRASP_INIT_STRING);
        random_GRASP(inst);
    }

    //parameters
    double temperature = 100;                                                   //Initial temperature
    double alpha = 0.999;                                                       //Decreasing factor
    double min_temp = 0.001;                                                    //Stopping criteria
    int partition_count = (int)(log(min_temp / temperature) / log(alpha)) + 1;  //# of operation on temperature
    double inner_time = inst->timelimit / partition_count;                      //time per temperature
    if (VERBOSE >= 20) printf("Time per temperature: %f\n", inner_time);

    //variable for 2-opt operation
    int swap[2];
    double c1, c2, c3, c4, delta;
    struct timespec start, end;
    double diff_time;
    while (temperature > min_temp) {
        diff_time = 0;
        clock_gettime(CLOCK_MONOTONIC, &start);
        while (diff_time < inner_time) {
            //pick 2 different index at random
            int done = 0;
            while (!done) {
                done = 1;
                swap[0] = rand() % (inst->nnodes);
                if (swap[0] == 0)
                    done = 0;
            }
            done = 0;
            while (!done) {
                done = 1;
                swap[1] = rand() % (inst->nnodes);
                if (swap[0] == swap[1] || swap[1] == 0)
                    done = 0;
            }
            //sorting
            if (swap[1] < swap[0]) {
                int tmp = swap[1];
                swap[1] = swap[0];
                swap[0] = tmp;
            }

            //check distance computation
            if (inst->distance[triangle_xpos(swap[0] - 1, swap[0], inst)] == -1)
                inst->distance[triangle_xpos(swap[0], -1, inst)] = dist(inst->nodes_tour[swap[0] - 1], inst->nodes_tour[swap[0]], inst);
            if (inst->distance[triangle_xpos(swap[1], (swap[1] + 1) % (inst->nnodes), inst)] == -1)
                inst->distance[triangle_xpos(swap[1], (swap[1] + 1) % (inst->nnodes), inst)] = dist(inst->nodes_tour[swap[1]], inst->nodes_tour[(swap[1] + 1) % (inst->nnodes)], inst);
            if (inst->distance[triangle_xpos(swap[0] - 1, swap[1], inst)] == -1)
                inst->distance[triangle_xpos(swap[0] - 1, swap[1], inst)] = dist(inst->nodes_tour[swap[0] - 1], inst->nodes_tour[swap[1]], inst);
            if (inst->distance[triangle_xpos(swap[0], (swap[1] + 1) % (inst->nnodes), inst)] == -1)
                inst->distance[triangle_xpos(swap[0] - 1, swap[1], inst)] = dist(inst->nodes_tour[swap[0]], inst->nodes_tour[(swap[1] + 1) % (inst->nnodes)], inst);

            c1 = inst->distance[triangle_xpos(swap[0] - 1, swap[0], inst)];
            c2 = inst->distance[triangle_xpos(swap[1], (swap[1] + 1) % (inst->nnodes), inst)];
            c3 = inst->distance[triangle_xpos(swap[0] - 1, swap[1], inst)];
            c4 = inst->distance[triangle_xpos(swap[0], (swap[1] + 1) % (inst->nnodes), inst)];
            delta = (c3 + c4) - (c1 + c2);  //new edges - old edges cost
            if (c1 == -1 || c2 == -1 || c3 == -1 || c4 == -1)
                print_error(ERROR_DISTANCE_STRING);

            if (delta < 0)
                opt_2_swap(inst->nodes_tour, swap[0], swap[1]);  //cost is decreasing
            else {
                if ((double)rand() / (double)RAND_MAX < exp(-delta / temperature)) {
                    opt_2_swap(inst->nodes_tour, swap[0], swap[1]);  //cost is increasing
                }
            }
            clock_gettime(CLOCK_MONOTONIC, &end);
            diff_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
        }
        //update temperature
        temperature = temperature * alpha;
    }
    //compute cost afer SA
    if (VERBOSE >= 50) {
        double SA_sol = 0;
        for (int i = 0; i < inst->nnodes; i++) {
            SA_sol += dist(inst->nodes_tour[i], inst->nodes_tour[(i + 1) % (inst->nnodes)], inst);
        }
        printf("After exec SA: %f\n", SA_sol);
    }
    //optimize solution
    opt_2(inst, inst->nodes_tour, 0);

    // plot solution
    FILE *gnuplotPipe = NULL;
    FILE *temp = NULL;
    if (VERBOSE >= VERBOSE_TO_PLOT) {
        gnuplotPipe = popen(GNUPLOT_OPEN_COMMAND, "w");
        char *commands[] = {"set title \"Simulated Annealing\""};
        sendCommands(gnuplotPipe, commands, 1);
        //temporary file
        temp = fopen("points", "w");
        createHeurPointsFile(inst, temp);
        char *plot[] = {"plot 'points' with lines"};
        sendCommands(gnuplotPipe, plot, 1);
        fclose(temp);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }
}

/**
 * @brief build CPLEX solution array
 * 
 * @param inst Pointer to instance struct
 * @param x_value Array to store CPLEX solution
 */
void build_cplex_sol(instance *inst, double *x_value) {
    for (int i = 0; i < inst->nnodes; i++) {
        x_value[triangle_xpos(inst->nodes_tour[i], inst->nodes_tour[(i + 1) % (inst->nnodes)], inst)] = 1.0;
        // printf("%d %d\n", inst->nodes_tour[i], inst->nodes_tour[(i + 1) % (inst->nnodes)]);
    }
}