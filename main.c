#include "compact_models_TSP.h"
#include "heuristic_models_TSP.h"
#include "heuristic_solver_TSP.h"
#include "mystring.h"
#include "tsp.h"
#include "utils.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s -help for help\n", argv[0]);
        return 0;
    }
    //struct initialization
    instance inst;

    //read commands and input
    parse_command_line(argc, argv, &inst);
    read_input(&inst);
    printf("Numbers of node: %d\n", inst.nnodes);

    // init_heuristic_inst(&inst);
    // VNS(&inst);
    // simulated_annealing(&inst);
    //nearest_neighborhood(&inst);
    // random_GRASP(&inst);
    // return 0;

    //solve problem
    int opt_code = 0;
    switch (inst.model_type) {
        case TSP_LOOP:
            printf("Solving problem using loop\n");
            opt_code = TSPopt(&inst);
            break;
        case TSP_LAZY_CALLBACK:
            printf("Solving problem using lazy callback\n");
            opt_code = TSPopt_callback(&inst);
            break;
        case COMPACT_TSP_MTZ:
            printf("Solving problem using Compact MTZ model\n");
            opt_code = compact_TSPopt(&inst, 0, 0);
            break;
        case COMPACT_TSP_LAZY_MTZ:
            printf("Solving problem using Compact MTZ_lazy model\n");
            opt_code = compact_TSPopt(&inst, 0, 1);
            break;
        case COMPACT_TSP_FLOW1:
            printf("Solving problem using Compact FLOW1 model\n");
            opt_code = compact_TSPopt(&inst, 1, 0);
            break;
        case COMPACT_TSP_LAZY_FLOW1:
            printf("Solving problem using Compact FLOW1_lazy model\n");
            opt_code = compact_TSPopt(&inst, 1, 1);
            break;
        case HARD_FIXING_TSP:
            printf("Solving problem using Heuristic Hard-Fixing\n");
            opt_code = hard_fixing_TSP(&inst);
            break;
        case SOFT_FIXING_TSP:
            printf("Solving problem using Heuristic Soft-Fixing\n");
            opt_code = soft_fixing_TSP(&inst);
            break;
        case OPT_2_TSP:
            init_heuristic_inst(&inst);
            printf("Solving problem using 2-Opt\n");
            //init
            if (inst.initialization == 0)
                nearest_neighborhood(&inst);
            else
                random_GRASP(&inst);
            opt_2(&inst, inst.nodes_tour, 1);
            break;
        case VNS_TSP:
            init_heuristic_inst(&inst);
            printf("Solving problem using Variable Neighborhood Search\n");
            VNS(&inst, inst.initialization);
            break;
        case SA_TSP:
            init_heuristic_inst(&inst);
            printf("Solving problem using Simulated Annealing\n");
            simulated_annealing(&inst, inst.initialization);
            break;
        default:
            printf("Model Type is not define correctly");
            return 0;
            break;
    }

    if (opt_code) printf("Solution problem.\nOpt_code: %d\n", opt_code);

    //free instance struct
    free_instance(&inst);

    return 0;
}
