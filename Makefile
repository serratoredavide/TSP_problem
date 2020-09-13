OBJS = main.o tsp.o utils.o compact_models_TSP.o heuristic_models_TSP.o heuristic_solver_TSP.o
HEADERS = tsp.h utils.h mystring.h compact_models_TSP.h heuristic_models_TSP.h heuristic_solver_TSP.h
HEADERS_DIR = /include/ 
EXE = exec_tsp
all: $(EXE) 
setting = -1   
OS := $(shell uname)

# ---------------------------------------------------------------------
# Compiler selection and libraries for Cplex -- windows
# ---------------------------------------------------------------------
#
#CC = c:/usr/mingw/bin/gcc -mno-cygwin
#LIBS = -Lc:/usr/cplex -lcplex
#INC = -Ic:/usr/cplex
#
# ---------------------------------------------------------------------
# Compiler selection and libraries for Cplex 
# ---------------------------------------------------------------------
# 

ifeq ($(OS),Darwin)
	setting = 0
	CPLEX_HOME = /Applications/CPLEX_Studio1210/cplex
	CC = clang -Qunused-arguments
	AR = ar rc
	LIBS = -Wl,-no_compact_unwind -L${CPLEX_HOME}/lib/x86-64_osx/static_pic -L. -lcplex -lm -lpthread 
	INC = -I. -I${CPLEX_HOME}/include/ilcplex -I./$(HEADERS_DIR) 
endif

ifeq ($(OS),Linux)
	setting = 1
	CPLEX_HOME = /opt/ibm/ILOG/CPLEX_Studio1210/cplex
	CC = gcc
	AR = ar rc
	LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
	INC = -I. -I${CPLEX_HOME}/include/ilcplex -I./$(HEADERS_DIR)
endif


# ---------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------
CFLAGS = -Wall -O3
RM = rm -rf

VPATH = ./include:./src

.SUFFIXES: .o .c
.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(EXE): $(OBJS) $(LIBUTILS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(OBJS) : $(HEADERS)

clean:
	$(RM) $(OBJS)
	$(RM) $(EXE) 
	$(RM) points
	$(RM) model.lp
	
again:                                                               
	make clean
	make    

who:
	@echo "you are user $(USER) with uname `uname` (OS = $(OS)) and you working with compiler setting $(setting)" 

