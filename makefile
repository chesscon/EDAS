
CCOMP=gcc
C_FLAGS=-lm -lgsl -lgslcblas

VFLAGS=--leak-check=full --show-leak-kinds=all

SRC_DIR=src

SRC_Main=testEDAs.c EDAS.c individuos.c funciones_objetivo.c

testEDAs: $(addprefix $(SRC_DIR)/, $(SRC_Main) )
	$(CCOMP) -O2 -o $@ $(addprefix $(SRC_DIR)/, $(SRC_Main) ) $(C_FLAGS) 

valgrindT1: $(addprefix $(SRC_DIR)/, $(SRC_Main) )
	$(CCOMP) -g -o testEDAs $^ $(C_FLAGS) 	
	valgrind $(V_FLAGS) ./testEDAs 1 3 1000 -600 600 10 0.5 1
	rm -f testEDAs

VAL_POR_SELECC=0.4
VAL_NUM_REP=30
VAL_TAM_POB=200
VAL_DIM=10
VAL_MAX_EVAL=200000

# Test Function 1: SPHERE
runtest_1: testEDAs
	./testEDAs 1 $(VAL_DIM) $(VAL_MAX_EVAL) -600 600 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

# Test Function 2: ELLIPSOID
runtest_2: testEDAs
	./testEDAs 2 $(VAL_DIM) $(VAL_MAX_EVAL) -20 20 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

# Test Function 3: Zakharov
runtest_3: testEDAs
	./testEDAs 3 $(VAL_DIM) $(VAL_MAX_EVAL) -20 20 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

# Test Function 4: Rosenbrock
runtest_4: testEDAs
	./testEDAs 4 $(VAL_DIM) $(VAL_MAX_EVAL) -20 20 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

# Test Function 5: Ackley
runtest_5: testEDAs
	./testEDAs 5 $(VAL_DIM) $(VAL_MAX_EVAL) -20 20 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

# Test Function 6: Griewangk
runtest_6: testEDAs
	./testEDAs 6 $(VAL_DIM) $(VAL_MAX_EVAL) -600 600 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

# Test Function 7: Rastrigin
runtest_7: testEDAs
	./testEDAs 7 $(VAL_DIM) $(VAL_MAX_EVAL) -20 20 $(VAL_TAM_POB) $(VAL_POR_SELECC) $(VAL_NUM_REP)

clean:
	rm -f testEDAs