CC=g++
CFLAGS=-Wall -Wextra -O3 -ftree-vectorize -msse3
LDFLAGS=
SOURCES=linear_solver/blas/blas.c linear_solver/conjugate_solver.c working_set.cpp eq_con_lcqp.cpp lcqp.cpp quadprog/src/QuadProg++.cc main.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -o $@ -lm

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) *.o linear_solver/*.o linear_solver/blas/*.o