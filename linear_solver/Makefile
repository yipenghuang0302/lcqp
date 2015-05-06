CC=gcc
CFLAGS=-Wall -Wextra -Werror -O3 -ftree-vectorize -msse3
LDFLAGS=
SOURCES=main.c conjugate_solver.c steepest_solver.c blas/blas.c
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -o $@ -lm

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) *.o blas/*.o
