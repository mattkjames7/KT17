# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC -std=c++17 -g

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj 

obj:
	$(CCo) magnetopause.cc -o magnetopause.o
	$(CCo) withinmp.cc -o withinmp.o

clean:
	rm -v *.o

