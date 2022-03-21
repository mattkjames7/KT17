# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC -std=c++17 -g

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj 

obj:
	$(CCo) conttrace.cc -o conttrace.o
	$(CCo) interptraceclosestpos.cc -o interptraceclosestpos.o
	$(CCo) latlt.cc -o latlt.o
	$(CCo) fieldlinernorm.cc -o fieldlinernorm.o
	$(CCo) trace.cc -o trace.o



clean:
	rm -v *.o

