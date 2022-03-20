# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC -std=c++17 -g

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj 

obj:
	$(CCo) dipole.cc -o dipole.o
	$(CCo) disk.cc -o disk.o
	$(CCo) qhsheet.cc -o qhsheet.o
	$(CCo) shield.cc -o shield.o
	$(CCo) kt14.cc -o kt14.o
	$(CCo) kt17.cc -o kt17.o
	$(CCo) ktmodel.cc -o ktmodel.o
	$(CCo) model.cc -o model.o


clean:
	rm -v *.o

