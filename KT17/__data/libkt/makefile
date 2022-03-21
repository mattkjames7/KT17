# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC -std=c++17 -g

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib

obj:
	cd magnetopause; make all
	cd model; make all
	cd spline; make all
	cd trace; make all
	$(CCo) libkt.cc -o libkt.o

lib:
	$(CC) -shared -o libkt.so *.o \
			magnetopause/*.o \
			model/*.o \
			spline/*.o \
			trace/*.o
	
	rm -v *.o
	rm -v magnetopause/*.o
	rm -v model/*.o
	rm -v spline/*.o
	rm -v trace/*.o

clean:
	rm -v *.so

